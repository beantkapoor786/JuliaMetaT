const BASE_A = Int8(1)
const BASE_T = Int8(2)
const BASE_G = Int8(3)
const BASE_C = Int8(4)
const BASE_N = Int8(0)

@inline function encode_base(c::UInt8)::Int8
    c == UInt8('A') || c == UInt8('a') ? BASE_A :
    c == UInt8('T') || c == UInt8('t') ? BASE_T :
    c == UInt8('G') || c == UInt8('g') ? BASE_G :
    c == UInt8('C') || c == UInt8('c') ? BASE_C :
                                         BASE_N
end

@inline function decode_base(b::Int8)::Char
    b == BASE_A ? 'A' :
    b == BASE_T ? 'T' :
    b == BASE_G ? 'G' :
    b == BASE_C ? 'C' :
                  'N'
end

function encode_reads(records::AbstractVector{<:FASTX.FASTQ.Record};
                      drop_with_n::Bool = true)
    n = length(records)
    n == 0 && return (zeros(Int8, 0, 0), Int32[], 0)

    keep = trues(n)
    max_len = 0
    @inbounds for (j, r) in enumerate(records)
        s = FASTX.sequence(r)
        L = length(s)
        sb = codeunits(s)
        contains_n = false
        for i in 1:L
            if encode_base(sb[i]) == BASE_N
                contains_n = true
                break
            end
        end
        if drop_with_n && contains_n
            keep[j] = false
        else
            L > max_len && (max_len = L)
        end
    end

    n_kept = count(keep)
    n_dropped = n - n_kept
    n_kept == 0 && return (zeros(Int8, 0, 0), Int32[], n_dropped)

    seqs    = zeros(Int8, max_len, n_kept)
    lengths = zeros(Int32, n_kept)

    col = 0
    @inbounds for j in 1:n
        keep[j] || continue
        col += 1
        s = FASTX.sequence(records[j])
        sb = codeunits(s)
        L = length(s)
        lengths[col] = L
        for i in 1:L
            seqs[i, col] = encode_base(sb[i])
        end
    end

    return seqs, lengths, n_dropped
end

function decode_read(seqs::AbstractMatrix{Int8},
                     lengths::AbstractVector{<:Integer},
                     j::Integer)
    L = Int(lengths[j])
    chars = Vector{Char}(undef, L)
    @inbounds for i in 1:L
        chars[i] = decode_base(seqs[i, j])
    end
    return String(chars)
end

"""
    load_fastq(paths...; chunk_size=1_000_000, drop_with_n=true) -> (seqs, lengths, n_dropped)

Read one or more FASTQ files (e.g., R1 and R2 for paired-end) and
encode all reads into a single combined matrix. Files ending in `.gz`
are transparently decompressed. Streams chunks internally to avoid
holding all FASTX records in memory simultaneously; downstream API is
unchanged from earlier versions.
"""
function load_fastq(paths::AbstractString...;
                    chunk_size::Int = 1_000_000,
                    drop_with_n::Bool = true)
    chunks = NamedTuple[]
    total_dropped = 0
    max_len = 0
    total_reads = 0

    for chunk in stream_fastq(paths...; chunk_size, drop_with_n)
        push!(chunks, chunk)
        total_dropped += chunk.n_dropped
        total_reads += length(chunk.lengths)
        max_len = max(max_len, size(chunk.seqs, 1))
    end

    if total_reads == 0
        return zeros(Int8, 0, 0), Int32[], total_dropped
    end

    # Allocate the final matrix once at known size; copy chunks in.
    seqs    = zeros(Int8, max_len, total_reads)
    lengths = Vector{Int32}(undef, total_reads)
    offset  = 0
    for chunk in chunks
        n_chunk = length(chunk.lengths)
        chunk_max_len = size(chunk.seqs, 1)
        @views seqs[1:chunk_max_len, offset+1:offset+n_chunk] .= chunk.seqs
        lengths[offset+1:offset+n_chunk] .= chunk.lengths
        offset += n_chunk
    end

    return seqs, lengths, total_dropped
end

"""
    stream_fastq(paths...; chunk_size=1_000_000, drop_with_n=true) -> Channel

Stream FASTQ files in chunks. Each yielded value is a NamedTuple
`(seqs, lengths, n_dropped)` where `seqs` is an `Int8` matrix of shape
`(max_len_in_chunk × n_reads_in_chunk)`, `lengths` is a vector of read
lengths, and `n_dropped` is the number of N-containing reads dropped
from this chunk.

Files ending in `.gz` are transparently decompressed. Streams are
closed when the channel closes (normally or via exception).
"""
function stream_fastq(paths::AbstractString...;
                      chunk_size::Int = 1_000_000,
                      drop_with_n::Bool = true)
    return Channel{NamedTuple}() do ch
        readers_and_streams = Tuple{FASTX.FASTQ.Reader,Any}[]
        try
            for p in paths
                if endswith(lowercase(p), ".gz")
                    stream = CodecZlib.GzipDecompressorStream(open(p))
                    push!(readers_and_streams, (FASTX.FASTQ.Reader(stream), stream))
                else
                    stream = open(p)
                    push!(readers_and_streams, (FASTX.FASTQ.Reader(stream), stream))
                end
            end

            buffer = FASTX.FASTQ.Record[]
            sizehint!(buffer, chunk_size)

            # Round-robin reading across paired files so chunks contain
            # interleaved mates (matches existing load_fastq semantics:
            # pair info isn't retained beyond this layer).
            done = false
            while !done
                done = true
                for (reader, _) in readers_and_streams
                    if !eof(reader)
                        done = false
                        rec = FASTX.FASTQ.Record()
                        read!(reader, rec)
                        push!(buffer, rec)
                        if length(buffer) >= chunk_size
                            seqs, lengths, n_dropped = encode_reads(buffer; drop_with_n)
                            put!(ch, (; seqs, lengths, n_dropped))
                            empty!(buffer)
                        end
                    end
                end
            end

            # Flush the final partial chunk.
            if !isempty(buffer)
                seqs, lengths, n_dropped = encode_reads(buffer; drop_with_n)
                put!(ch, (; seqs, lengths, n_dropped))
                empty!(buffer)
            end
        finally
            for (_, stream) in readers_and_streams
                close(stream)
            end
        end
    end
end
