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
    load_fastq(paths...; drop_with_n=true) -> (seqs, lengths, n_dropped)

Read one or more FASTQ files (e.g., R1 and R2 for paired-end) and
encode all reads into a single combined matrix. Files ending in `.gz`
are transparently decompressed. Pair info is not retained at this layer.
"""
function load_fastq(paths::AbstractString...; drop_with_n::Bool = true)
    records = FASTX.FASTQ.Record[]
    for p in paths
        if endswith(lowercase(p), ".gz")
            stream = CodecZlib.GzipDecompressorStream(open(p))
            try
                reader = FASTX.FASTQ.Reader(stream)
                for rec in reader
                    push!(records, rec)
                end
            finally
                close(stream)
            end
        else
            open(FASTX.FASTQ.Reader, p) do r
                for rec in r
                    push!(records, rec)
                end
            end
        end
    end
    return encode_reads(records; drop_with_n)
end
