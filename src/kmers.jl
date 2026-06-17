# Remap storage encoding (A=1,C=4,G=3,T=2) to 2-bit (A=0,C=1,G=2,T=3).
@inline function _pack2(b::Int8)::UInt64
    b == Int8(1) ? UInt64(0) :   # A
    b == Int8(4) ? UInt64(1) :   # C
    b == Int8(3) ? UInt64(2) :   # G
                   UInt64(3)     # T (b == 2)
end

# Reverse complement of a packed k-mer of length k.
@inline function _revcomp(kmer::UInt64, k::Int32)::UInt64
    mask = (UInt64(1) << (UInt64(2) * UInt64(k))) - UInt64(1)
    x = (~kmer) & mask
    r = UInt64(0)
    @inbounds for i in Int32(0):(k - Int32(1))
        pair = (x >> (UInt64(2) * UInt64(i))) & UInt64(0x3)
        r |= pair << (UInt64(2) * UInt64(k - Int32(1) - i))
    end
    return r
end

@inline function _canonical(kmer::UInt64, k::Int32)::UInt64
    rc = _revcomp(kmer, k)
    return kmer < rc ? kmer : rc
end

# Kernel: one thread per output k-mer slot g.
function _kmer_kernel!(g, seqs, lengths, read_of, pos_of, out_kmers, k::Int32)
    j = read_of[g]
    p = pos_of[g]
    L = lengths[j]
    if p < Int32(1) || p + k - Int32(1) > L
        out_kmers[g] = typemax(UInt64)
        return nothing
    end
    kmer = UInt64(0)
    @inbounds for i in Int32(0):(k - Int32(1))
        b = seqs[p + i, j]
        kmer = (kmer << 2) | _pack2(b)
    end
    out_kmers[g] = _canonical(kmer, k)
    return nothing
end

"""
    extract_kmers(seqs_h, lengths_h; k=31) -> Vector{UInt64}

Run the GPU kernel to produce a flat vector of canonical packed
k-mers, one per (read, valid-position) slot. Returns a host vector.
"""
function extract_kmers(seqs_h::AbstractMatrix{Int8},
                       lengths_h::AbstractVector{<:Integer};
                       k::Int = 31)
    n_reads = length(lengths_h)
    n_reads == 0 && return UInt64[]

    per_read = [max(0, Int(L) - k + 1) for L in lengths_h]
    total = sum(per_read)
    total == 0 && return UInt64[]

    read_of = Vector{Int32}(undef, total)
    pos_of  = Vector{Int32}(undef, total)
    g = 0
    @inbounds for j in 1:n_reads
        for p in 1:per_read[j]
            g += 1
            read_of[g] = Int32(j)
            pos_of[g]  = Int32(p)
        end
    end

    seqs_d    = JACC.to_device(seqs_h)
    lengths_d = JACC.to_device(Int32.(lengths_h))
    read_of_d = JACC.to_device(read_of)
    pos_of_d  = JACC.to_device(pos_of)
    kmers_d   = JACC.to_device(zeros(UInt64, total))

    JACC.parallel_for(total, _kmer_kernel!,
        seqs_d, lengths_d, read_of_d, pos_of_d, kmers_d, Int32(k))

    return JACC.to_host(kmers_d)
end

"""
    count_kmers(kmers; min_count=2) -> (unique_kmers, counts)

Sort the flat k-mer array and run-length encode. Drops k-mers with
count < min_count.
"""
function count_kmers(kmers::Vector{UInt64}; min_count::Int = 2)
    isempty(kmers) && return (UInt64[], Int32[])
    sort!(kmers)

    last_real = searchsortedlast(kmers, typemax(UInt64) - UInt64(1))
    last_real == 0 && return (UInt64[], Int32[])
    kmers_view = view(kmers, 1:last_real)

    unique_kmers = UInt64[]
    counts       = Int32[]
    sizehint!(unique_kmers, length(kmers_view) ÷ 4)
    sizehint!(counts, length(kmers_view) ÷ 4)

    i = 1
    n = length(kmers_view)
    @inbounds while i <= n
        cur = kmers_view[i]
        j = i
        while j <= n && kmers_view[j] == cur
            j += 1
        end
        c = j - i
        if c >= min_count
            push!(unique_kmers, cur)
            push!(counts, Int32(c))
        end
        i = j
    end
    return unique_kmers, counts
end

# Test/debug helpers.

function kmer_to_string(kmer::UInt64, k::Int)
    bytes = Vector{UInt8}(undef, k)
    @inbounds for i in 0:(k-1)
        bits = (kmer >> (2 * (k - 1 - i))) & UInt64(0x3)
        bytes[i+1] = bits == 0x0 ? UInt8('A') :
                     bits == 0x1 ? UInt8('C') :
                     bits == 0x2 ? UInt8('G') :
                                   UInt8('T')
    end
    return String(bytes)
end

function string_to_canonical_kmer(s::AbstractString, k::Int)
    @assert length(s) == k
    kmer = UInt64(0)
    for c in s
        bits = c == 'A' || c == 'a' ? UInt64(0) :
               c == 'C' || c == 'c' ? UInt64(1) :
               c == 'G' || c == 'g' ? UInt64(2) :
               c == 'T' || c == 't' ? UInt64(3) :
                                       error("Invalid base: $c")
        kmer = (kmer << 2) | bits
    end
    return _canonical(kmer, Int32(k))
end
