# Read mapping by k-mer pseudo-alignment + abundance estimation.

struct ContigAbundance
    contig_id::Int          # index into the contigs vector
    contig_length::Int      # bp
    raw_count::Int          # number of mapped reads
    rpkm::Float64
    tpm::Float64
    mean_coverage::Float32  # from contig assembly
end

# Build a sorted lookup table mapping every canonical k-mer that appears
# in any contig to its contig index. K-mers shared across contigs (rare
# for distinct transcripts) get assigned to the first contig that owns them.
"""
    build_contig_kmer_index(contigs, k) -> (sorted_kmers, contig_of_kmer)

Returns two parallel host arrays:
- sorted_kmers: Vector{UInt64} of canonical k-mers, sorted ascending
- contig_of_kmer: Vector{Int32} same length, giving the contig index
  each k-mer belongs to.
"""
function build_contig_kmer_index(contigs::AbstractVector{Contig}, k::Int)
    pairs = Tuple{UInt64,Int32}[]
    for (ci, c) in enumerate(contigs)
        s = c.sequence
        length(s) >= k || continue
        for i in 1:(length(s) - k + 1)
            km = string_to_canonical_kmer(s[i:i+k-1], k)
            push!(pairs, (km, Int32(ci)))
        end
    end
    # Sort by k-mer; dedup by keeping the first contig assignment for
    # each unique k-mer.
    sort!(pairs; by = p -> p[1])

    n = length(pairs)
    if n == 0
        return UInt64[], Int32[]
    end

    sorted_kmers = UInt64[]
    contig_of_kmer = Int32[]
    sizehint!(sorted_kmers, n)
    sizehint!(contig_of_kmer, n)
    prev = typemax(UInt64)
    @inbounds for (km, ci) in pairs
        if km != prev
            push!(sorted_kmers, km)
            push!(contig_of_kmer, ci)
            prev = km
        end
    end
    return sorted_kmers, contig_of_kmer
end

# --- GPU kernel: one thread per read ---
#
# Each read scans its k-mers, binary-searches the contig index, and
# tallies hits per contig. Then it picks the contig with the most hits
# (majority vote). Output: out_assignment[j] = contig index (0 = unmapped).
#
# We use a per-thread small array of counts up to n_contigs. To avoid
# allocating large per-thread buffers, we cap n_contigs at COUNT_CAP and
# linear-scan for the max. For an MVP with ~3-30 contigs this is fine.

const COUNT_CAP = Int32(64)   # max contigs supported by this kernel

function _map_kernel!(j, seqs, lengths, sorted_kmers, contig_of_kmer,
                     out_assignment, out_hits,
                     k::Int32, n_contigs::Int32, min_hits::Int32)
    L = lengths[j]
    if L < k
        out_assignment[j] = Int32(0)
        out_hits[j]       = Int32(0)
        return nothing
    end

    # Per-thread tally array (stack-allocated via mutable local).
    # We cap at COUNT_CAP — kernel users must ensure n_contigs <= COUNT_CAP.
    counts = MVector{Int32, Int(COUNT_CAP)}(zeros(Int32, Int(COUNT_CAP)))
    # ... well, MVector requires StaticArrays. Fall back to a simple
    # approach: pre-allocate per-thread arrays on host and pass them in.
    # See the higher-level wrapper.

    # NOTE: this kernel signature accepts pre-allocated `counts_buf` of
    # shape (COUNT_CAP, n_reads). See actual signature below.
    return nothing
end

# Real kernel (with pre-allocated per-thread counts buffer):
function _map_kernel_v2!(j, seqs, lengths, sorted_kmers, contig_of_kmer,
                          counts_buf, out_assignment, out_hits,
                          n_kmers_in_index::Int32,
                          k::Int32, n_contigs::Int32, min_hits::Int32)
    L = lengths[j]

    # Zero the count column for thread j.
    @inbounds for c in Int32(1):n_contigs
        counts_buf[c, j] = Int32(0)
    end

    if L < k
        out_assignment[j] = Int32(0)
        out_hits[j]       = Int32(0)
        return nothing
    end

    # Compute and look up each k-mer in this read.
    @inbounds for p in Int32(1):(L - k + Int32(1))
        # Build packed k-mer + its rev-comp; take canonical.
        kmer = UInt64(0)
        valid = true
        for i in Int32(0):(k - Int32(1))
            b = seqs[p + i, j]
            # remap A=1,C=4,G=3,T=2 -> 0..3
            bb = b == Int8(1) ? UInt64(0) :
                 b == Int8(4) ? UInt64(1) :
                 b == Int8(3) ? UInt64(2) :
                 b == Int8(2) ? UInt64(3) :
                                UInt64(99)   # N or padding
            if bb == UInt64(99)
                valid = false
                break
            end
            kmer = (kmer << 2) | bb
        end
        valid || continue

        # Canonicalize.
        mask = (UInt64(1) << (UInt64(2) * UInt64(k))) - UInt64(1)
        x = (~kmer) & mask
        rc = UInt64(0)
        for i in Int32(0):(k - Int32(1))
            pair = (x >> (UInt64(2) * UInt64(i))) & UInt64(0x3)
            rc |= pair << (UInt64(2) * UInt64(k - Int32(1) - i))
        end
        canon = kmer < rc ? kmer : rc

        # Binary search in sorted_kmers.
        lo = Int32(1)
        hi = n_kmers_in_index
        found_at = Int32(0)
        while lo <= hi
            mid = (lo + hi) >> 1
            v = sorted_kmers[mid]
            if v == canon
                found_at = mid; break
            elseif v < canon
                lo = mid + Int32(1)
            else
                hi = mid - Int32(1)
            end
        end
        if found_at > 0
            ci = contig_of_kmer[found_at]
            counts_buf[ci, j] += Int32(1)
        end
    end

    # Majority vote.
    best_c = Int32(0)
    best_n = Int32(0)
    @inbounds for c in Int32(1):n_contigs
        v = counts_buf[c, j]
        if v > best_n
            best_n = v
            best_c = c
        end
    end

    out_assignment[j] = best_n >= min_hits ? best_c : Int32(0)
    out_hits[j]       = best_n
    return nothing
end

"""
    map_reads(contigs, seqs_h, lengths_h; k=31, min_hits=3) -> (assignments, hits)

Map each read to its best-matching contig by canonical k-mer pseudo-
alignment. Reads with fewer than `min_hits` matches go to contig 0
(unmapped). Returns host arrays length n_reads.
"""
function map_reads(contigs::AbstractVector{Contig},
                   seqs_h::AbstractMatrix{Int8},
                   lengths_h::AbstractVector{<:Integer};
                   k::Int = 31,
                   min_hits::Int = 3)
    n_reads = length(lengths_h)
    n_contigs = length(contigs)
    @assert n_contigs <= COUNT_CAP "more than $COUNT_CAP contigs not supported in this kernel"
    n_reads == 0 && return (Int32[], Int32[])
    n_contigs == 0 && return (zeros(Int32, n_reads), zeros(Int32, n_reads))

    sorted_kmers, contig_of_kmer = build_contig_kmer_index(contigs, k)
    n_idx = length(sorted_kmers)

    if n_idx == 0
        return (zeros(Int32, n_reads), zeros(Int32, n_reads))
    end

    # Move to device.
    seqs_d            = JACC.to_device(seqs_h)
    lengths_d         = JACC.to_device(Int32.(lengths_h))
    sorted_kmers_d    = JACC.to_device(sorted_kmers)
    contig_of_kmer_d  = JACC.to_device(contig_of_kmer)
    counts_buf_d      = JACC.to_device(zeros(Int32, Int(COUNT_CAP), n_reads))
    out_assignment_d  = JACC.to_device(zeros(Int32, n_reads))
    out_hits_d        = JACC.to_device(zeros(Int32, n_reads))

    JACC.parallel_for(n_reads, _map_kernel_v2!,
        seqs_d, lengths_d, sorted_kmers_d, contig_of_kmer_d,
        counts_buf_d, out_assignment_d, out_hits_d,
        Int32(n_idx), Int32(k), Int32(n_contigs), Int32(min_hits))

    return JACC.to_host(out_assignment_d), JACC.to_host(out_hits_d)
end

"""
    compute_abundance(contigs, assignments) -> Vector{ContigAbundance}

Aggregate per-contig raw counts, RPKM, and TPM from read assignments.
"""
function compute_abundance(contigs::AbstractVector{Contig},
                           assignments::AbstractVector{<:Integer})
    n_contigs = length(contigs)
    raw_counts = zeros(Int, n_contigs)
    for a in assignments
        a == 0 && continue
        raw_counts[a] += 1
    end
    total_mapped = sum(raw_counts)
    total_reads_M = max(total_mapped, 1) / 1e6   # avoid div by zero

    # RPKM = count / (length_kb * total_reads_M)
    # TPM normalizes per-length rates to sum to 1e6
    per_kb_rate = [raw_counts[c] / (length(contigs[c].sequence) / 1000.0)
                   for c in 1:n_contigs]
    sum_rate = sum(per_kb_rate)
    sum_rate == 0 && (sum_rate = 1.0)   # avoid div by zero

    out = ContigAbundance[]
    for c in 1:n_contigs
        len = length(contigs[c].sequence)
        rpkm = raw_counts[c] / (len / 1000.0) / total_reads_M
        tpm  = (per_kb_rate[c] / sum_rate) * 1e6
        push!(out, ContigAbundance(c, len, raw_counts[c],
                                   rpkm, tpm,
                                   contigs[c].mean_coverage))
    end
    return out
end
