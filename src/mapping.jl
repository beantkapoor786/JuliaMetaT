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
@inline function _char_to_2bit(c::Char)::UInt64
    c == 'A' || c == 'a' ? UInt64(0) :
    c == 'C' || c == 'c' ? UInt64(1) :
    c == 'G' || c == 'g' ? UInt64(2) :
                           UInt64(3)
end

function build_contig_kmer_index(contigs::AbstractVector{Contig}, k::Int)
    k32 = Int32(k)
    mask = (UInt64(1) << (UInt64(2) * UInt64(k))) - UInt64(1)
    total_kmers = sum(max(0, length(c.sequence) - k + 1) for c in contigs; init=0)
    pairs = Vector{Tuple{UInt64,Int32}}(undef, total_kmers)
    pi = 0
    for (ci, c) in enumerate(contigs)
        s = c.sequence
        n = length(s)
        n >= k || continue
        # Build first k-mer
        kmer = UInt64(0)
        @inbounds for i in 1:k
            kmer = (kmer << 2) | _char_to_2bit(s[i])
        end
        rc = _revcomp_swar(kmer, k32)
        pairs[pi += 1] = (kmer < rc ? kmer : rc, Int32(ci))
        # Roll subsequent k-mers: drop leftmost base, shift in new base
        @inbounds for i in (k+1):n
            kmer = ((kmer << 2) | _char_to_2bit(s[i])) & mask
            rc = _revcomp_swar(kmer, k32)
            pairs[pi += 1] = (kmer < rc ? kmer : rc, Int32(ci))
        end
    end
    resize!(pairs, pi)
    # Sort by k-mer (lex order on Tuple puts UInt64 first); dedup keeping
    # the first contig assignment for each unique k-mer.
    sort!(pairs)

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

# SWAR 5-step 2-bit-group reversal: reverses order of k 2-bit groups in low 2k bits.
@inline function _revcomp_swar(kmer::UInt64, k::Int32)::UInt64
    mask = (UInt64(1) << (UInt64(2) * UInt64(k))) - UInt64(1)
    x = (~kmer) & mask
    # Reverse order of 2-bit groups via 5-step word-level shuffle.
    x = ((x & UInt64(0x3333333333333333)) << 2)  | ((x >> 2)  & UInt64(0x3333333333333333))
    x = ((x & UInt64(0x0F0F0F0F0F0F0F0F)) << 4)  | ((x >> 4)  & UInt64(0x0F0F0F0F0F0F0F0F))
    x = ((x & UInt64(0x00FF00FF00FF00FF)) << 8)  | ((x >> 8)  & UInt64(0x00FF00FF00FF00FF))
    x = ((x & UInt64(0x0000FFFF0000FFFF)) << 16) | ((x >> 16) & UInt64(0x0000FFFF0000FFFF))
    x = (x << 32) | (x >> 32)
    x >> (UInt64(64) - UInt64(2) * UInt64(k))
end

# Streaming majority-vote kernel with rolling hash and SWAR revcomp.
# Pass 1: Boyer-Moore majority candidate via O(L) rolling scan.
# Pass 2: verify candidate hit count >= min_hits.
function _map_kernel_streaming!(j, seqs, lengths, sorted_kmers, contig_of_kmer,
                                out_assignment, out_hits,
                                n_kmers_in_index::Int32,
                                k::Int32, min_hits::Int32)
    L = lengths[j]
    if L < k
        out_assignment[j] = Int32(0)
        out_hits[j]       = Int32(0)
        return nothing
    end

    mask = (UInt64(1) << (UInt64(2) * UInt64(k))) - UInt64(1)

    # --- Pass 1: Boyer-Moore majority candidate (rolling hash) ---
    # Also tracks candidate_hits: raw count of k-mers matching the *current*
    # candidate (resets when BM evicts the candidate). If candidate_hits >=
    # min_hits at the end, Pass 2 is skipped — any read that maps cleanly
    # without a mid-read BM candidate switch already has its answer here.
    candidate       = Int32(0)
    counter         = Int32(0)
    candidate_hits  = Int32(0)
    kmer            = UInt64(0)
    run             = Int32(0)

    @inbounds for pos in Int32(1):L
        b  = seqs[pos, j]
        bb = b == Int8(1) ? UInt64(0) :
             b == Int8(4) ? UInt64(1) :
             b == Int8(3) ? UInt64(2) :
             b == Int8(2) ? UInt64(3) :
                            UInt64(99)
        if bb == UInt64(99)
            run  = Int32(0)
            kmer = UInt64(0)
            continue
        end
        kmer = ((kmer << 2) | bb) & mask
        run += Int32(1)
        run < k && continue

        rc    = _revcomp_swar(kmer, k)
        canon = kmer < rc ? kmer : rc

        lo = Int32(1); hi = n_kmers_in_index
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
        found_at == Int32(0) && continue

        ci = contig_of_kmer[found_at]
        if counter == Int32(0)
            candidate      = ci
            counter        = Int32(1)
            candidate_hits = Int32(1)
        elseif candidate == ci
            counter        += Int32(1)
            candidate_hits += Int32(1)
        else
            counter -= Int32(1)
            # candidate_hits unchanged: BM evicted one vote but the raw
            # count for the surviving candidate is still valid.
        end
    end

    if candidate == Int32(0)
        out_assignment[j] = Int32(0)
        out_hits[j]       = Int32(0)
        return nothing
    end

    # Fast path: candidate never changed mid-read (or changed but
    # accumulated enough hits). candidate_hits equals the full-read hit
    # count only when BM never evicted the candidate — i.e., candidate_hits
    # == counter (all votes were for the same contig). Otherwise run Pass 2
    # to get the true full-read count.
    if candidate_hits == counter
        # BM counter == candidate_hits → no evictions occurred → single
        # candidate throughout. candidate_hits IS the true hit count.
        if candidate_hits >= min_hits
            out_assignment[j] = candidate
            out_hits[j]       = candidate_hits
        else
            out_assignment[j] = Int32(0)
            out_hits[j]       = candidate_hits
        end
        return nothing
    end

    # --- Pass 2: full-read hit count for reads with mid-read candidate
    # switches (ambiguous or chimeric reads). Minority of reads in practice.
    hits = Int32(0)
    kmer = UInt64(0)
    run  = Int32(0)

    @inbounds for pos in Int32(1):L
        b  = seqs[pos, j]
        bb = b == Int8(1) ? UInt64(0) :
             b == Int8(4) ? UInt64(1) :
             b == Int8(3) ? UInt64(2) :
             b == Int8(2) ? UInt64(3) :
                            UInt64(99)
        if bb == UInt64(99)
            run  = Int32(0)
            kmer = UInt64(0)
            continue
        end
        kmer = ((kmer << 2) | bb) & mask
        run += Int32(1)
        run < k && continue

        rc    = _revcomp_swar(kmer, k)
        canon = kmer < rc ? kmer : rc

        lo = Int32(1); hi = n_kmers_in_index
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
        found_at == Int32(0) && continue

        contig_of_kmer[found_at] == candidate && (hits += Int32(1))
    end

    if hits >= min_hits
        out_assignment[j] = candidate
        out_hits[j]       = hits
    else
        out_assignment[j] = Int32(0)
        out_hits[j]       = hits
    end
    return nothing
end

function map_reads(contigs::AbstractVector{Contig},
                   seqs_h::AbstractMatrix{Int8},
                   lengths_h::AbstractVector{<:Integer};
                   k::Int = 31,
                   min_hits::Int = 3,
                   reads_per_batch::Int = 1_000_000,
                   verbose::Bool = false,
                   log_io = nothing)
    n_reads = length(lengths_h)
    n_contigs = length(contigs)
    n_reads == 0 && return (Int32[], Int32[])
    n_contigs == 0 && return (zeros(Int32, n_reads), zeros(Int32, n_reads))

    sorted_kmers, contig_of_kmer = build_contig_kmer_index(contigs, k)
    n_idx = length(sorted_kmers)
    n_idx == 0 && return (zeros(Int32, n_reads), zeros(Int32, n_reads))

    # Upload the k-mer index once; it stays resident in VRAM across all batches.
    # seqs_h is (max_read_len × n_reads) — too large for a single upload at HPC scale
    # (4B reads × 150 bp ≈ 600 GB). Process reads in column-wise batches instead.
    sorted_kmers_d   = JACC.to_device(sorted_kmers)
    contig_of_kmer_d = JACC.to_device(contig_of_kmer)

    out_assignment = zeros(Int32, n_reads)
    out_hits       = zeros(Int32, n_reads)

    batch_sz  = min(reads_per_batch, n_reads)
    n_batches = cld(n_reads, batch_sz)

    # Log progress every ~10% of batches (at least every 10, at most every 50).
    log_every = clamp(n_batches ÷ 10, 10, 50)
    t_map_start = time()

    for b in 1:n_batches
        r_start = (b - 1) * batch_sz + 1
        r_end   = min(b * batch_sz, n_reads)
        n_batch = r_end - r_start + 1

        seqs_d    = JACC.to_device(@view seqs_h[:, r_start:r_end])
        lengths_d = JACC.to_device(Int32.(lengths_h[r_start:r_end]))
        assign_d  = JACC.to_device(zeros(Int32, n_batch))
        hits_d    = JACC.to_device(zeros(Int32, n_batch))

        JACC.parallel_for(n_batch, _map_kernel_streaming!,
            seqs_d, lengths_d, sorted_kmers_d, contig_of_kmer_d,
            assign_d, hits_d,
            Int32(n_idx), Int32(k), Int32(min_hits))

        out_assignment[r_start:r_end] = JACC.to_host(assign_d)
        out_hits[r_start:r_end]       = JACC.to_host(hits_d)

        if verbose && (b % log_every == 0 || b == n_batches)
            elapsed = round(time() - t_map_start, digits=1)
            pct     = round(100 * r_end / n_reads, digits=1)
            log_println(log_io,
                "[mapping] batch $b/$n_batches  reads=$r_end/$n_reads ($(pct)%)  $(elapsed)s elapsed")
        end
    end

    return out_assignment, out_hits
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
