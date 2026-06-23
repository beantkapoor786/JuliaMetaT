# KMC3 backend for k-mer counting at scale.
# DESIGN.md §III.3 — replaces in-memory sort-and-tally with KMC's
# external-memory minimizer-bucketed counter. The in-memory
# count_kmers in kmers.jl is retained as a reference implementation.

using Mmap

const KMC_BINARY       = "kmc"
const KMC_TOOLS_BINARY = "kmc_tools"

"""
    check_kmc()

Verify KMC and kmc_tools are on PATH. Returns the kmc version string,
or throws ErrorException with an install hint.
"""
function check_kmc()
    for bin in (KMC_BINARY, KMC_TOOLS_BINARY)
        if Sys.which(bin) === nothing
            error("""
                  $bin not found on PATH.
                  Install KMC3:
                    macOS:  brew install kmc
                    Linux:  conda install -c bioconda kmc
                  """)
        end
    end
    out = read(pipeline(`$KMC_BINARY`, stderr=devnull), String)
    m = match(r"KMC\)?\s+ver\.\s+(\S+)"i, out)
    return m === nothing ? "unknown" : m.captures[1]
end

# 2-bit packing matching kmers.jl: A=00, C=01, G=10, T=11.
# Note: this is a DIFFERENT mapping than the encoded-read alphabet
# (which uses A=1,T=2,G=3,C=4). The string→canonical conversion below
# produces a UInt64 in the same packing as count_kmers' output, so
# downstream consumers (build_graph) don't care which counter produced it.
@inline function _base2bit(b::UInt8)
    b == UInt8('A') ? UInt64(0) :
    b == UInt8('C') ? UInt64(1) :
    b == UInt8('G') ? UInt64(2) :
    b == UInt8('T') ? UInt64(3) :
                      UInt64(99)   # N or other → signal invalid
end

"""
    parse_kmer_ascii(buf, start, k) -> UInt64 or nothing

Read k bases from `buf` starting at byte index `start` (1-indexed),
return the canonical (lex-min of forward vs reverse complement)
packed UInt64. Returns `nothing` if any base is non-ACGT.
"""
@inline function parse_kmer_ascii(buf::AbstractVector{UInt8}, start::Int, k::Int)
    fwd = UInt64(0)
    @inbounds for i in 0:k-1
        b = _base2bit(buf[start + i])
        b == UInt64(99) && return nothing
        fwd = (fwd << 2) | b
    end
    mask = (UInt64(1) << (2 * k)) - UInt64(1)
    x = (~fwd) & mask
    rc = UInt64(0)
    @inbounds for i in 0:k-1
        pair = (x >> (2 * i)) & UInt64(0x3)
        rc |= pair << (2 * (k - 1 - i))
    end
    return fwd < rc ? fwd : rc
end

# Find the next newline at or after `pos`. Returns index of '\n', or
# `lastindex(buf) + 1` if none.
@inline function _next_newline(buf::AbstractVector{UInt8}, pos::Int)
    n = length(buf)
    @inbounds while pos <= n && buf[pos] != UInt8('\n')
        pos += 1
    end
    return pos
end

# Parse one "ACGT...<TAB><count>\n" line starting at `pos`.
# Returns (kmer::UInt64, count::Int32, next_pos::Int).
@inline function _parse_line(buf::AbstractVector{UInt8}, pos::Int, k::Int)
    kmer = parse_kmer_ascii(buf, pos, k)
    # KMC dump format: "<kmer>\t<count>\n"
    tab = pos + k
    # Defensive: skip if format unexpected
    @inbounds if tab > length(buf) || buf[tab] != UInt8('\t')
        nl = _next_newline(buf, pos)
        return (UInt64(0), Int32(0), nl + 1, false)
    end
    cnt = Int32(0)
    i = tab + 1
    n = length(buf)
    @inbounds while i <= n && buf[i] != UInt8('\n')
        d = buf[i] - UInt8('0')
        cnt = cnt * Int32(10) + Int32(d)
        i += 1
    end
    return (kmer === nothing ? UInt64(0) : kmer, cnt, i + 1, kmer !== nothing)
end

"""
    _parse_dump_chunk(buf, byte_start, byte_end, k) -> (uniq, counts)

Parse the byte range [byte_start, byte_end] of `buf`, aligning to
line boundaries. Returns sorted (uniq, counts) for this chunk.

Boundary handling:
- If byte_start > 1, skip the partial line at the start (the previous
  chunk owns it).
- Always read to the end of the line containing byte_end (so the next
  chunk doesn't have to handle it).
"""
function _parse_dump_chunk(buf::Vector{UInt8}, byte_start::Int, byte_end::Int, k::Int)
    pos = byte_start
    if byte_start > 1
        # Skip partial line; previous chunk handled it
        nl = _next_newline(buf, byte_start)
        pos = nl + 1
    end
    stop = _next_newline(buf, byte_end)
    # stop is the newline AT or after byte_end; we parse up to and
    # including the line that ends at `stop`.

    # KMC produces sorted output, so within a chunk k-mers come out sorted
    # by their ASCII representation. The canonicalization in
    # parse_kmer_ascii may reorder some of them (when forward differs from
    # rev-comp in lex order), so we sort the chunk at the end to be safe.
    uniq = UInt64[]
    counts = Int32[]
    sizehint!(uniq,   max(1, (stop - pos) ÷ 40))   # ~40 bytes/line at k=31
    sizehint!(counts, length(uniq))

    while pos <= stop
        kmer, cnt, next_pos, ok = _parse_line(buf, pos, k)
        if ok
            push!(uniq, kmer)
            push!(counts, cnt)
        end
        pos = next_pos
    end

    # Sort by k-mer (re-canonicalization may have permuted things)
    if !issorted(uniq)
        p = sortperm(uniq)
        uniq   = uniq[p]
        counts = counts[p]
    end
    return uniq, counts
end

# Binary min-heap sift-down for (kmer, chunk_idx) pairs ordered by kmer.
@inline function _heap_sift_down!(heap::Vector{Tuple{UInt64,Int}}, i::Int)
    n = length(heap)
    @inbounds while true
        left  = 2i
        right = 2i + 1
        smallest = i
        left  <= n && heap[left][1]  < heap[smallest][1] && (smallest = left)
        right <= n && heap[right][1] < heap[smallest][1] && (smallest = right)
        smallest == i && break
        heap[i], heap[smallest] = heap[smallest], heap[i]
        i = smallest
    end
end

# K-way merge of N sorted (uniq, counts) chunks using a binary min-heap.
# O(N log K) vs the naive O(N×K) linear scan. Same-key entries are summed
# (canonicalization can merge a forward and rev-comp k-mer from different chunks).
function _merge_sorted_chunks(chunks::Vector{Tuple{Vector{UInt64},Vector{Int32}}})
    n_chunks = length(chunks)
    n_chunks == 0 && return UInt64[], Int32[]
    n_chunks == 1 && return chunks[1]

    total = sum(length(c[1]) for c in chunks)
    uniq   = Vector{UInt64}(undef, total)
    counts = Vector{Int32}(undef, total)

    idx  = ones(Int, n_chunks)
    lens = [length(c[1]) for c in chunks]

    # Seed heap with first element from each non-empty chunk.
    heap = Tuple{UInt64,Int}[]
    sizehint!(heap, n_chunks)
    for c in 1:n_chunks
        lens[c] > 0 && push!(heap, (chunks[c][1][1], c))
    end
    for i in (length(heap) ÷ 2):-1:1
        _heap_sift_down!(heap, i)
    end

    out_i = 0
    @inbounds while !isempty(heap)
        best_val = heap[1][1]
        cnt_sum  = Int32(0)

        # Drain all heap entries with this kmer value (usually just one;
        # duplicates arise only when canonicalization maps two chunks to
        # the same canonical k-mer).
        while !isempty(heap) && heap[1][1] == best_val
            c = heap[1][2]
            cnt_sum += chunks[c][2][idx[c]]
            idx[c] += 1
            if idx[c] <= lens[c]
                heap[1] = (chunks[c][1][idx[c]], c)
                _heap_sift_down!(heap, 1)
            else
                heap[1] = heap[end]
                pop!(heap)
                !isempty(heap) && _heap_sift_down!(heap, 1)
            end
        end

        out_i += 1
        uniq[out_i]   = best_val
        counts[out_i] = cnt_sum
    end

    resize!(uniq,   out_i)
    resize!(counts, out_i)
    return uniq, counts
end

"""
    count_kmers_kmc(r1, r2=nothing; k=31, min_count=2, max_count=1_000_000,
                    threads=Threads.nthreads(), memory_gb=8,
                    workdir=mktempdir(), verbose=false)
        -> (uniq::Vector{UInt64}, counts::Vector{Int32})

Count canonical k-mers via KMC3. Output is sorted by canonical k-mer
and is byte-compatible with `count_kmers` from kmers.jl, so the
returned arrays can be passed directly to `build_graph`.

Performance: KMC handles counting (external-memory, parallelized).
Parsing of KMC's text dump is parallelized across `threads`.

Arguments:
- `r1`, `r2`: paths to FASTQ files (gzip OK; KMC handles .gz natively).
- `k`: k-mer length. Must match what build_graph expects.
- `min_count`: discard k-mers with count < this. Passed to KMC's -ci.
- `max_count`: counter ceiling. Set high enough that no real k-mer
  hits it. KMC silently caps counts at this value if hit.
- `threads`: parallelism for both KMC and dump parsing.
- `memory_gb`: KMC's RAM budget. Passed to KMC's -m.
- `workdir`: scratch directory. Cleaned up after.
- `verbose`: print stage timings to stdout.
"""
function count_kmers_kmc(r1::AbstractString, r2::Union{AbstractString,Nothing}=nothing;
                         k::Int = 31,
                         min_count::Int = 2,
                         max_count::Int = 1_000_000,
                         threads::Int = Threads.nthreads(),
                         memory_gb::Int = 8,
                         workdir::AbstractString = mktempdir(),
                         verbose::Bool = false,
                         log_io = nothing)
    check_kmc()
    mkpath(workdir)

    try
        # Step 1: build file-of-filenames
        fof = joinpath(workdir, "inputs.txt")
        open(fof, "w") do io
            println(io, r1)
            r2 === nothing || println(io, r2)
        end

        db = joinpath(workdir, "db")
        dump_path = joinpath(workdir, "dump.txt")

        # KMC's own stdout/stderr are captured to the pipeline log (rather than
        # discarded) so a process failure (e.g. exit code 1 from a full scratch
        # disk) is diagnosable from the log instead of just "ProcessExited(1)".
        kmc_out = log_io === nothing ? devnull : log_io
        kmc_err = log_io === nothing ? devnull : log_io
        log_io !== nothing && flush(log_io)

        # Step 2: run KMC
        t_count = @elapsed run(pipeline(`$KMC_BINARY
            -k$k -ci$min_count -cs$max_count
            -t$threads -m$memory_gb -fq
            @$fof $db $workdir`, stdout=kmc_out, stderr=kmc_err))
        verbose && log_println(log_io, "[kmc]     count: $(round(t_count, digits=2))s")

        # Step 3: dump to text, sorted
        log_io !== nothing && flush(log_io)
        t_dump = @elapsed run(pipeline(`$KMC_TOOLS_BINARY transform $db dump -s $dump_path`,
                                       stdout=kmc_out, stderr=kmc_err))
        verbose && log_println(log_io, "[kmc]     dump:  $(round(t_dump, digits=2))s")

        # Step 4: parallel parse via mmap
        t_parse = @elapsed begin
            buf = open(dump_path, "r") do io
                Mmap.mmap(io, Vector{UInt8})
            end

            n_bytes = length(buf)
            if n_bytes == 0
                return UInt64[], Int32[]
            end

            # Divide into thread chunks
            chunk_size = max(1, n_bytes ÷ threads)
            chunks_results = Vector{Tuple{Vector{UInt64},Vector{Int32}}}(undef, threads)

            Threads.@threads for t in 1:threads
                bstart = (t - 1) * chunk_size + 1
                bend   = t == threads ? n_bytes : t * chunk_size
                chunks_results[t] = _parse_dump_chunk(buf, bstart, bend, k)
            end

            uniq, counts = _merge_sorted_chunks(chunks_results)
        end
        verbose && log_println(log_io, "[kmc]     parse: $(round(t_parse, digits=2))s ($threads threads)")

        # Sanity: detect max_count clipping
        if !isempty(counts) && maximum(counts) == max_count
            @warn "Some k-mer counts hit max_count=$max_count; coverage measurements may be capped. Consider raising max_count."
        end

        return uniq, counts
    finally
        # Always clean up the workdir
        try
            rm(workdir; recursive=true, force=true)
        catch
        end
    end
end
