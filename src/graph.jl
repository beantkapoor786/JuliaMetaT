# Doubled-directed de Bruijn graph in CSR layout.
#
# Each canonical (k-1)-mer becomes TWO oriented nodes (forward and reverse).
# Node IDs are 1-based; twin(v) flips between forward/reverse:
#   forward node for canonical index i: 2i - 1
#   reverse node for canonical index i: 2i
#   twin(v) = isodd(v) ? v+1 : v-1
#
# Each k-mer with count >= min_count produces TWO directed edges between
# ORIENTED nodes (not canonical nodes). The graph satisfies a twin
# invariant: for every edge u -> v with label b, there exists a twin
# edge twin(v) -> twin(u) with label complement(b) and same weight.
#
# This representation gives every TRUE interior node clean
# in_degree = out_degree = 1, enabling proper unitig compaction.

struct DeBruijnGraph
    k::Int
    canonical_kmers::Vector{UInt64}    # canonical (k-1)-mer per canonical index
    edge_offsets::Vector{Int64}        # length = 2*n_canonical + 1
    edge_targets::Vector{Int64}        # destination oriented node ID
    edge_weights::Vector{Int32}        # k-mer count
    edge_labels::Vector{UInt8}         # base 0..3 emitted by traversing
    edge_twins::Vector{Int64}          # index of this edge's twin edge
    edge_alive::BitVector
end

# Node/edge IDs use Int64: at HPC scale (billions of unique k-mers), the
# number of canonical (k-1)-mers and the number of edges both exceed the
# Int32 range (2^31-1) — see 5B-read-pair run with 17.9B unique k-mers.
@inline twin(v::Integer) = isodd(v) ? Int64(v + 1) : Int64(v - 1)
@inline forward_id(canonical_idx::Integer) = Int64(2 * canonical_idx - 1)
@inline reverse_id(canonical_idx::Integer) = Int64(2 * canonical_idx)
@inline is_forward(v::Integer) = isodd(v)
@inline canonical_idx(v::Integer) = Int64((v + 1) ÷ 2)

n_nodes(g::DeBruijnGraph) = 2 * length(g.canonical_kmers)

"""Edges leaving oriented node v as a UnitRange of edge indices."""
@inline function out_edges(g::DeBruijnGraph, v::Integer)
    g.edge_offsets[v]:(g.edge_offsets[v + 1] - 1)
end

function out_degree(g::DeBruijnGraph, v::Integer)
    cnt = 0
    for e in out_edges(g, v)
        g.edge_alive[e] && (cnt += 1)
    end
    cnt
end

# --- (k-1)-mer helpers ---

@inline function _split_kmer(kmer::UInt64, k::Int)
    km1 = k - 1
    mask_km1 = (UInt64(1) << (UInt64(2) * UInt64(km1))) - UInt64(1)
    first_base = UInt8((kmer >> (UInt64(2) * UInt64(km1))) & UInt64(0x3))
    last_base  = UInt8(kmer & UInt64(0x3))
    prefix = (kmer >> 2) & mask_km1
    suffix = kmer & mask_km1
    return prefix, suffix, first_base, last_base
end

# Canonicalize a (k-1)-mer. Returns (canonical_kmer, is_forward_orientation).
# is_forward_orientation = true if the input was already canonical
# (i.e., the smaller of itself and its rev-comp).
@inline function _canonicalize_km1(km1mer::UInt64, km1::Int)
    rc = _revcomp(km1mer, Int32(km1))
    if km1mer <= rc
        return km1mer, true
    else
        return rc, false
    end
end

@inline _complement_base(b::UInt8) = UInt8(b ⊻ 0x3)

# Map a raw (k-1)-mer (in some orientation) to its oriented node ID.
@inline function _oriented_id(raw_km1mer::UInt64, km1::Int,
                              canonical_kmers::Vector{UInt64})
    canon, is_fwd = _canonicalize_km1(raw_km1mer, km1)
    cidx = searchsortedfirst(canonical_kmers, canon)
    return is_fwd ? forward_id(cidx) : reverse_id(cidx)
end

"""
    build_graph(uniq_kmers, counts; k=31) -> DeBruijnGraph

Construct a doubled-directed de Bruijn graph in CSR form from canonical
k-mer counts. Each k-mer yields exactly two directed edges (forward and
its twin in reverse orientation) connecting oriented nodes.
"""
function build_graph(uniq_kmers::Vector{UInt64},
                     counts::Vector{Int32};
                     k::Int = 31)
    @assert length(uniq_kmers) == length(counts)
    km1 = k - 1

    # --- Pass 1: collect canonical (k-1)-mer set via sort+dedup. ---
    # Flat preallocated array + sort is ~3-5x faster than Set{UInt64} at
    # this scale due to cache-friendly sequential access.
    n_kmers = length(uniq_kmers)
    km1mers = Vector{UInt64}(undef, 2 * n_kmers)
    @inbounds for i in 1:n_kmers
        prefix, suffix, _, _ = _split_kmer(uniq_kmers[i], k)
        km1mers[2i-1] = _canonicalize_km1(prefix, km1)[1]
        km1mers[2i]   = _canonicalize_km1(suffix, km1)[1]
    end
    sort!(km1mers)
    # Manual dedup in-place, then copy down to an exactly-sized array —
    # resize! alone would keep the 2*n_kmers buffer reserved underneath.
    n_canonical = 0
    @inbounds for i in 1:(2 * n_kmers)
        if i == 1 || km1mers[i] != km1mers[i-1]
            n_canonical += 1
            km1mers[n_canonical] = km1mers[i]
        end
    end
    canonical_kmers = km1mers[1:n_canonical]
    km1mers = UInt64[]   # release the 2*n_kmers oversized buffer before Pass 2
    GC.gc()

    # --- Pass 2: CSR construction without materializing per-edge unsorted
    # arrays or a permutation. Each k-mer i contributes exactly one forward
    # edge and its twin reverse edge; since both are scattered within the
    # same loop iteration, their final (sorted) positions are known to each
    # other directly — no inv_perm/perm/twin_temp bookkeeping required.
    # This trades the prior loop's thread-parallelism for ~5 fewer
    # n_edges-sized arrays alive at once, which is the dominant memory cost
    # at HPC scale (billions of k-mers).
    n_edges    = 2 * n_kmers
    n_oriented = 2 * n_canonical

    # Histogram pass: out-degree per oriented source node.
    edge_offsets = zeros(Int64, n_oriented + 1)
    @inbounds for i in 1:n_kmers
        prefix, suffix, _, _ = _split_kmer(uniq_kmers[i], k)
        src_fwd = _oriented_id(prefix, km1, canonical_kmers)
        dst_fwd = _oriented_id(suffix, km1, canonical_kmers)
        src_rev = twin(dst_fwd)
        edge_offsets[src_fwd + 1] += Int64(1)
        edge_offsets[src_rev + 1] += Int64(1)
    end
    edge_offsets[1] = Int64(1)
    @inbounds for i in 2:(n_oriented + 1)
        edge_offsets[i] += edge_offsets[i - 1]
    end

    # Scatter pass: write both directions of each k-mer's edge pair directly
    # into their final sorted slots, linking twins inline.
    edge_dst_s   = Vector{Int64}(undef, n_edges)
    edge_w_s     = Vector{Int32}(undef, n_edges)
    edge_lbl_s   = Vector{UInt8}(undef, n_edges)
    edge_twins_s = Vector{Int64}(undef, n_edges)
    # Use edge_offsets itself as the mutable write-cursor (no separate `pos`
    # copy — at HPC scale that's another n_oriented-sized Int64 array, the
    # difference between fitting in RAM and OOMing). After the scatter,
    # edge_offsets[v] no longer holds node v's start offset — it holds the
    # position just past node v's last written edge, i.e. node (v+1)'s
    # original start. We restore proper 1-based CSR starts below via an
    # in-place backward shift, using that exact relationship.
    @inbounds for i in 1:n_kmers
        prefix, suffix, first_base, last_base = _split_kmer(uniq_kmers[i], k)
        src_fwd = _oriented_id(prefix, km1, canonical_kmers)
        dst_fwd = _oriented_id(suffix, km1, canonical_kmers)
        src_rev = twin(dst_fwd)
        dst_rev = twin(src_fwd)

        p_fwd = edge_offsets[src_fwd]; edge_offsets[src_fwd] = p_fwd + Int64(1)
        p_rev = edge_offsets[src_rev]; edge_offsets[src_rev] = p_rev + Int64(1)

        edge_dst_s[p_fwd]   = dst_fwd
        edge_w_s[p_fwd]     = counts[i]
        edge_lbl_s[p_fwd]   = last_base
        edge_twins_s[p_fwd] = p_rev

        edge_dst_s[p_rev]   = dst_rev
        edge_w_s[p_rev]     = counts[i]
        edge_lbl_s[p_rev]   = _complement_base(first_base)
        edge_twins_s[p_rev] = p_fwd
    end

    # Restore proper CSR start offsets: edge_offsets[v] currently holds
    # original_start[v+1] (the cursor ran exactly outdeg[v] times past it).
    # Shift right by one in place, highest index first to avoid clobbering
    # values not yet read.
    @inbounds for v in n_oriented:-1:1
        edge_offsets[v + 1] = edge_offsets[v]
    end
    edge_offsets[1] = Int64(1)

    edge_alive = trues(n_edges)

    return DeBruijnGraph(k, canonical_kmers, edge_offsets,
                         edge_dst_s, edge_w_s, edge_lbl_s,
                         edge_twins_s, edge_alive)
end

# --- Edge simplification ---
#
# Twin-paired pruning: when we kill an edge, we kill its twin too.
# This preserves the twin invariant at every step.

"""
    remove_tips!(g; min_edge_weight=2, relative_threshold=0.05, ...) -> n_removed

Prune low-coverage edges, in twin pairs. Returns count of edges marked dead
(includes both members of every twin pair, so the count is even).

Threshold: `floor = max(min_edge_weight, ceil(relative_threshold * max_weight))`.
Edges with weight < floor are killed (along with their twins).

`tip_length` is accepted for API stability but unused.
"""
function remove_tips!(g::DeBruijnGraph;
                      tip_length::Int = 2 * g.k,
                      min_edge_weight::Int = 2,
                      relative_threshold::Float64 = 0.05,
                      verbose::Bool = false,
                      log_io = nothing)
    n   = n_nodes(g)
    n_e = length(g.edge_targets)

    # Upload static arrays once — weights, twins, offsets never change after
    # build_graph. Only edge_alive round-trips each pass.
    weights_d = JACC.to_device(g.edge_weights)
    twins_d   = JACC.to_device(g.edge_twins)
    offsets_d = JACC.to_device(g.edge_offsets)
    # Scratch buffers for relative prune (reused across convergence iterations).
    src_d     = JACC.to_device(zeros(Int64, n_e))
    max_d     = JACC.to_device(zeros(Int32, n))

    # Floor prune is idempotent: edge weights never change, so one pass is enough.
    t_floor = @elapsed removed_floor = _floor_prune!(g, min_edge_weight)
    total_removed = removed_floor

    # Upload alive once after floor prune. Relative prune keeps it on device
    # across iterations — no re-upload per iteration, only one to_host per iter.
    alive_d = JACC.to_device(UInt8.(g.edge_alive))

    # Relative prune loops until convergence.
    t_rel = 0.0
    iters = 0
    while true
        tr = @elapsed r = _relative_prune!(g, relative_threshold,
                                           offsets_d, weights_d, twins_d,
                                           src_d, max_d, alive_d)
        t_rel += tr
        iters += 1
        total_removed += r
        r == 0 && break
    end
    if verbose
        log_println(log_io, "[graph/prune] floor=$(round(t_floor,digits=3))s  " *
                "relative=$(round(t_rel,digits=3))s  iters=$iters")
    end
    return total_removed
end

# JACC kernel: one thread per oriented node v. Scans v's own out-edge
# range (disjoint across nodes — no races) to record each edge's source
# and the strongest live out-edge weight at v.
function _node_pass_kernel!(v, edge_offsets, edge_weights, edge_alive, edge_src, max_out_w)
    lo = edge_offsets[v]
    hi = edge_offsets[v + Int64(1)] - Int64(1)
    m = Int32(0)
    @inbounds for e in lo:hi
        edge_src[e] = v
        if edge_alive[e] != UInt8(0)
            w = edge_weights[e]
            w > m && (m = w)
        end
    end
    max_out_w[v] = m
    return nothing
end

# JACC kernel: one thread per edge. Kills edges (and twins) whose weight
# falls below relative_threshold × the local max out-weight at their
# source node. Same idempotent twin-write safety as _floor_prune_kernel!.
function _relative_prune_kernel!(e, edge_weights, edge_twins, edge_src,
                                  max_out_w, relative_threshold::Float32, edge_alive)
    @inbounds begin
        edge_alive[e] == UInt8(0) && return nothing
        local_max = max_out_w[edge_src[e]]
        local_max == Int32(0) && return nothing
        # Float32 throughout — Metal shaders reject Float64 entirely.
        if Float32(edge_weights[e]) < relative_threshold * Float32(local_max)
            edge_alive[e] = UInt8(0)
            edge_alive[edge_twins[e]] = UInt8(0)
        end
    end
    return nothing
end

# Run the local relative-threshold pruning pass entirely on-device:
# (1) per-node pass computes edge_src + local max out-weight,
# (2) per-edge pass kills edges below relative_threshold × local max.
# alive_d is owned by remove_tips! and persists across iterations — this
# function reads/writes it in-place without re-uploading from host.
# Writes results back into g.edge_alive; returns count of newly-dead edges.
function _relative_prune!(g::DeBruijnGraph, relative_threshold::Float64,
                          offsets_d, weights_d, twins_d, src_d, max_d, alive_d)
    n   = n_nodes(g)
    n_e = length(g.edge_targets)

    JACC.parallel_for(n, _node_pass_kernel!, offsets_d, weights_d, alive_d, src_d, max_d)
    # Pass Float32 explicitly — Metal shaders reject Float64 kernel args.
    JACC.parallel_for(n_e, _relative_prune_kernel!,
        weights_d, twins_d, src_d, max_d, Float32(relative_threshold), alive_d)

    alive_host = JACC.to_host(alive_d)
    killed = 0
    @inbounds for e in 1:n_e
        was = g.edge_alive[e]
        now = alive_host[e] != UInt8(0)
        g.edge_alive[e] = now
        was & !now && (killed += 1)
    end
    return killed
end

# Floor prune runs entirely on CPU threads — no Metal dispatch overhead,
# no device transfers, no BitVector word-level races.
# Strategy: mark kills into a Vector{UInt8} scratch buffer (one byte per
# edge, thread-safe for idempotent zeroing), then fold back into edge_alive.
function _floor_prune!(g::DeBruijnGraph, min_edge_weight::Int)
    n_e = length(g.edge_targets)
    w   = g.edge_weights
    tw  = g.edge_twins
    mw  = Int32(min_edge_weight)

    # Scratch buffer: 1 = alive, 0 = dead. Start fully alive.
    scratch = ones(UInt8, n_e)

    # Each thread owns a disjoint set of edges (by index). The only shared
    # writes are to twin slots, but twin pairs write the same value (0),
    # so concurrent writes are idempotent — no race can corrupt state.
    Threads.@threads for e in 1:n_e
        @inbounds if w[e] < mw
            scratch[e]      = UInt8(0)
            scratch[tw[e]]  = UInt8(0)
        end
    end

    killed = 0
    @inbounds for e in 1:n_e
        if scratch[e] == UInt8(0) && g.edge_alive[e]
            g.edge_alive[e] = false
            killed += 1
        end
    end
    return killed
end


# --- Unitig compaction ---
#
# On a doubled-directed graph, a "linear interior" oriented node has
# in_degree == out_degree == 1 considering only live edges.
# Compaction walks forward from edges whose source is not linear,
# accumulating bases until reaching a non-linear destination.

struct CompactedGraph
    k::Int
    canonical_kmers::Vector{UInt64}     # same canonical (k-1)-mer table as parent
    edge_sources::Vector{Int64}         # source oriented node ID per unitig
    edge_targets::Vector{Int64}         # destination oriented node ID per unitig
    edge_weights::Vector{Float32}       # mean coverage along the unitig
    edge_sequences::Vector{Vector{UInt8}}  # appended bases (0..3) per unitig
    edge_twins::Vector{Int32}           # index of this unitig's twin (unitig count stays Int32-bounded)
end

n_canonical(g::CompactedGraph) = length(g.canonical_kmers)
n_unitigs(g::CompactedGraph) = length(g.edge_targets)

"""
    compact_unitigs(g::DeBruijnGraph) -> CompactedGraph

Collapse maximal non-branching paths in the doubled-directed graph.
Each unitig has a twin unitig representing the reverse-complement path.
"""
# JACC kernel: one thread per oriented node v. Scans v's own out-edge
# range (disjoint across nodes — no races) to label edge sources and
# count live out-degree.
function _node_degree_pass_kernel!(v, edge_offsets, edge_alive, edge_src, outdeg)
    lo = edge_offsets[v]
    hi = edge_offsets[v + Int64(1)] - Int64(1)
    cnt = Int32(0)
    @inbounds for e in lo:hi
        edge_src[e] = v
        edge_alive[e] != UInt8(0) && (cnt += Int32(1))
    end
    outdeg[v] = cnt
    return nothing
end

# JACC kernel: in-degree via the twin invariant — every live edge ending
# at v has a twin edge starting at twin(v) (twin-paired pruning keeps
# aliveness in sync), so indeg(v) == outdeg(twin(v)). This turns the
# in-degree scatter (which would race across source nodes) into a
# race-free elementwise gather.
function _indeg_from_twin_kernel!(v, outdeg, indeg)
    @inbounds indeg[v] = outdeg[twin(v)]
    return nothing
end

# Compute (indeg, outdeg, edge_src) for the doubled-directed graph as
# JACC parallel_for dispatches — see kernels above for the race-free
# decomposition (per-node out-edge ranges + twin-invariant gather).
function _compute_degrees_and_edge_src(g::DeBruijnGraph)
    n   = n_nodes(g)
    n_e = length(g.edge_targets)
    edge_alive_u8 = UInt8.(g.edge_alive)

    offsets_d = JACC.to_device(g.edge_offsets)
    alive_d   = JACC.to_device(edge_alive_u8)
    src_d     = JACC.to_device(zeros(Int64, n_e))
    outdeg_d  = JACC.to_device(zeros(Int32, n))
    JACC.parallel_for(n, _node_degree_pass_kernel!, offsets_d, alive_d, src_d, outdeg_d)

    indeg_d = JACC.to_device(zeros(Int32, n))
    JACC.parallel_for(n, _indeg_from_twin_kernel!, outdeg_d, indeg_d)

    return JACC.to_host(indeg_d), JACC.to_host(outdeg_d), JACC.to_host(src_d)
end

function compact_unitigs(g::DeBruijnGraph; verbose::Bool = false, log_io = nothing)
    n_e = length(g.edge_targets)

    t_deg = @elapsed indeg, outdeg, edge_src = _compute_degrees_and_edge_src(g)

    # Phase 1: O(E) scan for non-cycle chain starts.
    # A chain start is any alive edge whose source node is non-linear
    # (indeg != 1 or outdeg != 1). Chains from distinct non-linear sources
    # are structurally disjoint (a linear interior node has exactly one live
    # predecessor, placing it in exactly one maximal chain), so no sequential
    # chain-marking is needed to avoid duplicate starts.
    t_starts = @elapsed begin
        starts = Int64[]
        sizehint!(starts, n_e ÷ 4)
        for e in 1:n_e
            g.edge_alive[e] || continue
            src = edge_src[e]
            (indeg[src] != 1 || outdeg[src] != 1) || continue
            push!(starts, Int64(e))
        end
    end
    n_nc = length(starts)

    # Phase 2: parallel walks from non-cycle starts.
    # Structurally disjoint chains -> threads write only to their own slot.
    walk_src    = Vector{Int64}(undef, n_nc)
    walk_dst    = Vector{Int64}(undef, n_nc)
    walk_weight = Vector{Float32}(undef, n_nc)
    walk_seqs   = Vector{Vector{UInt8}}(undef, n_nc)
    walk_edges  = Vector{Vector{Int64}}(undef, n_nc)

    t_walks = @elapsed Threads.@threads for i in 1:n_nc
        e_start = starts[i]
        bases         = UInt8[g.edge_labels[e_start]]
        edges_in_walk = Int64[e_start]
        weight_sum    = Int32(g.edge_weights[e_start])
        weight_count  = Int32(1)
        cur_dst = g.edge_targets[e_start]
        # local_visited not needed: linear chains (indeg==outdeg==1) can never
        # loop back — a back-edge to any prior node would give that node indeg>1,
        # breaking the while condition before we could revisit it.
        while indeg[cur_dst] == 1 && outdeg[cur_dst] == 1
            next_e = Int64(0)
            for e in out_edges(g, cur_dst)
                g.edge_alive[e] && (next_e = Int64(e); break)
            end
            next_e == 0 && break
            push!(edges_in_walk, next_e)
            push!(bases, g.edge_labels[next_e])
            weight_sum   += g.edge_weights[next_e]
            weight_count += Int32(1)
            cur_dst = g.edge_targets[next_e]
        end
        walk_src[i]    = edge_src[e_start]
        walk_dst[i]    = cur_dst
        walk_weight[i] = Float32(weight_sum) / Float32(weight_count)
        walk_seqs[i]   = bases
        walk_edges[i]  = edges_in_walk
    end

    # Phase 3: pure-cycle detection (sequential; rare in prokaryotic DBGs).
    # Build a covered set from Phase 2 walks, then sweep for uncovered alive
    # edges — each is a cycle start. Walk it to full coverage sequentially.
    t_cycles = @elapsed begin
        covered = falses(n_e)
        for i in 1:n_nc
            for e in walk_edges[i]
                @inbounds covered[e] = true
            end
        end
        for e in 1:n_e
            @inbounds g.edge_alive[e] || continue
            @inbounds covered[e] && continue
            cy_bases      = UInt8[g.edge_labels[e]]
            cy_edges      = Int64[e]
            cy_weight_sum = Int32(g.edge_weights[e])
            cy_weight_cnt = Int32(1)
            covered[e] = true
            cur_dst = g.edge_targets[e]
            while indeg[cur_dst] == 1 && outdeg[cur_dst] == 1
                next_e = 0
                for ee in out_edges(g, cur_dst)
                    g.edge_alive[ee] && !covered[ee] || continue
                    next_e = ee; break
                end
                next_e == 0 && break
                covered[next_e] = true
                push!(cy_edges, Int64(next_e))
                push!(cy_bases, g.edge_labels[next_e])
                cy_weight_sum += g.edge_weights[next_e]
                cy_weight_cnt += Int32(1)
                cur_dst = g.edge_targets[next_e]
            end
            push!(starts, Int64(e))
            push!(walk_src, edge_src[e])
            push!(walk_dst, cur_dst)
            push!(walk_weight, Float32(cy_weight_sum) / Float32(cy_weight_cnt))
            push!(walk_seqs, cy_bases)
            push!(walk_edges, cy_edges)
        end
    end
    n_walks = length(starts)

    t_twin = @elapsed begin
        edge_to_unitig = zeros(Int32, n_e)
        for i in 1:n_walks
            for e in walk_edges[i]
                edge_to_unitig[e] = Int32(i)
            end
        end
        unitig_twins = zeros(Int32, n_walks)
        @inbounds for u in 1:n_walks
            for e in walk_edges[u]
                t = g.edge_twins[e]
                if edge_to_unitig[t] != 0
                    unitig_twins[u] = edge_to_unitig[t]
                    break
                end
            end
        end
    end

    if verbose
        log_println(log_io, "[graph/compact] degrees=$(round(t_deg,digits=3))s  " *
                "starts=$(round(t_starts,digits=3))s  " *
                "walks=$(round(t_walks,digits=3))s  " *
                "cycles=$(round(t_cycles,digits=3))s  " *
                "twin=$(round(t_twin,digits=3))s  " *
                "n_unitigs=$n_walks")
    end

    return CompactedGraph(g.k, copy(g.canonical_kmers),
                          walk_src, walk_dst, walk_weight, walk_seqs,
                          unitig_twins)
end

# --- Sequence reconstruction ---

# Decode an oriented node back to a (k-1)-mer string.
function node_to_string(g::DeBruijnGraph, v::Integer)
    cidx = canonical_idx(v)
    canonical = g.canonical_kmers[cidx]
    raw = is_forward(v) ? canonical : _revcomp(canonical, Int32(g.k - 1))
    kmer_to_string(raw, g.k - 1)
end
node_to_string(g::CompactedGraph, v::Integer) = begin
    cidx = canonical_idx(v)
    canonical = g.canonical_kmers[cidx]
    raw = is_forward(v) ? canonical : _revcomp(canonical, Int32(g.k - 1))
    kmer_to_string(raw, g.k - 1)
end

"""
    unitig_sequence(g::CompactedGraph, u::Integer) -> String

Reconstruct the full DNA sequence of a unitig: the source node's
(k-1)-mer followed by the unitig's appended bases. Length = (k-1) + length(bases).
"""
function unitig_sequence(g::CompactedGraph, u::Integer)
    src = g.edge_sources[u]
    prefix = node_to_string(g, src)
    suffix_bytes = Vector{UInt8}(undef, length(g.edge_sequences[u]))
    @inbounds for (i, b) in enumerate(g.edge_sequences[u])
        suffix_bytes[i] = b == 0x00 ? UInt8('A') :
                          b == 0x01 ? UInt8('C') :
                          b == 0x02 ? UInt8('G') :
                                      UInt8('T')
    end
    return prefix * String(suffix_bytes)
end
