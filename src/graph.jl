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
    edge_offsets::Vector{Int32}        # length = 2*n_canonical + 1
    edge_targets::Vector{Int32}        # destination oriented node ID
    edge_weights::Vector{Int32}        # k-mer count
    edge_labels::Vector{UInt8}         # base 0..3 emitted by traversing
    edge_twins::Vector{Int32}          # index of this edge's twin edge
    edge_alive::BitVector
end

@inline twin(v::Integer) = isodd(v) ? Int32(v + 1) : Int32(v - 1)
@inline forward_id(canonical_idx::Integer) = Int32(2 * canonical_idx - 1)
@inline reverse_id(canonical_idx::Integer) = Int32(2 * canonical_idx)
@inline is_forward(v::Integer) = isodd(v)
@inline canonical_idx(v::Integer) = Int32((v + 1) ÷ 2)

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
                              node_id_map::Dict{UInt64,Int32})
    canon, is_fwd = _canonicalize_km1(raw_km1mer, km1)
    cidx = node_id_map[canon]
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

    # --- Pass 1: collect canonical (k-1)-mer set. ---
    node_set = Set{UInt64}()
    sizehint!(node_set, 2 * length(uniq_kmers))
    for kmer in uniq_kmers
        prefix, suffix, _, _ = _split_kmer(kmer, k)
        push!(node_set, _canonicalize_km1(prefix, km1)[1])
        push!(node_set, _canonicalize_km1(suffix, km1)[1])
    end
    canonical_kmers = sort!(collect(node_set))
    n_canonical = length(canonical_kmers)
    node_id_map = Dict{UInt64,Int32}(km => Int32(i)
                                     for (i, km) in enumerate(canonical_kmers))

    # --- Pass 2: build edge list with twin tracking. ---
    # Each k-mer -> 2 edges. Even-indexed in our temp arrays are forward,
    # odd-indexed are reverse twins (or vice versa); we'll record the
    # mapping explicitly so the twin pointer survives sorting.
    n_edges = 2 * length(uniq_kmers)
    edge_src   = Vector{Int32}(undef, n_edges)
    edge_dst   = Vector{Int32}(undef, n_edges)
    edge_w     = Vector{Int32}(undef, n_edges)
    edge_lbl   = Vector{UInt8}(undef, n_edges)
    twin_temp  = Vector{Int32}(undef, n_edges)   # original index of twin

    e = 0
    for (kmer, count) in zip(uniq_kmers, counts)
        prefix, suffix, first_base, last_base = _split_kmer(kmer, k)

        # Forward: oriented node of `prefix` -> oriented node of `suffix`
        src_fwd = _oriented_id(prefix, km1, node_id_map)
        dst_fwd = _oriented_id(suffix, km1, node_id_map)

        # Reverse twin: oriented node of rc(suffix) -> oriented node of rc(prefix)
        # Note: _oriented_id of rc(x) = twin of _oriented_id of x.
        src_rev = twin(dst_fwd)
        dst_rev = twin(src_fwd)

        e += 1
        i_fwd = e
        edge_src[i_fwd] = src_fwd
        edge_dst[i_fwd] = dst_fwd
        edge_w[i_fwd]   = count
        edge_lbl[i_fwd] = last_base

        e += 1
        i_rev = e
        edge_src[i_rev] = src_rev
        edge_dst[i_rev] = dst_rev
        edge_w[i_rev]   = count
        edge_lbl[i_rev] = _complement_base(first_base)

        twin_temp[i_fwd] = i_rev
        twin_temp[i_rev] = i_fwd
    end

    # Sort by source. We need to permute everything in lockstep AND
    # update twin pointers (since edge indices change).
    perm = sortperm(edge_src)
    edge_src_s = edge_src[perm]
    edge_dst_s = edge_dst[perm]
    edge_w_s   = edge_w[perm]
    edge_lbl_s = edge_lbl[perm]

    # inv_perm[old_idx] = new_idx
    inv_perm = Vector{Int32}(undef, n_edges)
    @inbounds for new_idx in 1:n_edges
        inv_perm[perm[new_idx]] = Int32(new_idx)
    end
    edge_twins_s = Vector{Int32}(undef, n_edges)
    @inbounds for new_idx in 1:n_edges
        old_idx = perm[new_idx]
        edge_twins_s[new_idx] = inv_perm[twin_temp[old_idx]]
    end

    # Build CSR offsets over 2*n_canonical oriented nodes.
    n_oriented = 2 * n_canonical
    edge_offsets = zeros(Int32, n_oriented + 1)
    for s in edge_src_s
        edge_offsets[s + 1] += 1
    end
    edge_offsets[1] = 1
    for i in 2:(n_oriented + 1)
        edge_offsets[i] += edge_offsets[i - 1]
    end

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
                      relative_threshold::Float64 = 0.05)
    total_removed = 0
    while true
        removed = _simplify_pass!(g, min_edge_weight, relative_threshold)
        total_removed += removed
        removed == 0 && break
    end
    return total_removed
end

function _simplify_pass!(g::DeBruijnGraph,
                         min_edge_weight::Int,
                         relative_threshold::Float64)
    # Global maximum live edge weight.
    max_w = Int32(0)
    for e in 1:length(g.edge_targets)
        g.edge_alive[e] || continue
        g.edge_weights[e] > max_w && (max_w = g.edge_weights[e])
    end
    max_w == 0 && return 0

    floor_w = max(Int32(min_edge_weight),
                  Int32(ceil(relative_threshold * Float64(max_w))))

    removed = 0
    for e in 1:length(g.edge_targets)
        g.edge_alive[e] || continue
        if g.edge_weights[e] < floor_w
            t = g.edge_twins[e]
            g.edge_alive[e] = false
            if g.edge_alive[t]
                g.edge_alive[t] = false
                removed += 1
            end
            removed += 1
        end
    end
    return removed
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
    edge_sources::Vector{Int32}         # source oriented node ID per unitig
    edge_targets::Vector{Int32}         # destination oriented node ID per unitig
    edge_weights::Vector{Float32}       # mean coverage along the unitig
    edge_sequences::Vector{Vector{UInt8}}  # appended bases (0..3) per unitig
    edge_twins::Vector{Int32}           # index of this unitig's twin
end

n_canonical(g::CompactedGraph) = length(g.canonical_kmers)
n_unitigs(g::CompactedGraph) = length(g.edge_targets)

"""
    compact_unitigs(g::DeBruijnGraph) -> CompactedGraph

Collapse maximal non-branching paths in the doubled-directed graph.
Each unitig has a twin unitig representing the reverse-complement path.
"""
function compact_unitigs(g::DeBruijnGraph)
    n_or = n_nodes(g)   # = 2 * n_canonical
    n_e  = length(g.edge_targets)

    # Compute in/out degrees over LIVE edges.
    indeg  = zeros(Int32, n_or)
    outdeg = zeros(Int32, n_or)
    edge_src = zeros(Int32, n_e)
    for v in 1:n_or
        for e in out_edges(g, v)
            edge_src[e] = Int32(v)
            g.edge_alive[e] || continue
            outdeg[v] += Int32(1)
            indeg[g.edge_targets[e]] += Int32(1)
        end
    end

    visited = falses(n_e)
    new_src     = Int32[]
    new_dst     = Int32[]
    new_weight  = Float32[]
    new_seqs    = Vector{Vector{UInt8}}()
    edge_to_unitig = zeros(Int32, n_e)   # which output unitig an edge ended up in

    function walk_forward(e_start::Integer)
        bases   = UInt8[g.edge_labels[e_start]]
        weights = Int32[g.edge_weights[e_start]]
        edges_in_walk = Int32[Int32(e_start)]
        visited[e_start] = true
        cur_dst = g.edge_targets[e_start]
        while indeg[cur_dst] == 1 && outdeg[cur_dst] == 1
            next_e = Int32(0)
            for e in out_edges(g, cur_dst)
                if g.edge_alive[e] && !visited[e]
                    next_e = Int32(e); break
                end
            end
            next_e == 0 && break
            visited[next_e] = true
            push!(edges_in_walk, next_e)
            push!(bases, g.edge_labels[next_e])
            push!(weights, g.edge_weights[next_e])
            cur_dst = g.edge_targets[next_e]
        end
        return bases, weights, cur_dst, edges_in_walk
    end

    function emit_walk(e_start::Integer)
        src = edge_src[e_start]
        bases, weights, final_dst, edges_in_walk = walk_forward(e_start)
        push!(new_src, src)
        push!(new_dst, final_dst)
        push!(new_weight, Float32(sum(weights)) / Float32(length(weights)))
        push!(new_seqs, bases)
        u_idx = Int32(length(new_seqs))
        for e in edges_in_walk
            edge_to_unitig[e] = u_idx
        end
    end

    # Pass 1: start walks at edges whose source is non-linear.
    for e_start in 1:n_e
        g.edge_alive[e_start] || continue
        visited[e_start] && continue
        src = edge_src[e_start]
        (indeg[src] != 1 || outdeg[src] != 1) || continue
        emit_walk(e_start)
    end

    # Pass 2: anything still unvisited belongs to a pure cycle of
    # linear nodes — pick an arbitrary edge as the start.
    for e_start in 1:n_e
        g.edge_alive[e_start] || continue
        visited[e_start] && continue
        emit_walk(e_start)
    end

    # Build twin mapping for unitigs.
    # The twin of a unitig is the unitig containing the twin of any of
    # its edges (they should all map to the same twin unitig if the
    # invariant holds).
    n_u = length(new_seqs)
    unitig_twins = Vector{Int32}(undef, n_u)
    @inbounds for u in 1:n_u
        # Find any edge in this unitig and look up where its twin landed.
        # We don't have edges_in_walk stored, but edge_to_unitig is the
        # inverse — scan for any edge belonging to u.
        twin_u = Int32(0)
        for e in 1:n_e
            edge_to_unitig[e] == u || continue
            t = g.edge_twins[e]
            if edge_to_unitig[t] != 0
                twin_u = edge_to_unitig[t]
                break
            end
        end
        unitig_twins[u] = twin_u
    end

    return CompactedGraph(g.k, copy(g.canonical_kmers),
                          new_src, new_dst, new_weight, new_seqs,
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
    suffix_chars = Vector{Char}(undef, length(g.edge_sequences[u]))
    @inbounds for (i, b) in enumerate(g.edge_sequences[u])
        suffix_chars[i] = b == 0 ? 'A' :
                          b == 1 ? 'C' :
                          b == 2 ? 'G' :
                                   'T'
    end
    return prefix * String(suffix_chars)
end
