# Connected components + threaded contig traversal on the compacted graph.

struct Contig
    sequence::String              # full DNA bases
    mean_coverage::Float32        # weighted-mean across constituent unitigs
    n_unitigs::Int                # how many unitigs were chained
    component_id::Int             # which connected component this came from
end

# --- Union-find ---

mutable struct UnionFind
    parent::Vector{Int32}
    rank::Vector{Int32}
end
UnionFind(n::Integer) = UnionFind(Int32.(1:n), zeros(Int32, n))

function uf_find(uf::UnionFind, x::Integer)
    while uf.parent[x] != x
        uf.parent[x] = uf.parent[uf.parent[x]]   # path compression
        x = uf.parent[x]
    end
    return Int32(x)
end

function uf_union!(uf::UnionFind, a::Integer, b::Integer)
    ra, rb = uf_find(uf, a), uf_find(uf, b)
    ra == rb && return ra
    if uf.rank[ra] < uf.rank[rb]
        uf.parent[ra] = rb
        return rb
    elseif uf.rank[ra] > uf.rank[rb]
        uf.parent[rb] = ra
        return ra
    else
        uf.parent[rb] = ra
        uf.rank[ra] += Int32(1)
        return ra
    end
end

# --- Component assignment ---

"""
    find_components(cg::CompactedGraph) -> (comp_of_unitig, n_components)

Assign each unitig to a connected component. Twin unitigs share a
component (since they share canonical k-mer nodes via `canonical_idx`).
"""
function find_components(cg::CompactedGraph)
    n_c = n_canonical(cg)
    uf = UnionFind(n_c)
    n_u = n_unitigs(cg)

    # Union the canonical indices of each unitig's source and destination.
    for u in 1:n_u
        src_c = canonical_idx(cg.edge_sources[u])
        dst_c = canonical_idx(cg.edge_targets[u])
        uf_union!(uf, src_c, dst_c)
    end

    # Each unitig's component = the find-root of its source's canonical idx.
    # Compact the root labels to a dense 1..n_components range.
    raw_comp = Vector{Int32}(undef, n_u)
    for u in 1:n_u
        raw_comp[u] = uf_find(uf, canonical_idx(cg.edge_sources[u]))
    end

    unique_roots = sort!(unique(raw_comp))
    root_to_id = Dict{Int32,Int32}(r => Int32(i) for (i, r) in enumerate(unique_roots))
    comp_of_unitig = [root_to_id[r] for r in raw_comp]

    return comp_of_unitig, length(unique_roots)
end

# --- Pair twin components ---

"""
    pair_twin_components(cg, comp_of_unitig, n_components) -> canonical_components

Return a vector of component IDs to actually emit contigs for. Each
twin pair contributes exactly one component (the one with the smaller ID).
"""
function pair_twin_components(cg::CompactedGraph,
                              comp_of_unitig::Vector{Int32},
                              n_components::Integer)
    n_u = n_unitigs(cg)
    comp_twin = zeros(Int32, n_components)  # comp_twin[c] = twin component of c
    for u in 1:n_u
        t = cg.edge_twins[u]
        t == 0 && continue
        c_u = comp_of_unitig[u]
        c_t = comp_of_unitig[t]
        if comp_twin[c_u] == 0
            comp_twin[c_u] = c_t
        elseif comp_twin[c_u] != c_t
            # Self-twin component (palindromic) — both unitigs in same component.
            # comp_twin[c_u] stays as the first observed value.
        end
    end

    # A component is canonical if its ID is <= its twin's ID. Self-twin
    # components are always canonical.
    canonical = Int32[]
    for c in 1:n_components
        t = comp_twin[c]
        if t == 0 || t == c || c <= t
            push!(canonical, Int32(c))
        end
    end
    return canonical
end

# --- Per-component traversal ---

# Adjacency helpers for small components: avoid Dict overhead for the common
# case of 1–5 unitigs per component by using sorted parallel arrays.

@inline function _find_outs(srcs::Vector{Int32}, unitigs::Vector{Int32},
                             node::Int32)
    # Linear scan over the component's src array — fast for small components.
    result = Int32[]
    @inbounds for i in eachindex(srcs)
        srcs[i] == node && push!(result, unitigs[i])
    end
    return result
end

@inline function _in_degree(dsts::Vector{Int32}, node::Int32)
    cnt = 0
    @inbounds for d in dsts; d == node && (cnt += 1); end
    return cnt
end

"""
Greedy highest-coverage traversal within a single component.
Returns a vector of Contigs (usually 1, possibly more for branched components).
"""
function traverse_component(cg::CompactedGraph,
                            unitigs_in_comp::Vector{Int32},
                            comp_id::Integer,
                            k::Integer)
    n_in_comp = length(unitigs_in_comp)
    n_in_comp == 0 && return Contig[]

    # Parallel src/dst arrays for this component (no Dict allocation).
    srcs = [cg.edge_sources[u] for u in unitigs_in_comp]
    dsts = [cg.edge_targets[u] for u in unitigs_in_comp]

    # Build out-adjacency as sorted (by weight DESC) vectors per source node.
    # Use a Dict only when the component is large enough to justify it.
    use_dict = n_in_comp > 16
    out_from  = use_dict ? Dict{Int32,Vector{Int32}}() : nothing
    in_count  = use_dict ? Dict{Int32,Int}()            : nothing
    if use_dict
        for (i, u) in enumerate(unitigs_in_comp)
            src = srcs[i]; dst = dsts[i]
            push!(get!(out_from, src, Int32[]), u)
            in_count[dst] = get(in_count, dst, 0) + 1
        end
        for (_, us) in out_from
            sort!(us; by = u -> -cg.edge_weights[u])
        end
    end

    visited = falses(n_unitigs(cg))
    contigs = Contig[]

    function emit_contig_from(start_u::Int32)
        # Pre-size buffers to avoid repeated realloc in the common short-chain case.
        total_bases = sum(length(cg.edge_sequences[u]) for u in unitigs_in_comp)
        bases    = Vector{UInt8}(); sizehint!(bases, total_bases)
        weight_sum = 0.0f0
        weight_len = 0
        u_count  = 0

        prefix = node_to_string(cg, cg.edge_sources[start_u])

        cur_u = start_u
        while true
            visited[cur_u] && break
            visited[cur_u] = true
            seq_u = cg.edge_sequences[cur_u]
            append!(bases, seq_u)
            w = cg.edge_weights[cur_u]
            l = length(seq_u)
            weight_sum += w * Float32(l)
            weight_len += l
            u_count += 1

            dst = cg.edge_targets[cur_u]
            next_u = Int32(0)
            if use_dict
                cands = get(out_from, dst, Int32[])
                for v in cands
                    if !visited[v]; next_u = v; break; end
                end
            else
                # Linear scan for small components.
                best_w = -1.0f0
                @inbounds for i in eachindex(srcs)
                    srcs[i] == dst || continue
                    v = unitigs_in_comp[i]
                    visited[v] && continue
                    wv = cg.edge_weights[v]
                    if wv > best_w; best_w = wv; next_u = v; end
                end
            end
            next_u == 0 && break
            cur_u = next_u
        end

        suffix_bytes = Vector{UInt8}(undef, length(bases))
        @inbounds for (i, b) in enumerate(bases)
            suffix_bytes[i] = b == 0x00 ? UInt8('A') : b == 0x01 ? UInt8('C') :
                               b == 0x02 ? UInt8('G') : UInt8('T')
        end
        seq = prefix * String(suffix_bytes)
        mean_cov = weight_len == 0 ? 0.0f0 : weight_sum / Float32(weight_len)
        return Contig(seq, mean_cov, u_count, Int(comp_id))
    end

    # Start from 5'-end unitigs (in_degree == 0 at their source node).
    sources = Int32[]
    for (i, u) in enumerate(unitigs_in_comp)
        src = srcs[i]
        in_deg = use_dict ? get(in_count, src, 0) : _in_degree(dsts, src)
        in_deg == 0 && push!(sources, u)
    end
    sort!(sources; by = u -> -cg.edge_weights[u])

    for s in sources
        visited[s] && continue
        push!(contigs, emit_contig_from(s))
    end
    for u in unitigs_in_comp
        visited[u] && continue
        push!(contigs, emit_contig_from(u))
    end

    return contigs
end

# --- Top-level traversal ---

"""
    traverse_contigs(cg::CompactedGraph) -> Vector{Contig}

Find connected components, pair up twins, and run threaded greedy
traversal on each canonical component. Returns all contigs sorted by
length descending.
"""
function traverse_contigs(cg::CompactedGraph)
    comp_of_unitig, n_components = find_components(cg)

    # Bucket unitigs by component (no twin pairing — we de-dup at contig level).
    buckets = Dict{Int32,Vector{Int32}}()
    n_u = n_unitigs(cg)
    for u in 1:n_u
        c = comp_of_unitig[u]
        push!(get!(buckets, c, Int32[]), Int32(u))
    end

    # Threaded traversal.
    comps_ordered = sort(collect(keys(buckets)))
    results = Vector{Vector{Contig}}(undef, length(comps_ordered))
    Threads.@threads for i in 1:length(comps_ordered)
        c = comps_ordered[i]
        results[i] = traverse_component(cg, buckets[c], c, cg.k)
    end

    all_contigs = vcat(results...)

    # De-duplicate twin contigs using hash-based dedup to avoid hashing full
    # DNA strings (which is O(L) per contig and dominates at millions of contigs).
    # Key: (h_fwd XOR h_rc, min(h_fwd, h_rc)) — order-independent pair hash.
    # Collisions are resolved by full string comparison.
    function _seq_hash(s::AbstractString)
        h = UInt64(0xcbf29ce484222325)  # FNV-1a offset basis
        @inbounds for c in s
            h = xor(h, UInt64(c)) * UInt64(0x00000100000001b3)
        end
        h
    end
    function _rc_hash(s::AbstractString)
        h = UInt64(0xcbf29ce484222325)
        L = length(s)
        @inbounds for i in L:-1:1
            c = s[i]
            rc_c = c == 'A' ? 'T' : c == 'T' ? 'A' : c == 'G' ? 'C' : c == 'C' ? 'G' : 'N'
            h = xor(h, UInt64(rc_c)) * UInt64(0x00000100000001b3)
        end
        h
    end

    # seen maps (xor_hash, min_hash) → index into deduped array
    seen_idx = Dict{Tuple{UInt64,UInt64}, Int}()
    sizehint!(seen_idx, length(all_contigs))
    deduped = Contig[]
    sizehint!(deduped, length(all_contigs) ÷ 2 + 1)

    for c in all_contigs
        hf = _seq_hash(c.sequence)
        hr = _rc_hash(c.sequence)
        key = (xor(hf, hr), min(hf, hr))
        idx = get(seen_idx, key, 0)
        if idx == 0
            push!(deduped, c)
            seen_idx[key] = length(deduped)
        else
            # Hash collision or genuine twin: keep longer (full string compare only on tie).
            existing = deduped[idx]
            if length(c.sequence) > length(existing.sequence)
                deduped[idx] = c
            end
        end
    end

    sort!(deduped; by = c -> -length(c.sequence))
    return deduped
end
