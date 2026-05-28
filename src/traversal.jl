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

    # Build local adjacency: out-unitigs from each oriented node within this component.
    # Use a Dict because node IDs can be sparse.
    out_from = Dict{Int32,Vector{Int32}}()   # node -> list of unitig indices starting there
    in_count = Dict{Int32,Int}()             # node -> number of unitigs ending there
    for u in unitigs_in_comp
        src = cg.edge_sources[u]
        dst = cg.edge_targets[u]
        push!(get!(out_from, src, Int32[]), u)
        in_count[dst] = get(in_count, dst, 0) + 1
    end

    # Sort each adjacency list by weight DESC, so popping picks highest first.
    for (_, us) in out_from
        sort!(us; by = u -> -cg.edge_weights[u])
    end

    visited = Set{Int32}()
    contigs = Contig[]

    function emit_contig_from(start_u::Int32)
        bases    = UInt8[]
        weights  = Float32[]
        u_count  = 0
        chain    = Int32[]

        # Source prefix: the (k-1)-mer at start_u's source node.
        prefix = node_to_string(cg, cg.edge_sources[start_u])

        cur_u = start_u
        while true
            (cur_u in visited) && break
            push!(visited, cur_u)
            push!(chain, cur_u)
            append!(bases, cg.edge_sequences[cur_u])
            push!(weights, cg.edge_weights[cur_u])
            u_count += 1

            # Move to next unitig: pick highest-weight unvisited out from this dst.
            dst = cg.edge_targets[cur_u]
            cands = get(out_from, dst, Int32[])
            next_u = Int32(0)
            for v in cands
                if !(v in visited)
                    next_u = v; break
                end
            end
            next_u == 0 && break
            cur_u = next_u
        end

        # Compose contig sequence: prefix + decoded bases.
        suffix_chars = Vector{Char}(undef, length(bases))
        @inbounds for (i, b) in enumerate(bases)
            suffix_chars[i] = b == 0 ? 'A' : b == 1 ? 'C' :
                              b == 2 ? 'G' : 'T'
        end
        seq = prefix * String(suffix_chars)

        # Coverage: weighted mean by unitig length (longer unitigs dominate).
        total_w = 0.0f0
        total_len = 0
        for (i, u) in enumerate(chain)
            ulen = length(cg.edge_sequences[u])
            total_w += weights[i] * Float32(ulen)
            total_len += ulen
        end
        mean_cov = total_len == 0 ? 0.0f0 : total_w / Float32(total_len)

        return Contig(seq, mean_cov, u_count, Int(comp_id))
    end

    # Pick starting unitigs: ones whose source has in_count == 0 (5' ends).
    # If none, fall back to any unvisited unitig.
    sources = Int32[]
    for u in unitigs_in_comp
        src = cg.edge_sources[u]
        if get(in_count, src, 0) == 0
            push!(sources, u)
        end
    end

    # Sort sources by descending weight so we emit the strongest contig first.
    sort!(sources; by = u -> -cg.edge_weights[u])

    for s in sources
        s in visited && continue
        push!(contigs, emit_contig_from(s))
    end

    # Any leftover unvisited unitigs (cycles or branches we missed) → emit too.
    for u in unitigs_in_comp
        u in visited && continue
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

    # De-duplicate twin contigs: for each contig, compute its canonical
    # form (lex-min of forward vs rev-comp). Keep one per canonical form,
    # preferring the longer contig (and the one whose sequence == canonical
    # form if lengths tie, for determinism).
    function _rc(s::AbstractString)
        out = Vector{Char}(undef, length(s))
        L = length(s)
        @inbounds for (i, c) in enumerate(s)
            out[L - i + 1] =
                c == 'A' ? 'T' : c == 'T' ? 'A' :
                c == 'G' ? 'C' : c == 'C' ? 'G' : 'N'
        end
        String(out)
    end

    seen = Dict{String,Contig}()
    for c in all_contigs
        rev = _rc(c.sequence)
        key = c.sequence <= rev ? c.sequence : rev
        if !haskey(seen, key) ||
           length(c.sequence) > length(seen[key].sequence)
            seen[key] = c
        end
    end

    deduped = collect(values(seen))
    sort!(deduped; by = c -> -length(c.sequence))
    return deduped
end
