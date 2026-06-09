module JuliaMetaT

using FASTX
using JACC
using Random
using CodecZlib

include("backend.jl")
include("encoding.jl")
include("synthetic.jl")
include("kmers.jl")
include("graph.jl")
include("traversal.jl")
include("mapping.jl")
include("io.jl")
include("pipeline.jl")
include("kmc_backend.jl")

# Round 0
export set_backend, current_backend
export encode_reads, decode_read, load_fastq, stream_fastq
export generate_synthetic_paired_fastq

# Round 1
export extract_kmers, count_kmers, kmer_to_string, string_to_canonical_kmer

# KMC backend (Round 1, alternative path)
export count_kmers_kmc, check_kmc

# Round 2
export DeBruijnGraph, CompactedGraph
export build_graph, remove_tips!, compact_unitigs
export n_nodes, out_edges, out_degree
export n_canonical, n_unitigs
export node_to_string, unitig_sequence
export twin, forward_id, reverse_id, is_forward, canonical_idx

# Round 3
export Contig, find_components, traverse_contigs

# Round 4
export ContigAbundance
export build_contig_kmer_index, map_reads, compute_abundance
export write_contigs_fasta, write_abundance_tsv

# Pipeline
export run_pipeline

end # module
