# Top-level orchestrator: runs all 5 stages, emits FASTA + TSV + log.

"""
    run_pipeline(r1_path, r2_path=nothing; k=31, min_count=2, ...,
                 out_fasta="contigs.fasta", out_tsv="abundance.tsv",
                 out_dir="output",
                 reads_per_batch=1_000_000, min_contig_length=100, verbose=true)

End-to-end metatranscriptomics assembly. Reads FASTQ (paired-end if
both paths given, single-end if only r1), runs the full pipeline,
writes FASTA + TSV outputs. Returns (contigs, abundances).

A log file is always written to `out_dir/log/pipeline.log` in append
mode with a timestamped header. Stage timings, counts, and any errors
(including stacktraces) are captured there. Reads are never fully
materialised in RAM — mapping streams FASTQ in `reads_per_batch`-sized
chunks so peak host memory is O(reads_per_batch), not O(total_reads).
"""
function run_pipeline(r1_path::AbstractString,
                      r2_path::Union{AbstractString,Nothing} = nothing;
                      k::Int = 31,
                      min_count::Int = 2,
                      min_edge_weight::Int = 2,
                      relative_threshold::Float64 = 0.05,
                      min_hits::Int = 3,
                      reads_per_batch::Int = 1_000_000,
                      out_fasta::AbstractString = "contigs.fasta",
                      out_tsv::AbstractString   = "abundance.tsv",
                      out_dir::AbstractString = "output",
                      min_contig_length::Int = 100,
                      verbose::Bool = true)
    set_backend()

    log_io = open_log(out_dir)

    try
        _run_pipeline_inner(r1_path, r2_path;
            k, min_count, min_edge_weight, relative_threshold,
            min_hits, reads_per_batch, out_fasta, out_tsv, min_contig_length,
            verbose, log_io)
    catch err
        msg = sprint(showerror, err, catch_backtrace())
        println(log_io, "\nERROR\n", msg); flush(log_io)
        rethrow()
    finally
        close(log_io)
    end
end

function _run_pipeline_inner(r1_path, r2_path;
                             k, min_count, min_edge_weight, relative_threshold,
                             min_hits, reads_per_batch, out_fasta, out_tsv,
                             min_contig_length, verbose, log_io)
    emit(msg) = verbose && log_println(log_io, msg)

    emit("r1=$(r1_path)")
    r2_path !== nothing && emit("r2=$(r2_path)")
    emit("params: k=$k  min_count=$min_count  min_edge_weight=$min_edge_weight  relative_threshold=$relative_threshold  min_hits=$min_hits  min_contig_length=$min_contig_length")

    paths = r2_path === nothing ? (r1_path,) : (r1_path, r2_path)

    t_kmer = @elapsed begin
        uniq, cnts = count_kmers_kmc(r1_path, r2_path;
            k = k, min_count = min_count, verbose = verbose, log_io = log_io)
    end
    emit("[kmers]   unique≥$(min_count)=$(length(uniq))  ($(round(t_kmer, digits=3))s)")

    t_graph = @elapsed begin
        g = build_graph(uniq, cnts; k = k)
        n_removed = remove_tips!(g; min_edge_weight, relative_threshold, verbose, log_io)
        cg = compact_unitigs(g; verbose, log_io)
    end
    emit("[graph]   edges_pruned=$n_removed unitigs=$(n_unitigs(cg))  ($(round(t_graph, digits=3))s)")

    t_traverse = @elapsed begin
        contigs = traverse_contigs(cg)
        contigs = filter(c -> length(c.sequence) >= min_contig_length, contigs)
    end
    emit("[contigs] n=$(length(contigs)) longest=$(isempty(contigs) ? 0 : length(contigs[1].sequence))  ($(round(t_traverse, digits=3))s)")

    if isempty(contigs)
        emit("[done]    no contigs >= $min_contig_length bp; skipping mapping")
        return contigs, ContigAbundance[]
    end

    t_map = @elapsed begin
        result    = map_reads_streaming(contigs, paths...; k = k, min_hits = min_hits,
                                        reads_per_batch = reads_per_batch,
                                        verbose = verbose, log_io = log_io)
        abundances = compute_abundance(contigs, result.raw_counts)
    end
    emit("[mapping] reads_mapped=$(result.n_mapped)/$(result.n_reads) ($(round(100*result.n_mapped/max(result.n_reads,1), digits=1))%) dropped=$(result.n_dropped)  ($(round(t_map, digits=3))s)")

    write_contigs_fasta(out_fasta, contigs; min_length = min_contig_length)
    write_abundance_tsv(out_tsv, abundances; min_length = min_contig_length)
    emit("[output]  $out_fasta  $out_tsv")

    return contigs, abundances
end
