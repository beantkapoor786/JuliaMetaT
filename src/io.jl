# FASTA + TSV output writers, and pipeline log helpers.

"""
    open_log(out_dir) -> IO

Open (or create) a plain-text log file at `out_dir/log/pipeline.log` in append
mode. Writes a timestamped run header so multiple runs are separated cleanly.
"""
function open_log(out_dir::AbstractString)
    log_dir = joinpath(out_dir, "log")
    mkpath(log_dir)
    path = joinpath(log_dir, "pipeline.log")
    io = open(path, "a")
    println(io, "\n", "="^72)
    println(io, "RUN  ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
    println(io, "="^72)
    flush(io)
    return io
end

"""
    log_println(log_io, msg)

Print `msg` to stdout and to `log_io` (if not nothing), then flush.
"""
@inline function log_println(log_io, msg::AbstractString)
    println(msg)
    if log_io !== nothing
        println(log_io, msg)
        flush(log_io)
    end
end

"""
    write_contigs_fasta(path, contigs; min_length=0)

Write contigs to FASTA, one record per contig.
"""
function write_contigs_fasta(path::AbstractString,
                             contigs::AbstractVector{Contig};
                             min_length::Int = 0)
    open(FASTX.FASTA.Writer, path) do w
        for (i, c) in enumerate(contigs)
            length(c.sequence) >= min_length || continue
            id = "contig_$(i)"
            desc = "length=$(length(c.sequence)) coverage=$(round(c.mean_coverage, digits=2))"
            rec = FASTX.FASTA.Record("$id $desc", c.sequence)
            write(w, rec)
        end
    end
end

"""
    write_abundance_tsv(path, abundances; min_length=0)

Write per-contig abundance metrics to a tab-separated file.
Columns: contig_id, length, raw_count, rpkm, tpm, mean_coverage.
"""
function write_abundance_tsv(path::AbstractString,
                             abundances::AbstractVector{ContigAbundance};
                             min_length::Int = 0)
    open(path, "w") do io
        println(io, "contig_id\tlength\traw_count\trpkm\ttpm\tmean_coverage")
        for a in abundances
            a.contig_length >= min_length || continue
            println(io, "contig_$(a.contig_id)\t$(a.contig_length)\t",
                    "$(a.raw_count)\t",
                    "$(round(a.rpkm, digits=4))\t",
                    "$(round(a.tpm, digits=4))\t",
                    "$(round(a.mean_coverage, digits=2))")
        end
    end
end
