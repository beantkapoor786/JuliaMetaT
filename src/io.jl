# FASTA + TSV output writers.

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
