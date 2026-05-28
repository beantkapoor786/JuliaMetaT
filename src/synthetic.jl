"""
    generate_synthetic_paired_fastq(r1_path, r2_path; ...)

Generate 3 random transcripts and simulate Illumina-style paired-end
reads. Writes R1 and R2 FASTQ files. Returns truth data.
"""
function generate_synthetic_paired_fastq(r1_path::AbstractString,
                                         r2_path::AbstractString;
                                         seed::Int = 42,
                                         transcript_len::Int = 500,
                                         coverages::NTuple{3,Int} = (10, 50, 200),
                                         read_len::Int = 151,
                                         insert_mean::Int = 300,
                                         insert_sd::Int = 50,
                                         error_rate::Float64 = 0.01)
    rng = Random.MersenneTwister(seed)
    bases = ('A', 'T', 'G', 'C')

    transcripts = [String(rand(rng, bases, transcript_len)) for _ in 1:3]

    r1_records = FASTX.FASTQ.Record[]
    r2_records = FASTX.FASTQ.Record[]
    pair_id = 0
    qual_str = String(fill('I', read_len))

    for (t_idx, (tx, cov)) in enumerate(zip(transcripts, coverages))
        n_pairs = cld(transcript_len * cov, 2 * read_len)
        for _ in 1:n_pairs
            insert = round(Int, insert_mean + insert_sd * randn(rng))
            insert = clamp(insert, 2 * read_len, transcript_len)
            max_start = transcript_len - insert + 1
            max_start < 1 && continue
            start = rand(rng, 1:max_start)

            r1_seq = collect(tx[start : start + read_len - 1])
            _add_errors!(r1_seq, rng, bases, error_rate)

            r2_src = tx[start + insert - read_len : start + insert - 1]
            r2_seq = collect(_revcomp_string(r2_src))
            _add_errors!(r2_seq, rng, bases, error_rate)

            pair_id += 1
            id1 = "read_$(pair_id)_tx$(t_idx)/1"
            id2 = "read_$(pair_id)_tx$(t_idx)/2"
            push!(r1_records, FASTX.FASTQ.Record(id1, String(r1_seq), qual_str))
            push!(r2_records, FASTX.FASTQ.Record(id2, String(r2_seq), qual_str))
        end
    end

    perm = randperm(rng, length(r1_records))
    r1_records = r1_records[perm]
    r2_records = r2_records[perm]

    open(FASTX.FASTQ.Writer, r1_path) do w
        for r in r1_records; write(w, r); end
    end
    open(FASTX.FASTQ.Writer, r2_path) do w
        for r in r2_records; write(w, r); end
    end

    return (transcripts = transcripts,
            n_pairs = length(r1_records),
            n_reads_total = 2 * length(r1_records),
            coverages = coverages,
            read_len = read_len)
end

function _add_errors!(chars::Vector{Char}, rng, bases, rate)
    @inbounds for i in eachindex(chars)
        if rand(rng) < rate
            orig = chars[i]
            chars[i] = rand(rng, filter(b -> b != orig, collect(bases)))
        end
    end
    return chars
end

function _revcomp_string(s::AbstractString)
    out = Vector{Char}(undef, length(s))
    L = length(s)
    @inbounds for (i, c) in enumerate(s)
        out[L - i + 1] =
            c == 'A' ? 'T' :
            c == 'T' ? 'A' :
            c == 'G' ? 'C' :
            c == 'C' ? 'G' :
                       'N'
    end
    return String(out)
end
