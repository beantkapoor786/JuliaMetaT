using Test
using JuliaMetaT
using FASTX

@testset "Round 0: backend + encoding + synthetic" begin
    backend = JuliaMetaT.set_backend()
    @test backend in (:cuda, :amdgpu, :metal, :threads)
    @test JuliaMetaT.current_backend() == backend

    rec = FASTX.FASTQ.Record("r1", "ATGCAT", "IIIIII")
    seqs, lengths, n_dropped = JuliaMetaT.encode_reads([rec])
    @test size(seqs) == (6, 1)
    @test seqs[:, 1] == Int8[1, 2, 3, 4, 1, 2]
    @test lengths[1] == 6
    @test n_dropped == 0
    @test JuliaMetaT.decode_read(seqs, lengths, 1) == "ATGCAT"

    rec_ok  = FASTX.FASTQ.Record("ok",  "ATGC", "IIII")
    rec_bad = FASTX.FASTQ.Record("bad", "ATNC", "IIII")
    seqs2, lengths2, n_dropped2 = JuliaMetaT.encode_reads([rec_ok, rec_bad])
    @test size(seqs2, 2) == 1
    @test n_dropped2 == 1
    @test JuliaMetaT.decode_read(seqs2, lengths2, 1) == "ATGC"

    r1 = tempname() * "_R1.fastq"
    r2 = tempname() * "_R2.fastq"
    truth = JuliaMetaT.generate_synthetic_paired_fastq(r1, r2)
    @test length(truth.transcripts) == 3
    @test all(length.(truth.transcripts) .== 500)
    @test truth.read_len == 151
    @test truth.coverages == (10, 50, 200)
    @test 400 <= truth.n_pairs <= 432
    @test truth.n_reads_total == 2 * truth.n_pairs

    seqs3, lengths3, n_dropped3 = JuliaMetaT.load_fastq(r1, r2)
    @test size(seqs3, 2) == truth.n_reads_total
    @test all(lengths3 .== 151)
    @test n_dropped3 == 0

    rm(r1); rm(r2)
    
    # --- Gzip input support ---
    using CodecZlib
    r1 = tempname() * "_R1.fastq"
    r2 = tempname() * "_R2.fastq"
    truth = JuliaMetaT.generate_synthetic_paired_fastq(r1, r2)

    # Gzip-compress copies.
    r1gz = r1 * ".gz"
    r2gz = r2 * ".gz"
    for (src, dst) in ((r1, r1gz), (r2, r2gz))
        open(CodecZlib.GzipCompressorStream, dst, "w") do out
            write(out, read(src))
        end
    end

    # Load plain vs gzipped — should produce identical matrices.
    seqs_plain, lens_plain, _ = JuliaMetaT.load_fastq(r1, r2)
    seqs_gz, lens_gz, _       = JuliaMetaT.load_fastq(r1gz, r2gz)
    @test seqs_plain == seqs_gz
    @test lens_plain == lens_gz

    rm(r1); rm(r2); rm(r1gz); rm(r2gz)
end

@testset "Round 1: k-mer extraction & canonical counting" begin
    JuliaMetaT.set_backend()

    # --- Test A: hand-checkable case, k=4 ---
    seqs = reshape(Int8[1, 4, 3, 2, 1, 4, 3, 2], 8, 1)
    lens = Int32[8]
    flat = JuliaMetaT.extract_kmers(seqs, lens; k = 4)
    @test length(flat) == 5
    decoded = sort([JuliaMetaT.kmer_to_string(x, 4) for x in flat])
    @test decoded == sort(["ACGT", "CGTA", "GTAC", "CGTA", "ACGT"])

    uniq, cnts = JuliaMetaT.count_kmers(copy(flat); min_count = 1)
    pairs = Dict(JuliaMetaT.kmer_to_string(u, 4) => c
                 for (u, c) in zip(uniq, cnts))
    @test pairs["ACGT"] == 2
    @test pairs["CGTA"] == 2
    @test pairs["GTAC"] == 1

    uniq2, cnts2 = JuliaMetaT.count_kmers(copy(flat); min_count = 2)
    pairs2 = Dict(JuliaMetaT.kmer_to_string(u, 4) => c
                  for (u, c) in zip(uniq2, cnts2))
    @test haskey(pairs2, "ACGT") && pairs2["ACGT"] == 2
    @test haskey(pairs2, "CGTA") && pairs2["CGTA"] == 2
    @test !haskey(pairs2, "GTAC")

    # --- Test B: a k-mer and its revcomp must collide ---
    seqs_b = Int8[1 2; 1 2; 1 2; 1 2; 1 2]
    lens_b = Int32[5, 5]
    flat_b = JuliaMetaT.extract_kmers(seqs_b, lens_b; k = 5)
    @test length(flat_b) == 2
    @test all(x -> JuliaMetaT.kmer_to_string(x, 5) == "AAAAA", flat_b)
    uniq_b, cnts_b = JuliaMetaT.count_kmers(copy(flat_b); min_count = 1)
    @test length(uniq_b) == 1
    @test cnts_b[1] == 2

    # --- Test C: synthetic paired-end end-to-end at k=31 ---
    r1 = tempname() * "_R1.fastq"
    r2 = tempname() * "_R2.fastq"
    truth = JuliaMetaT.generate_synthetic_paired_fastq(r1, r2)
    seqs_h, lengths_h, _ = JuliaMetaT.load_fastq(r1, r2)

    flat_c = JuliaMetaT.extract_kmers(seqs_h, lengths_h; k = 31)
    @test length(flat_c) == 121 * size(seqs_h, 2)

    uniq_c, cnts_c = JuliaMetaT.count_kmers(copy(flat_c); min_count = 2)

    true_kmers = Set{UInt64}()
    for tx in truth.transcripts
        for i in 1:(length(tx) - 31 + 1)
            push!(true_kmers,
                  JuliaMetaT.string_to_canonical_kmer(tx[i:i+30], 31))
        end
    end
    @test length(true_kmers) >= 1400 && length(true_kmers) <= 1420

    recovered = count(k -> k in true_kmers, uniq_c)
    @test recovered >= 1300
    @test recovered / length(true_kmers) >= 0.92

    true_mask = [k in true_kmers for k in uniq_c]
    true_counts = cnts_c[true_mask]
    @test maximum(true_counts) >= 140
    @test maximum(true_counts) <= 250

    false_counts = cnts_c[.!true_mask]
    @test maximum(false_counts) <= 10

    rm(r1); rm(r2)
end

@testset "Round 2: doubled-directed de Bruijn graph" begin
    JuliaMetaT.set_backend()

    # --- Test A: small hand-checkable graph ---
    # Build a graph from a single short transcript and verify structure.
    tx = "AAACGTACGTAAAACGTACGT"   # 21bp, k=4
    k = 4
    tx_kmers = UInt64[]
    by_kmer = Dict{UInt64,Int32}()
    for i in 1:(length(tx)-k+1)
        km = JuliaMetaT.string_to_canonical_kmer(tx[i:i+k-1], k)
        by_kmer[km] = get(by_kmer, km, Int32(0)) + Int32(1)
    end
    uniq_a = collect(keys(by_kmer))
    cnts_a = collect(values(by_kmer))

    g_a = JuliaMetaT.build_graph(uniq_a, cnts_a; k = k)
    # Doubled graph: 2 oriented nodes per canonical (k-1)-mer.
    @test JuliaMetaT.n_nodes(g_a) == 2 * length(g_a.canonical_kmers)
    # Each k-mer -> 2 directed edges (one per orientation).
    @test length(g_a.edge_targets) == 2 * length(uniq_a)
    @test all(g_a.edge_alive)
    # CSR sanity.
    @test g_a.edge_offsets[1] == 1
    @test g_a.edge_offsets[end] == length(g_a.edge_targets) + 1
    # Twin invariant.
    @test all(g_a.edge_twins[g_a.edge_twins[e]] == e
              for e in 1:length(g_a.edge_targets))

    # --- Test B: twin & orientation helpers ---
    @test JuliaMetaT.twin(JuliaMetaT.forward_id(7)) == JuliaMetaT.reverse_id(7)
    @test JuliaMetaT.twin(JuliaMetaT.reverse_id(7)) == JuliaMetaT.forward_id(7)
    @test JuliaMetaT.is_forward(JuliaMetaT.forward_id(3)) == true
    @test JuliaMetaT.is_forward(JuliaMetaT.reverse_id(3)) == false
    @test JuliaMetaT.canonical_idx(JuliaMetaT.forward_id(5)) == 5
    @test JuliaMetaT.canonical_idx(JuliaMetaT.reverse_id(5)) == 5

    # --- Test C: end-to-end synthetic, real assembly ---
    r1 = tempname() * "_R1.fastq"
    r2 = tempname() * "_R2.fastq"
    truth = JuliaMetaT.generate_synthetic_paired_fastq(r1, r2)
    seqs_h, lengths_h, _ = JuliaMetaT.load_fastq(r1, r2)
    flat = JuliaMetaT.extract_kmers(seqs_h, lengths_h; k = 31)
    uniq, cnts = JuliaMetaT.count_kmers(copy(flat); min_count = 2)

    g = JuliaMetaT.build_graph(uniq, cnts; k = 31)
    n_edges_before = count(g.edge_alive)

    n_removed = JuliaMetaT.remove_tips!(g;
        min_edge_weight = 2, relative_threshold = 0.05)
    @test n_removed > 0
    @test count(g.edge_alive) == n_edges_before - n_removed
    # Pruning removes the vast majority of error edges.
    @test n_removed >= 6000

    # Twin invariant preserved post-pruning.
    @test all(g.edge_alive[e] == g.edge_alive[g.edge_twins[e]]
              for e in 1:length(g.edge_targets))

    # After pruning, interior nodes should be cleanly linear.
    N = JuliaMetaT.n_nodes(g)
    indeg  = zeros(Int, N); outdeg = zeros(Int, N)
    for v in 1:N
        for e in g.edge_offsets[v]:(g.edge_offsets[v+1]-1)
            g.edge_alive[e] || continue
            outdeg[v] += 1
            indeg[g.edge_targets[e]] += 1
        end
    end
    n_linear = count(v -> indeg[v] == 1 && outdeg[v] == 1, 1:N)
    @test n_linear >= 1500   # ~470 interior nodes × 2 orientations × ~2 strong transcripts

    # --- Compact ---
    cg = JuliaMetaT.compact_unitigs(g)
    @test cg.k == 31
    @test JuliaMetaT.n_unitigs(cg) >= 6   # at least the 3 transcripts × 2 orientations
    @test JuliaMetaT.n_unitigs(cg) <= 50  # but not pathologically fragmented

    # --- Verify reconstruction against truth ---
    lens_bases = length.(cg.edge_sequences)
    longest = maximum(lens_bases)
    @test longest >= 400    # the 200× transcript assembles nearly end-to-end

    # Reverse-complement helper for strings.
    function rc(s::AbstractString)
        out = Vector{Char}(undef, length(s))
        L = length(s)
        for (i, c) in enumerate(s)
            out[L - i + 1] =
                c == 'A' ? 'T' : c == 'T' ? 'A' :
                c == 'G' ? 'C' : c == 'C' ? 'G' : 'N'
        end
        String(out)
    end

    # Reconstruct the top 6 unitigs and check they each match a true
    # transcript (forward or rev-comp) as a substring.
    top6 = sortperm(lens_bases; rev = true)[1:min(6, JuliaMetaT.n_unitigs(cg))]
    matches = 0
    for u in top6
        seq = JuliaMetaT.unitig_sequence(cg, u)
        seq_rc = rc(seq)
        for tx in truth.transcripts
            if occursin(seq, tx) || occursin(seq_rc, tx)
                matches += 1
                break
            end
        end
    end
    @test matches == length(top6)

    # The TWO highest-coverage transcripts (50× and 200×) should each
    # be reconstructed nearly end-to-end in BOTH orientations.
    # That gives us 4 unitigs of length >= 400.
    n_full = count(>=(400), lens_bases)
    @test n_full >= 4

    # Coverage signal: the highest-weight unitig should correspond to
    # the 200× transcript.
    max_weight = maximum(cg.edge_weights)
    @test max_weight >= 100    # 200x transcript edges average ~150-180

    rm(r1); rm(r2)
end

@testset "Round 3: connected components + threaded traversal" begin
    JuliaMetaT.set_backend()

    # --- Build a graph through Round 2 first ---
    r1 = tempname() * "_R1.fastq"
    r2 = tempname() * "_R2.fastq"
    truth = JuliaMetaT.generate_synthetic_paired_fastq(r1, r2)
    seqs_h, lengths_h, _ = JuliaMetaT.load_fastq(r1, r2)
    flat = JuliaMetaT.extract_kmers(seqs_h, lengths_h; k = 31)
    uniq, cnts = JuliaMetaT.count_kmers(copy(flat); min_count = 2)
    g = JuliaMetaT.build_graph(uniq, cnts; k = 31)
    JuliaMetaT.remove_tips!(g; min_edge_weight = 2, relative_threshold = 0.05)
    cg = JuliaMetaT.compact_unitigs(g)

    # --- Test A: find_components ---
    comp_of_unitig, n_components = JuliaMetaT.find_components(cg)
    @test length(comp_of_unitig) == JuliaMetaT.n_unitigs(cg)
    @test minimum(comp_of_unitig) >= 1
    @test maximum(comp_of_unitig) == n_components
    # We expect at most ~6 components (3 transcripts × 2 orientations).
    # Could be fewer if any transcript is palindromic (very unlikely on random).
    @test 3 <= n_components <= 10

    # --- Test B: traverse_contigs returns sensible contigs ---
    contigs = JuliaMetaT.traverse_contigs(cg)
    @test length(contigs) >= 2     # at least the 200x and 50x got reconstructed
    @test length(contigs) <= 20    # not pathologically fragmented

    # The longest two contigs should each be near full transcript length.
    @test length(contigs[1].sequence) >= 400
    @test length(contigs[2].sequence) >= 400

    # --- Test C: contigs match real transcripts ---
    function rc(s::AbstractString)
        out = Vector{Char}(undef, length(s))
        L = length(s)
        for (i, c) in enumerate(s)
            out[L - i + 1] =
                c == 'A' ? 'T' : c == 'T' ? 'A' :
                c == 'G' ? 'C' : c == 'C' ? 'G' : 'N'
        end
        String(out)
    end

    # Each contig should be a substring (in either orientation) of some
    # true transcript.
    matched = 0
    for c in contigs[1:min(2, end)]
        seq = c.sequence
        seq_rc = rc(seq)
        for tx in truth.transcripts
            if occursin(seq, tx) || occursin(seq_rc, tx)
                matched += 1
                break
            end
        end
    end
    @test matched == min(2, length(contigs))

    # --- Test D: coverage ordering ---
    # Among the long contigs (>=400bp), the one with highest coverage
    # should match the 200× transcript.
    long_contigs = filter(c -> length(c.sequence) >= 400, contigs)
    @test length(long_contigs) >= 2
    sort!(long_contigs; by = c -> -c.mean_coverage)
    @test long_contigs[1].mean_coverage >= 100   # 200× transcript
    @test long_contigs[1].mean_coverage >= long_contigs[2].mean_coverage

    # The top contig should match transcript 3 (the 200× one).
    top_seq = long_contigs[1].sequence
    @test occursin(top_seq, truth.transcripts[3]) ||
          occursin(rc(top_seq), truth.transcripts[3])

    # --- Test E: twin de-duplication (long contigs only) ---
    # Short fragments may be palindromic by chance; that's not a bug.
    # But no two LONG contigs should be reverse-complements of each other.
    long_only = filter(c -> length(c.sequence) >= 100, contigs)
    for i in 1:length(long_only), j in (i+1):length(long_only)
        si, sj = long_only[i].sequence, long_only[j].sequence
        abs(length(si) - length(sj)) > 5 && continue
        @test si != sj
        @test si != rc(sj)
    end
end

@testset "Round 4: read mapping + abundance + output" begin
    JuliaMetaT.set_backend()

    # Build contigs through Round 3 first.
    r1 = tempname() * "_R1.fastq"
    r2 = tempname() * "_R2.fastq"
    truth = JuliaMetaT.generate_synthetic_paired_fastq(r1, r2)
    seqs_h, lengths_h, _ = JuliaMetaT.load_fastq(r1, r2)
    flat = JuliaMetaT.extract_kmers(seqs_h, lengths_h; k=31)
    uniq, cnts = JuliaMetaT.count_kmers(copy(flat); min_count=2)
    g = JuliaMetaT.build_graph(uniq, cnts; k=31)
    JuliaMetaT.remove_tips!(g; min_edge_weight=2, relative_threshold=0.05)
    cg = JuliaMetaT.compact_unitigs(g)
    contigs = JuliaMetaT.traverse_contigs(cg)
    contigs = filter(c -> length(c.sequence) >= 100, contigs)
    @test length(contigs) >= 2

    # --- Test A: k-mer index ---
    sorted_kmers, contig_of_kmer = JuliaMetaT.build_contig_kmer_index(contigs, 31)
    @test length(sorted_kmers) == length(contig_of_kmer)
    @test issorted(sorted_kmers)
    @test allunique(sorted_kmers)
    @test all(1 .<= contig_of_kmer .<= length(contigs))

    # --- Test B: read mapping ---
    assignments, hits = JuliaMetaT.map_reads(contigs, seqs_h, lengths_h;
                                             k=31, min_hits=3)
    @test length(assignments) == size(seqs_h, 2)
    @test all(0 .<= assignments .<= length(contigs))
    # Most reads should map: 200x and 50x transcripts assembled cleanly.
    n_mapped = count(>(0), assignments)
    n_reads = size(seqs_h, 2)
    @test n_mapped / n_reads >= 0.8

    # --- Test C: abundance ---
    abundances = JuliaMetaT.compute_abundance(contigs, assignments)
    @test length(abundances) == length(contigs)
    @test sum(a.raw_count for a in abundances) == n_mapped

    # TPMs should sum to ~1e6 (off slightly due to integer counts).
    tpm_sum = sum(a.tpm for a in abundances)
    @test 999990 < tpm_sum < 1000010

    # Coverage ordering: the highest-coverage contig (200x transcript)
    # should have the most reads.
    sort!(abundances; by = a -> -a.raw_count)
    top = abundances[1]
    @test top.mean_coverage >= 100   # 200x transcript

    # Ratio sanity: 200x should be roughly 4x the 50x contig.
    # Find the two long, high-coverage abundances.
    long_abundances = filter(a -> a.contig_length >= 400, abundances)
    sort!(long_abundances; by = a -> -a.mean_coverage)
    if length(long_abundances) >= 2
        ratio = long_abundances[1].raw_count / max(long_abundances[2].raw_count, 1)
        @test 2.5 <= ratio <= 6.0   # expected ~4 (200/50)
    end

    # --- Test D: I/O writes valid files ---
    fasta_out = tempname() * ".fasta"
    tsv_out   = tempname() * ".tsv"
    JuliaMetaT.write_contigs_fasta(fasta_out, contigs; min_length=100)
    JuliaMetaT.write_abundance_tsv(tsv_out, abundances; min_length=100)

    @test isfile(fasta_out) && filesize(fasta_out) > 0
    @test isfile(tsv_out)   && filesize(tsv_out) > 0

    # Read FASTA back and check.
    written_contigs = open(FASTX.FASTA.Reader, fasta_out) do r
        collect(r)
    end
    @test length(written_contigs) == length(contigs)

    # TSV has correct header.
    tsv_lines = readlines(tsv_out)
    @test startswith(tsv_lines[1], "contig_id\tlength\traw_count")
    @test length(tsv_lines) == length(contigs) + 1   # header + 1 per contig

    rm(fasta_out); rm(tsv_out)

    # --- Test E: full pipeline entry point ---
    fasta2 = tempname() * ".fasta"
    tsv2   = tempname() * ".tsv"
    contigs_p, abundances_p = JuliaMetaT.run_pipeline(r1, r2;
        out_fasta = fasta2, out_tsv = tsv2, verbose = false)
    @test length(contigs_p) >= 2
    @test length(abundances_p) == length(contigs_p)
    @test isfile(fasta2) && isfile(tsv2)

    rm(fasta2); rm(tsv2); rm(r1); rm(r2)
end
