# JuliaMetaT

A de novo metatranscriptomic assembler for prokaryotic communities, written in Julia. It is built around a doubled-directed de Bruijn graph and uses [JACC.jl](https://github.com/JuliaParallel/JACC.jl) for backend-agnostic parallelism, so the same code runs on CPU threads, NVIDIA GPUs (CUDA), AMD GPUs (ROCm), and Apple Metal without modification.

The project targets real environmental metatranscriptomes at 2B+ read pairs, where reference genomes are unavailable and transcript abundance varies over a wide dynamic range.

---

## What it does

JuliaMetaT takes adapter-trimmed, ribodepleted paired-end FASTQ files and produces assembled contigs with per-contig abundance estimates.

The pipeline stages are:

1. **Load** reads into a packed 2-bit matrix
2. **Count k-mers** via KMC3 (external, fast, memory-bounded)
3. **Build a doubled-directed de Bruijn graph** from the filtered k-mer set
4. **Prune** low-coverage edges with a floor pass followed by a relative (per-node) pass that converges to a fixed point
5. **Compact** surviving paths into unitigs
6. **Traverse** unitig chains greedily by coverage to produce contigs
7. **Map** reads back to contigs via streaming Boyer-Moore majority vote over a sorted k-mer index
8. **Compute abundance** (RPKM and TPM) from the read assignments

The doubled-directed graph design (as opposed to a bidirected graph) eliminates a class of chimeric contig errors caused by ambiguous reverse-complement path merging. This is the single highest-leverage correctness property in the codebase.

---

## Requirements

- Julia 1.12+
- [KMC3](https://github.com/refresh-bio/KMC) on your `PATH`
- A JACC-compatible backend: CPU threads (default), CUDA, ROCm, or Metal

The Julia dependencies are declared in `Project.toml` and pinned in `Manifest.toml`. Running `Pkg.instantiate()` fetches everything.

---

## Quick start

```julia
using Pkg
Pkg.instantiate()

using JuliaMetaT

run_pipeline(
    "reads_R1.fastq.gz",
    "reads_R2.fastq.gz";
    k                  = 31,
    min_count          = 2,
    min_edge_weight    = 2,
    relative_threshold = 0.05,
    min_hits           = 3,
    min_contig_length  = 100,
    out_fasta          = "contigs.fasta",
    out_tsv            = "abundance.tsv",
    verbose            = true,
)
```

The input reads are expected to be **adapter-trimmed and ribodepleted**. The pipeline does not include a preprocessing stage. If your reads contain significant rRNA content, KMC3 will process them but the resulting graph will be dominated by rRNA k-mers and assembly quality will suffer.

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `k` | 31 | K-mer length for graph construction |
| `min_count` | 2 | Minimum k-mer frequency to include in the graph |
| `min_edge_weight` | 2 | Floor threshold for edge pruning |
| `relative_threshold` | 0.05 | Per-node relative pruning threshold (fraction of max outgoing weight) |
| `min_hits` | 3 | Minimum k-mer hits for a read to be assigned to a contig |
| `min_contig_length` | 100 | Minimum contig length to include in output |

---

## Running on a SLURM cluster

A ready-to-use SLURM script targeting GPU nodes is included at `run_juliametat.slurm`. It sets up the Julia depot, loads CUDA, and passes the correct thread count automatically. Fill in your account, paths, and module names at the top.

For GPU nodes, set `JACC_BACKEND=cuda`. For CPU-only nodes, set `JACC_BACKEND=threads` and drop the `--gres` and `module load cuda` lines.

On a first run, Julia will precompile packages into `$JULIA_DEPOT_PATH`. If compute nodes on your cluster do not have internet access, run `julia --project=. -e 'using Pkg; Pkg.instantiate()'` from a login node first to warm the depot before submitting.

---

## Backend selection

JACC selects a backend at runtime based on the `JACC_BACKEND` environment variable. If it is not set, JACC tries CUDA, then ROCm, then Metal, then falls back to threads.

```bash
export JACC_BACKEND=threads   # CPU, always available
export JACC_BACKEND=cuda      # NVIDIA GPUs
export JACC_BACKEND=rocm      # AMD GPUs
export JACC_BACKEND=metal     # Apple Silicon
```

The threads backend is the reference implementation. All correctness tests must pass on threads before a GPU backend is considered production-ready.

---

## Running tests

```bash
julia --project=. -t 4 test/runtests.jl
```

The test suite covers encoding, k-mer extraction, KMC3 agreement, graph construction and pruning, contig traversal, read mapping, and abundance output. It uses synthetic reads generated internally so no external data files are required.

---

## Performance at 1M read pairs (Apple M-series, warm run)

| Stage | Time |
|---|---|
| Load | 3.4s |
| KMC3 count + parse | 4.0s |
| Graph build + prune + compact | 5.4s |
| Contig traversal | 1.6s |
| Read mapping | 2.1s |
| **Total** | **~16.6s** |

Mapping rate on real ribodepleted data: ~91%. The graph and pruning stages use JACC kernels and scale to GPU backends for larger datasets. The scalability design target is 2B read pairs (4B total reads) on HPC nodes with 1 TB RAM; the primary engineering milestone before that scale is chunked GPU dispatch in the mapping stage (current implementation uploads the full read matrix to GPU VRAM in one shot).

---

## Project structure

```
src/
  pipeline.jl       top-level orchestrator
  graph.jl          de Bruijn graph construction and pruning
  traversal.jl      unitig compaction and contig traversal
  mapping.jl        read mapping and abundance computation
  kmc_backend.jl    KMC3 integration
  kmers.jl          in-memory k-mer counting (reference/test oracle)
  encoding.jl       2-bit base encoding and FASTQ loading
  io.jl             FASTA and TSV output
  backend.jl        JACC backend selection
  synthetic.jl      synthetic read generation for tests
test/
  runtests.jl
```
