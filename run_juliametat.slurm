#!/bin/bash
#SBATCH --job-name=juliametat
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:2
#SBATCH --mem=256G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# --- Fill in these absolute paths ---
R1=/path/to/your_R1.fastq.gz
R2=/path/to/your_R2.fastq.gz
OUTDIR=/path/to/output
PROJDIR=/path/to/julia-metat-experimental

# --- Load modules (adjust to your cluster's module names) ---
module load julia/1.12
module load kmc/3
module load cuda

# --- Julia depot / package cache (avoids re-precompiling every run) ---
export JULIA_DEPOT_PATH=$PROJDIR/.julia_depot:$JULIA_DEPOT_PATH
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# JACC CUDA backend — JACC will use the first visible GPU automatically.
export JACC_BACKEND=cuda

mkdir -p "$OUTDIR/fasta" "$OUTDIR/abundance"

cd "$PROJDIR"

# Instantiate (no-op if depot is warm)
julia --project=. -e 'using Pkg; Pkg.instantiate()'

julia --project=. \
      -t $SLURM_CPUS_PER_TASK \
      -e "
using JuliaMetaT
run_pipeline(
    \"$R1\", \"$R2\";
    k                  = 31,
    min_count          = 2,
    min_edge_weight    = 2,
    relative_threshold = 0.05,
    min_hits           = 3,
    min_contig_length  = 100,
    out_fasta          = \"$OUTDIR/fasta/contigs.fasta\",
    out_tsv            = \"$OUTDIR/abundance/abundance.tsv\",
    verbose            = true,
)
"
