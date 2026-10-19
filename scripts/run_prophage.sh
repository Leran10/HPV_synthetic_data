#!/bin/bash
#SBATCH --job-name=prophage
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/sahlab/leran/HPV_synthetic_data/data/bulk/logs/prophage_%j.log
#SBATCH --error=/scratch/sahlab/leran/HPV_synthetic_data/data/bulk/logs/prophage_%j.err

# prophage-pipeline driver for the jackalope BULK cohort.
#
# Usage:  sbatch run_prophage.sh <reads_dir> <outdir>
#
# This is the Snakemake DRIVER job: it stays resident and submits the real work
# as separate SLURM jobs via the pipeline's own profile (40 concurrent), so it
# needs little CPU/RAM itself but a long wall clock. Assembly + binning +
# GTDB-Tk across 100 samples is a multi-hour to multi-day job.
#
# Two environment requirements that are easy to miss:
#   1. Snakemake 8+ with the SLURM executor plugin. Hecatomb's env ships 7.18.2,
#      which cannot run this pipeline.
#   2. `conda` must be ON PATH, not merely installed. The profile sets
#      conda-prefix to the pre-built env directory, but Snakemake still shells
#      out to `conda info --json` and dies with exit 127 if conda is unreachable.
#
# Resumable: if this dies partway, fix the cause and resubmit -- completed rules
# are skipped. That is why we run all 100 directly rather than piloting first.

set -euo pipefail

READS=${1:?usage: sbatch run_prophage.sh <reads_dir> <outdir>}
OUTDIR=${2:?usage: sbatch run_prophage.sh <reads_dir> <outdir>}

PIPELINE=/scratch/sahlab/leran/HPV_synthetic_data/pipelines/prophage-pipeline
SMK=/ref/sahlab/software/miniforge3/envs/snakemake_v8.20.3/bin/snakemake

# conda on PATH -- see note 2 above
export PATH=/ref/sahlab/software/miniforge3/bin:/ref/sahlab/software/miniforge3/condabin:$PATH
command -v conda >/dev/null || { echo "FATAL: conda not on PATH" >&2; exit 1; }

mkdir -p "$OUTDIR" /scratch/sahlab/leran/HPV_synthetic_data/data/bulk/logs
cd "$PIPELINE/workflow"

n=$(ls "$READS"/*_R1.fastq.gz 2>/dev/null | wc -l)

echo "=================================================================="
echo " prophage-pipeline"
echo "   reads     : $READS  ($n samples)"
echo "   outdir    : $OUTDIR"
echo "   snakemake : $("$SMK" --version)"
echo "   fastp min : $(grep '^fastp_min_sequence_length' "$PIPELINE/config/config.yaml")"
echo "   tax scope : $(grep '^taxonomy_scope' "$PIPELINE/config/config.yaml")"
echo "   started   : $(date)"
echo "=================================================================="
echo

"$SMK" \
    --profile ../profile/slurm/ \
    --config reads="$READS" outdir="$OUTDIR" \
    --nolock

echo
echo " finished: $(date)"
echo "=== outputs ==="
ls "$OUTDIR" 2>/dev/null
