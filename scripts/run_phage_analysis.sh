#!/bin/bash
#SBATCH --job-name=phage_analysis
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/sahlab/leran/jackalope_bulk/logs/phage_%j.log
#SBATCH --error=/scratch/sahlab/leran/jackalope_bulk/logs/phage_%j.err

# phage-analysis driver for the jackalope BULK cohort.
#
# WHY THIS PIPELINE AND NOT prophage-pipeline:
# prophage-pipeline failed at 25% (job 48989998). It bins bacterial genomes and
# looks for prophages inside them, but ViroForge supplies bacterial content as
# short FRAGMENTS (bacterial_fragments.fasta) purely as realism ballast -- there
# are no whole genomes to bin. DAS Tool reported "No bins with bin-score >0.5
# found" on 230-306 kb assemblies, so filter_unbinned died on a missing
# unbinned.fa for 94 samples. Structural mismatch, not a tuning problem.
#
# phage-analysis works from assemblies and hunts viral contigs directly
# (Jaeger/geNomad/CheckV -> vOTU clustering -> taxonomy/lifestyle/host), with no
# binning step. That matches what ViroForge actually produces.
#
# INPUT: the 100 per-sample assemblies from the prophage run, pooled into one
# FASTA with sample-prefixed contig names (VIR_001__NODE_...). The pipeline takes
# a single assembly_file, and the prefix is what lets CoverM attribute abundance
# back to samples afterwards.
#
# Two environment requirements:
#   1. Snakemake 8+ with the SLURM executor plugin (Hecatomb's 7.18.2 cannot).
#   2. `conda` must be ON PATH -- Snakemake shells out to `conda info --json`
#      and exits 127 if it is unreachable, even though conda-prefix points at
#      the pre-built env directory.

set -euo pipefail

PIPELINE=/scratch/sahlab/leran/phage-analysis
SMK=/ref/sahlab/software/miniforge3/envs/snakemake_v8.20.3/bin/snakemake
CONFIG=${1:-$PIPELINE/config/my_config.yaml}

export PATH=/ref/sahlab/software/miniforge3/bin:/ref/sahlab/software/miniforge3/condabin:$PATH
command -v conda >/dev/null || { echo "FATAL: conda not on PATH" >&2; exit 1; }

mkdir -p /scratch/sahlab/leran/jackalope_bulk/logs
cd "$PIPELINE"

echo "=================================================================="
echo " phage-analysis"
echo "   config    : $CONFIG"
echo "   assembly  : $(grep '^assembly_file' "$CONFIG" | cut -d'"' -f2)"
echo "   outdir    : $(grep '^output_dir' "$CONFIG" | cut -d'"' -f2)"
echo "   snakemake : $("$SMK" --version)"
echo "   started   : $(date)"
echo "=================================================================="
echo

# The profile's `set-resources:` block is silently ignored by this Snakemake
# build, so the two mmseqs rules landed on the 4-CPU default despite declaring
# `threads: 24` and despite set-resources naming cpus_per_task. That made
# mmseqs_taxonomy the bottleneck (1h55m on step 4 of 47, no end in sight).
# CLI --set-resources DOES take effect, verified by --dry-run, so pin them here.
#
# Runtimes raised after run 49108760 reached 98% (866/881 steps) then lost two
# rules to TIMEOUT at the profile's 1440-minute default:
#   mmseqs_phage_taxonomy  - 72h. The first mmseqs search took 11h36m on 4780
#     contigs; this one has only 210 but adds -a and --tax-lineage 2 with a
#     22-field output, and the nr index is built with --splits 5 so every
#     search scans the database five times.
#   check_iphop_input_files - 12h. A file-existence test that should take
#     seconds; it timed out because it waits on the phold/iphop split files.
"$SMK" \
    --snakefile workflow/Snakefile \
    --profile profile/slurm/ \
    --configfile "$CONFIG" \
    --set-resources \
        mmseqs_taxonomy:cpus_per_task=24 \
        mmseqs_taxonomy:mem_mb=700000 \
        mmseqs_taxonomy:runtime=4320 \
        mmseqs_phage_taxonomy:cpus_per_task=24 \
        mmseqs_phage_taxonomy:mem_mb=700000 \
        mmseqs_phage_taxonomy:runtime=4320 \
        check_iphop_input_files:runtime=720 \
        iphop_single_prediction:runtime=1440 \
    --nolock

echo
echo " finished: $(date)"
