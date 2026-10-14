#!/bin/bash
#SBATCH --job-name=phage_vlp
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/sahlab/leran/jackalope_cohort/logs/phage_vlp_%j.log
#SBATCH --error=/scratch/sahlab/leran/jackalope_cohort/logs/phage_vlp_%j.err

# phage-analysis on the VLP cohort.
#
# WHY VLP: the bulk cohort is ~1% viral, so its ~9,700 viral reads per sample
# were split across 25-54 genomes -- only 9/25 genomes reached 5x and 14 were
# below 1x. A 50 kb phage got ~1x and could not assemble. Result: 211 viral
# contigs from 1,287 true viral genomes, 3.1 kb mean contig, no cross-sample
# clustering. The VLP cohort is ~98% viral, which should lift a 50 kb phage from
# ~1x to ~110x.
#
# ALL RESOURCES ARE PINNED ON THE CLI. The profile's set-resources block is
# silently ignored by this Snakemake build -- partly because it writes runtime
# as time STRINGS ("24:00:00") where Snakemake wants integer MINUTES. Every
# failure so far came from a rule falling back to a default:
#   run 1: mmseqs got 4 CPUs (profile said 24)        -> 1h55m for 1 step
#   run 2: two rules hit the 1440-min default         -> TIMEOUT at 98%
#   run 3: iphop got 50 GB (profile said 100G)        -> OOM at 880/881
# So pin cpus/mem/runtime explicitly for every expensive rule, in minutes.
#
# VLP will produce MORE viral contigs than bulk (98% vs 1% viral), so the two
# mmseqs searches -- already the bottleneck at 12-24h each -- get 72h and the
# downstream per-split rules get generous ceilings.

set -euo pipefail

PIPELINE=/scratch/sahlab/leran/phage-analysis
SMK=/ref/sahlab/software/miniforge3/envs/snakemake_v8.20.3/bin/snakemake
CONFIG=${1:?usage: sbatch run_phage_analysis_vlp.sh <configfile>}

export PATH=/ref/sahlab/software/miniforge3/bin:/ref/sahlab/software/miniforge3/condabin:$PATH
command -v conda >/dev/null || { echo "FATAL: conda not on PATH" >&2; exit 1; }

mkdir -p /scratch/sahlab/leran/jackalope_cohort/logs
cd "$PIPELINE"

echo "=================================================================="
echo " phage-analysis -- VLP cohort"
echo "   config    : $CONFIG"
echo "   assembly  : $(grep '^assembly_file' "$CONFIG" | cut -d'"' -f2)"
echo "   outdir    : $(grep '^output_dir'   "$CONFIG" | cut -d'"' -f2)"
echo "   snakemake : $("$SMK" --version)"
echo "   started   : $(date)"
echo "=================================================================="
echo

# runtime is in MINUTES. 4320 = 72h, 1440 = 24h, 720 = 12h, 240 = 4h.
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
        genomad_prediction:cpus_per_task=16 \
        genomad_prediction:mem_mb=100000 \
        genomad_prediction:runtime=1440 \
        jaeger_prediction:cpus_per_task=16 \
        jaeger_prediction:mem_mb=60000 \
        jaeger_prediction:runtime=1440 \
        checkv_assessment:cpus_per_task=16 \
        checkv_assessment:mem_mb=60000 \
        checkv_assessment:runtime=720 \
        checkv_final_assessment:cpus_per_task=16 \
        checkv_final_assessment:mem_mb=60000 \
        checkv_final_assessment:runtime=720 \
        cluster_phages:cpus_per_task=16 \
        cluster_phages:mem_mb=60000 \
        cluster_phages:runtime=720 \
        prodigal_orf_prediction:cpus_per_task=16 \
        prodigal_orf_prediction:mem_mb=40000 \
        prodigal_orf_prediction:runtime=720 \
        phold_single_prediction:cpus_per_task=8 \
        phold_single_prediction:mem_mb=60000 \
        phold_single_prediction:runtime=1440 \
        phabox_prediction:cpus_per_task=16 \
        phabox_prediction:mem_mb=100000 \
        phabox_prediction:runtime=2880 \
        vcontact3_taxonomy:cpus_per_task=16 \
        vcontact3_taxonomy:mem_mb=200000 \
        vcontact3_taxonomy:runtime=2880 \
        bacphlip_lifestyle:mem_mb=40000 \
        bacphlip_lifestyle:runtime=720 \
        check_iphop_input_files:runtime=720 \
        iphop_single_prediction:cpus_per_task=16 \
        iphop_single_prediction:mem_mb=150000 \
        iphop_single_prediction:runtime=2880 \
        iphop_aggregate_results:mem_mb=32000 \
        iphop_aggregate_results:runtime=240 \
    --nolock

echo
echo " finished: $(date)"
