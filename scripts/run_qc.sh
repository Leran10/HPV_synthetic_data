#!/bin/bash
#SBATCH --job-name=hecatomb_qc
#SBATCH --partition=general
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=${PROJECT_DIR}/data/vlp/logs/qc_%j.log
#SBATCH --error=${PROJECT_DIR}/data/vlp/logs/qc_%j.err

# Hecatomb QC (preprocessing stage only) for the jackalope cohort.
#
# This is the Snakemake DRIVER job -- it stays resident and submits the real
# work as separate SLURM jobs via the profile (up to 100 concurrent). So it
# needs little CPU/RAM itself but a long wall clock.
#
# Stage `preprocessing` stops after results/seqtable.fasta. It does NOT run the
# taxonomic cascade, assembly, or annotation -- those would take days and the
# phage pipelines do that work instead.

set -euo pipefail

COHORT=${PROJECT_DIR}/data/vlp
CONFIG=$COHORT/hecatomb.qc.yaml
OUT=$COHORT/qc_hecatomb

source ${CONDA_BASE}/bin/activate hecatomb_v1.1.0

mkdir -p "$OUT" "$COHORT/logs"
cd "$COHORT"

echo "=================================================================="
echo " Hecatomb QC -- jackalope cohort"
echo " started : $(date)"
echo " reads   : $COHORT/READS  ($(ls "$COHORT"/READS/*_R1.fastq.gz | wc -l) samples)"
echo " output  : $OUT"
echo " stage   : preprocessing"
echo "=================================================================="
echo

# --reads/--output/--preprocess/--host must be passed explicitly: the CLI
# validates them before reading the config file.
hecatomb run \
    --reads "$COHORT/READS" \
    --output "$OUT" \
    --preprocess paired \
    --host human \
    --search fast \
    --configfile "$CONFIG" \
    --profile slurm \
    --threads 16 \
    --snake-default="--rerun-incomplete --printshellcmds --nolock --show-failed-logs" \
    preprocessing

echo
echo "=================================================================="
echo " finished: $(date)"
echo "=================================================================="

echo
echo "=== results/ ==="
ls -la "$OUT/results/" 2>/dev/null || echo "(no results directory)"
