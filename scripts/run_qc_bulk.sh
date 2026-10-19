#!/bin/bash
#SBATCH --job-name=qc_bulk
#SBATCH --partition=general
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/sahlab/leran/HPV_synthetic_data/data/bulk/logs/qc_%j.log
#SBATCH --error=/scratch/sahlab/leran/HPV_synthetic_data/data/bulk/logs/qc_%j.err

# Hecatomb QC (preprocessing stage) for the BULK jackalope cohort.
# Snakemake driver job: stays resident and submits the real work via the
# slurm profile (up to 100 concurrent), so it needs little CPU/RAM itself.

set -euo pipefail

COHORT=/scratch/sahlab/leran/HPV_synthetic_data/data/bulk
CONFIG=$COHORT/hecatomb.qc.bulk.yaml
OUT=$COHORT/qc_hecatomb

source /ref/sahlab/software/anaconda3/bin/activate hecatomb_v1.1.0

mkdir -p "$OUT" "$COHORT/logs"
cd "$COHORT"

echo "=================================================================="
echo " Hecatomb QC -- jackalope BULK cohort (~40% host DNA)"
echo " started : $(date)"
echo " reads   : $COHORT/READS  ($(ls "$COHORT"/READS/*_R1.fastq.gz 2>/dev/null | wc -l) samples)"
echo " output  : $OUT"
echo "=================================================================="
echo

# --reads/--output/--preprocess/--host must be explicit: the CLI validates
# them before reading the config file.
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
echo " finished: $(date)"
ls -la "$OUT/results/" 2>/dev/null || echo "(no results directory)"
