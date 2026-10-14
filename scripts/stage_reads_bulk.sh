#!/bin/bash
#SBATCH --job-name=stage_bulk
#SBATCH --partition=general
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/scratch/sahlab/leran/jackalope_bulk/logs/stage_%j.log
#SBATCH --error=/scratch/sahlab/leran/jackalope_bulk/logs/stage_%j.err

# Stage the BULK jackalope cohort into one directory with unique,
# metadata-joinable sample names (VIR_001_R1.fastq.gz ...).
#
# The generated files are uncompressed and share filenames across samples
# (every HPV- sample is vaginal_virome_healthy_R1.fastq), which breaks any
# pipeline that infers sample identity from the filename.

set -euo pipefail

COHORT=/scratch/sahlab/leran/jackalope_bulk
READS=$COHORT/READS

mkdir -p "$READS" "$COHORT/logs"

gzip_one() {
    local src=$1 dest=$2
    [ -s "$dest" ] && { echo "SKIP  $(basename "$dest")"; return 0; }
    gzip -c "$src" > "$dest.tmp" && mv "$dest.tmp" "$dest"
    echo "DONE  $(basename "$dest")"
}
export -f gzip_one

: > /tmp/stage_bulk.$$
for d in "$COHORT"/VIR_*/; do
    s=$(basename "$d")
    for r in R1 R2; do
        src=$(ls "$d"/fastq/*_"$r".fastq 2>/dev/null | head -1)
        [ -z "$src" ] && { echo "WARN: no $r for $s" >&2; continue; }
        printf '%s\t%s\n' "$src" "$READS/${s}_${r}.fastq.gz" >> /tmp/stage_bulk.$$
    done
done

n=$(wc -l < /tmp/stage_bulk.$$)
echo "Staging $n files with ${SLURM_CPUS_PER_TASK:-8} workers -> $READS"
echo

xargs -a /tmp/stage_bulk.$$ -P "${SLURM_CPUS_PER_TASK:-8}" -I{} \
    bash -c 'IFS=$'"'"'\t'"'"' read -r src dest <<< "{}"; gzip_one "$src" "$dest"'

rm -f /tmp/stage_bulk.$$

echo
echo "=== staged: $(ls "$READS" | wc -l) files (expect 200) ==="
du -sh "$READS"
