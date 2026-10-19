#!/bin/bash
#SBATCH --job-name=stage_reads
#SBATCH --partition=general
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/sahlab/leran/HPV_synthetic_data/data/vlp/logs/stage_%j.log
#SBATCH --error=/scratch/sahlab/leran/HPV_synthetic_data/data/vlp/logs/stage_%j.err

# Stage the jackalope cohort FASTQs into a single directory with unique,
# metadata-joinable sample names (VIR_001_R1.fastq.gz ...).
#
# The source files are uncompressed and share filenames across samples
# (all HPV- samples are vaginal_virome_healthy_R1.fastq), which breaks any
# pipeline that infers sample identity from the filename. This gzips them
# into READS/ under their VIR_### id.

set -euo pipefail

COHORT=/scratch/sahlab/leran/HPV_synthetic_data/data/vlp
READS=$COHORT/READS

mkdir -p "$READS" "$COHORT/logs"

gzip_one() {
    local src=$1 dest=$2
    [ -s "$dest" ] && { echo "SKIP  $(basename "$dest") (exists)"; return 0; }
    gzip -c "$src" > "$dest.tmp" && mv "$dest.tmp" "$dest"
    echo "DONE  $(basename "$dest")"
}
export -f gzip_one

# Build the work list: one line per file as "source<TAB>destination"
: > /tmp/stage_list.$$
for d in "$COHORT"/VIR_*/; do
    s=$(basename "$d")
    for r in R1 R2; do
        src=$(ls "$d"/fastq/*_"$r".fastq 2>/dev/null | head -1)
        [ -z "$src" ] && { echo "WARN: no $r for $s" >&2; continue; }
        printf '%s\t%s\n' "$src" "$READS/${s}_${r}.fastq.gz" >> /tmp/stage_list.$$
    done
done

n=$(wc -l < /tmp/stage_list.$$)
echo "Staging $n files with ${SLURM_CPUS_PER_TASK:-8} parallel workers"
echo "Destination: $READS"
echo

# gzip in parallel, one file per core
xargs -a /tmp/stage_list.$$ -P "${SLURM_CPUS_PER_TASK:-8}" -I{} \
    bash -c 'IFS=$'"'"'\t'"'"' read -r src dest <<< "{}"; gzip_one "$src" "$dest"'

rm -f /tmp/stage_list.$$

echo
echo "=== staged files: $(ls "$READS" | wc -l) (expect 200) ==="
du -sh "$READS"
