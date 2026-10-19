#!/bin/bash
#SBATCH --job-name=asm_vlp
#SBATCH --partition=general
#SBATCH --array=1-100%25
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/scratch/sahlab/leran/jackalope_cohort/logs/asm_%A_%a.log
#SBATCH --error=/scratch/sahlab/leran/jackalope_cohort/logs/asm_%A_%a.err

# Assemble the VLP cohort. The bulk cohort's assemblies came free from the
# prophage-pipeline run (which died later, at binning); the VLP cohort has none,
# so assemble directly.
#
# Same MEGAHIT settings the prophage pipeline used, so the two cohorts'
# assemblies are comparable: --k-min 21 --k-max 141 --k-step 12, min contig
# 1000bp. Median runtime there was 226 s/sample (max 732 s) on bulk data; VLP is
# ~98% viral rather than ~1%, so expect more real signal and somewhat longer runs.

set -euo pipefail

C=/scratch/sahlab/leran/jackalope_cohort
READS=$C/READS
OUT=$C/assemblies

mkdir -p "$OUT" "$C/logs"

i=${SLURM_ARRAY_TASK_ID}
s=$(printf "VIR_%03d" "$i")

r1=$READS/${s}_R1.fastq.gz
r2=$READS/${s}_R2.fastq.gz
[ -s "$r1" ] && [ -s "$r2" ] || { echo "FATAL: missing reads for $s" >&2; exit 1; }

final=$OUT/${s}.contigs_1000bp.fasta
[ -s "$final" ] && { echo "SKIP: $s already assembled"; exit 0; }

export PATH=/ref/sahlab/software/miniforge3/envs/smk_envs_prophage_pipeline/23e6050a730b5b5a42945305b9ffed01_/bin:$PATH
command -v megahit >/dev/null || {
    # fall back to any megahit on the share
    p=$(ls -d /ref/sahlab/software/miniforge3/envs/*/bin/megahit 2>/dev/null | head -1)
    [ -n "$p" ] && export PATH="$(dirname "$p"):$PATH"
}
command -v megahit >/dev/null || { echo "FATAL: megahit not on PATH" >&2; exit 1; }
command -v seqkit  >/dev/null || export PATH=/ref/sahlab/software/miniforge3/envs/bbtools/bin:$PATH

wd=$OUT/$s
rm -rf "$wd"

echo "=== $s ==="; echo "started: $(date)"

megahit -1 "$r1" -2 "$r2" \
    -o "$wd" --out-prefix "$s" \
    --k-min 21 --k-max 141 --k-step 12 \
    --min-contig-len 1000 \
    -t "${SLURM_CPUS_PER_TASK:-16}"

# Flatten to one predictable path and prefix contigs with the sample id, so the
# pooled FASTA can attribute every contig back to its sample.
#
# Also STRIP EVERYTHING AFTER THE FIRST SPACE. MEGAHIT writes headers like
#   >k141_30 flag=1 multi=7.0000 len=1244
# and the embedded spaces break phage-analysis's cluster_phages rule: it runs
# `seqkit fx2tab --name` (which keeps the whole description as field 1) and then
# `sort -k1,1` / `join -1 1 -2 1`, both of which split on whitespace. The join
# keys never match, the rule dies after "Extracting sequence lengths...", and
# Snakemake removes the output dir so there is no trace of why.
sed -e "s/^>/>${s}__/" -e "s/^\(>[^ ]*\).*/\1/" "$wd/${s}.contigs.fa" > "$final"

echo "contigs >=1000bp: $(grep -c '^>' "$final")"
echo "finished: $(date)"

# keep the flattened fasta, drop the heavy intermediate graph files
rm -rf "$wd/intermediate_contigs"
