#!/bin/bash
#SBATCH --job-name=votu_abund
#SBATCH --partition=general
#SBATCH --array=1-100%25
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --output=${PROJECT_DIR}/data/bulk/logs/abund_%A_%a.log
#SBATCH --error=${PROJECT_DIR}/data/bulk/logs/abund_%A_%a.err

# Per-sample vOTU abundance for the jackalope BULK cohort.
#
# WHY THIS EXISTS: phage-analysis reached 98% but its CoverM abundance step is
# among the unfinished rules, and clusters.tsv only lists the 212 cluster
# REPRESENTATIVES (49 of 100 samples) -- not every contributing contig. So there
# is currently no way to say how abundant each vOTU is in each sample, which is
# exactly what the HPV+ vs HPV- comparison needs.
#
# This maps each sample's reads back to the 210 vOTU representatives and reports
# coverage, so the comparison does not depend on the pipeline finishing.
#
# Output per sample: coverage depth + breadth + mapped read count per vOTU.
# Breadth matters as much as depth here: a vOTU with high depth over 5% of its
# length is a conserved-region artefact, not a present genome.

set -euo pipefail

VOTU=${PROJECT_DIR}/data/bulk/phage_analysis/out/02_clustering/vOTU_repSeqs.fasta
READS=${PROJECT_DIR}/data/bulk/READS
OUT=${PROJECT_DIR}/data/bulk/votu_abundance

mkdir -p "$OUT" ${PROJECT_DIR}/data/bulk/logs

i=${SLURM_ARRAY_TASK_ID}
s=$(printf "VIR_%03d" "$i")

r1=$READS/${s}_R1.fastq.gz
r2=$READS/${s}_R2.fastq.gz
[ -s "$r1" ] && [ -s "$r2" ] || { echo "FATAL: missing reads for $s" >&2; exit 1; }

cov=$OUT/${s}.cov_stats
[ -s "$cov" ] && { echo "SKIP: $s already done"; exit 0; }

# bbmap.sh is NOT in the hecatomb env itself -- Hecatomb gets it from a per-rule
# conda env that Snakemake builds. Use the standalone bbtools env on the share.
export PATH=${CONDA_BASE}/envs/bbtools/bin:$PATH
command -v bbmap.sh >/dev/null || { echo "FATAL: bbmap.sh not on PATH" >&2; exit 1; }

echo "=== $s ==="
echo "started: $(date)"

# minid=90 matches the Hecatomb QC mapping threshold used elsewhere in this
# project, so coverage numbers are comparable across steps.
bbmap.sh \
    ref="$VOTU" \
    in="$r1" in2="$r2" \
    nodisk \
    covstats="$cov" \
    rpkm="$OUT/${s}.rpkm" \
    scafstats="$OUT/${s}.scafstats" \
    statsfile="$OUT/${s}.statsfile" \
    ambiguous=random \
    minid=90 \
    physcov=t \
    threads="${SLURM_CPUS_PER_TASK:-8}" \
    -Xmx20g

echo "finished: $(date)"
echo "vOTUs with any coverage: $(awk -F'\t' 'NR>1 && $5>0' "$cov" | wc -l)"
