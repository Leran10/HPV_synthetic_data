#!/bin/bash
#SBATCH --job-name=bulk_cohort
#SBATCH --partition=general
#SBATCH --array=1-100%25
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/sahlab/leran/jackalope_bulk/logs/gen_%A_%a.log
#SBATCH --error=/scratch/sahlab/leran/jackalope_bulk/logs/gen_%A_%a.err

# Regenerate the jackalope cohort as a BULK metagenome (no VLP enrichment).
#
# Why: the original cohort used VLP tangential-flow enrichment, whose
# nuclease_efficiency=0.98 stripped host DNA from 20% down to 0.67%. That made
# the QC host-removal step look like a no-op, which is the wrong lesson for a
# tutorial -- real bulk cervicovaginal metagenomes run 40-95% human
# (e.g. PMC9040584 measured 61%).
#
# --no-vlp --contamination-level heavy delivers, per pilot:
#     bacterial 48.8% | host 41.1% | rRNA 5.5% | fungal 3.6% | viral 0.98%
# which matches the published ~0.83% viral fraction for vaginal viromes and
# keeps the bacterial fraction that prophage-pipeline needs for binning.
#
# Sample identity, group assignment, seeds and per-sample coverage are taken
# from jackalope_clinical_metadata.csv, exactly as the original cohort was
# built (coverage scales 5x->20x with Antler_Size_mm for HPV+ samples;
# HPV- samples are flat 10x). So VIR_### means the same sample as before.
#
# Array is throttled to 25 concurrent (%25) to stay a good neighbour on the
# shared partition. Resources are deliberately generous: the pilot needed
# >40 min at 5x, and bulk mode simulates ~100x more non-viral material than
# the VLP run did, so the top-coverage samples are the long pole.

set -euo pipefail

VIROFORGE=/scratch/sahlab/leran/viroforge
VFENV=/scratch/sahlab/leran/miniforge3/envs/vf_benchmark

# Put the env's bin on PATH, not just call viroforge by absolute path:
# viroforge shells out to InSilicoSeq (`iss`), which it looks up via PATH.
# Calling the absolute path alone leaves iss unfindable and every sample dies
# at the read-simulation step with "InSilicoSeq (iss) not found in PATH".
export PATH="$VFENV/bin:$PATH"
VF="$VFENV/bin/viroforge"

command -v iss >/dev/null || { echo "FATAL: iss not on PATH" >&2; exit 1; }
METADATA=/scratch/sahlab/leran/jackalope_cohort/jackalope_clinical_metadata.csv
OUTBASE=/scratch/sahlab/leran/jackalope_bulk

# Same mapping the original cohort used (examples/batch_configs/jackalope_cohort.py)
MIN_SIZE=8; MAX_SIZE=63
MIN_COV=5; MAX_COV=20
COLL_NEG=16        # HPV-negative  (Vaginal Virome, healthy)
COLL_POS=24        # HPV-positive  (HPV-Infected Vaginal)

mkdir -p "$OUTBASE/logs"

i=${SLURM_ARRAY_TASK_ID}
sample=$(printf "VIR_%03d" "$i")

# Pull this sample's row out of the metadata
row=$(awk -F, -v s="$sample" 'NR>1 && $1==s {print; exit}' "$METADATA")
if [ -z "$row" ]; then
    echo "FATAL: $sample not found in $METADATA" >&2
    exit 1
fi

antler=$(echo "$row" | cut -d, -f6)
size=$(echo "$row" | cut -d, -f7)
seed=$((1000 + i))          # same seed rule as the original run

if [ "$antler" = "Yes" ]; then
    collection=$COLL_POS
    # coverage scales linearly with antler size
    coverage=$(awk -v s="$size" -v mn=$MIN_SIZE -v mx=$MAX_SIZE \
                   -v cl=$MIN_COV -v ch=$MAX_COV \
                   'BEGIN{printf "%.1f", cl + ((s-mn)/(mx-mn))*(ch-cl)}')
    group="HPV+"
else
    collection=$COLL_NEG
    coverage=10
    group="HPV-"
fi

outdir="$OUTBASE/$sample"

echo "=================================================================="
echo " $sample  ($group)"
echo "   collection : $collection"
echo "   coverage   : ${coverage}x"
echo "   seed       : $seed"
echo "   antler     : $antler ${size}mm"
echo "   output     : $outdir"
echo "   started    : $(date)"
echo "=================================================================="

# Idempotent: skip if this sample already produced both mates
if ls "$outdir"/fastq/*_R1.fastq* >/dev/null 2>&1 && \
   ls "$outdir"/fastq/*_R2.fastq* >/dev/null 2>&1; then
    echo "SKIP: $sample already has paired FASTQ output"
    exit 0
fi

cd "$VIROFORGE"

"$VF" generate \
    --collection-id "$collection" \
    --platform novaseq \
    --coverage "$coverage" \
    --seed "$seed" \
    --no-vlp \
    --contamination-level heavy \
    --output "$outdir"

echo
echo "=== delivered read composition (sampled) ==="
r1=$(ls "$outdir"/fastq/*_R1.fastq 2>/dev/null | head -1)
if [ -n "$r1" ]; then
    head -400000 "$r1" | awk 'NR%4==1' \
        | grep -o 'source=[a-z_]*' | sort | uniq -c | sort -rn
    echo "--- R1 size: $(du -h "$r1" | cut -f1)"
else
    echo "WARNING: no R1 FASTQ produced for $sample" >&2
    exit 1
fi

echo
echo "finished: $(date)"
