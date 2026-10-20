#!/usr/bin/env bash
# =============================================================================
# Step 1 -- QC a paired-end virome, one sample at a time.
#
# Runs the operations in the QC tool list in docs/TOOLS_WALKTHROUGH.md:
#   1. fastp     adapters, quality, length, poly-X, N, duplicates
#   1b. bbduk    PhiX + vector            <- we skipped this; you should not
#   2. minimap2  host removal
#   3. repair.sh re-pair orphaned reads   <- most-forgotten step
#   4. mmseqs    dereplication (OFF by default; see --derep)
#
# No cluster, no Snakemake, no Hecatomb. Just the tools.
#
# Usage:
#   ./01_qc.sh -1 R1.fastq.gz -2 R2.fastq.gz -s SAMPLE -o outdir [options]
#
#   -r FASTA    host genome to remove (e.g. GRCh38). Omit to skip host removal.
#   -p FASTA    PhiX/vector reference. Omit to skip (simulated data only!).
#   -t INT      threads (default 8)
#   --derep     also run MMseqs2 dereplication (only if a big DB search follows)
#
# Writes outdir/SAMPLE_R{1,2}.clean.fastq.gz and outdir/SAMPLE.readcounts.tsv
# =============================================================================
set -euo pipefail

THREADS=8
DEREP=0
HOST=""
PHIX=""
R1=""; R2=""; SAMPLE=""; OUTDIR=""

usage() { sed -n '2,25p' "$0" | sed 's/^# \?//'; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        -1) R1="$2"; shift 2 ;;
        -2) R2="$2"; shift 2 ;;
        -s) SAMPLE="$2"; shift 2 ;;
        -o) OUTDIR="$2"; shift 2 ;;
        -r) HOST="$2"; shift 2 ;;
        -p) PHIX="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        --derep) DEREP=1; shift ;;
        -h|--help) usage 0 ;;
        *) echo "unknown argument: $1" >&2; usage 1 ;;
    esac
done

[[ -n "$R1" && -n "$R2" && -n "$SAMPLE" && -n "$OUTDIR" ]] || {
    echo "ERROR: -1, -2, -s and -o are all required." >&2; usage 1; }
[[ -f "$R1" ]] || { echo "ERROR: no such file: $R1" >&2; exit 1; }
[[ -f "$R2" ]] || { echo "ERROR: no such file: $R2" >&2; exit 1; }

# ---- preflight: fail now, not 20 minutes in -------------------------------
need=(fastp)
[[ -n "$HOST" ]] && need+=(minimap2 samtools)
[[ -n "$PHIX" ]] && need+=(bbduk.sh)
[[ -n "$HOST" ]] && need+=(repair.sh)
(( DEREP )) && need+=(mmseqs)
# tool -> conda package (they differ: repair.sh/bbduk.sh both ship in bbmap)
pkg_for() {
    case "$1" in
        repair.sh|bbduk.sh) echo bbmap ;;
        mmseqs)             echo mmseqs2 ;;
        *)                  echo "$1" ;;
    esac
}
# Collect missing tools as a newline-delimited string, not an array: bash 3.2
# with `set -u` treats an empty array as unbound, which crashed fastp-only runs.
missing=""
for t in "${need[@]}"; do
    command -v "$t" >/dev/null 2>&1 || missing="${missing}${t}"$'\n'
done

# repair.sh is optional -- samtools can re-pair instead (see step 3)
REPAIR_MODE="bbmap"
if [[ "$missing" == *"repair.sh"* ]] && command -v samtools >/dev/null 2>&1; then
    REPAIR_MODE="samtools"
    missing=$(printf '%s' "$missing" | grep -v '^repair\.sh$' || true)
fi

missing=$(printf '%s' "$missing" | grep -v '^$' || true)
if [[ -n "$missing" ]]; then
    echo "ERROR: not on PATH: $(printf '%s' "$missing" | tr '\n' ' ')" >&2
    pkgs=$(printf '%s\n' "$missing" | while read -r t; do
               [[ -n "$t" ]] && pkg_for "$t"; done | sort -u | tr '\n' ' ')
    echo "Install with: conda install -c bioconda -c conda-forge $pkgs" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/qc_${SAMPLE}_XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

COUNTS="$OUTDIR/${SAMPLE}.readcounts.tsv"
printf 'step\toperation\ttool\treads_in\treads_out\treads_removed\n' > "$COUNTS"

# count reads in a fastq (handles .gz); 4 lines per record
countfq() {
    local f="$1"
    if [[ "$f" == *.gz ]]; then gzip -cd "$f" | wc -l; else wc -l < "$f"; fi \
        | awk '{print int($1/4)}'
}
log() { printf '[%s] %s\n' "$SAMPLE" "$*" >&2; }

RAW=$(( $(countfq "$R1") * 2 ))
log "raw reads: $RAW"

# ---- 1. fastp -------------------------------------------------------------
log "step 1: fastp"
fastp -i "$R1" -I "$R2" \
      -o "$WORK/f_R1.fastq" -O "$WORK/f_R2.fastq" \
      -j "$OUTDIR/${SAMPLE}.fastp.json" -h "$OUTDIR/${SAMPLE}.fastp.html" \
      --qualified_quality_phred 15 \
      --length_required 90 \
      --detect_adapter_for_pe \
      --cut_tail --cut_tail_window_size 25 --cut_tail_mean_quality 15 \
      --dedup --dup_calc_accuracy 4 \
      --trim_poly_x \
      --thread "$THREADS" 2> "$OUTDIR/${SAMPLE}.fastp.log"

AFTER_FASTP=$(( $(countfq "$WORK/f_R1.fastq") * 2 ))
printf '1\ttrim+filter+dedup\tfastp\t%d\t%d\t%d\n' \
    "$RAW" "$AFTER_FASTP" "$(( RAW - AFTER_FASTP ))" >> "$COUNTS"

# fastp's JSON holds the per-reason breakdown -- pull it out if jq is present
if command -v jq >/dev/null 2>&1; then
    J="$OUTDIR/${SAMPLE}.fastp.json"
    for key in low_quality_reads too_many_N_reads too_short_reads; do
        v=$(jq -r ".filtering_result.${key} // 0" "$J")
        printf '1.x\tdiscarded: %s\tfastp\t-\t-\t%s\n' "$key" "$v" >> "$COUNTS"
    done
    v=$(jq -r '.duplication.rate // 0' "$J")
    printf '1.x\tduplication rate\tfastp\t-\t-\t%s\n' "$v" >> "$COUNTS"
fi

CUR1="$WORK/f_R1.fastq"; CUR2="$WORK/f_R2.fastq"

# ---- 1b. PhiX / vector ----------------------------------------------------
if [[ -n "$PHIX" ]]; then
    log "step 1b: bbduk PhiX/vector"
    bbduk.sh in="$CUR1" in2="$CUR2" ref="$PHIX" \
             out="$WORK/x_R1.fastq" out2="$WORK/x_R2.fastq" \
             stats="$OUTDIR/${SAMPLE}.phix.stats" \
             k=31 hammingdistance=1 ordered=t ow=t \
             threads="$THREADS" 2> "$OUTDIR/${SAMPLE}.bbduk.log"
    AFTER_PHIX=$(( $(countfq "$WORK/x_R1.fastq") * 2 ))
    printf '1b\tremove PhiX+vector\tbbduk\t%d\t%d\t%d\n' \
        "$AFTER_FASTP" "$AFTER_PHIX" "$(( AFTER_FASTP - AFTER_PHIX ))" >> "$COUNTS"
    CUR1="$WORK/x_R1.fastq"; CUR2="$WORK/x_R2.fastq"
else
    log "step 1b: SKIPPED (no -p reference given)"
    printf '1b\tremove PhiX+vector\tbbduk\t-\t-\tSKIPPED\n' >> "$COUNTS"
fi

PREV=$(( $(countfq "$CUR1") * 2 ))

# ---- 2 + 3. host removal, then re-pair ------------------------------------
if [[ -n "$HOST" ]]; then
    IDX="$WORK/host.idx"
    log "step 2: indexing host (once per reference -- cache this for real runs)"
    minimap2 -t "$THREADS" -d "$IDX" "$HOST" 2> "$OUTDIR/${SAMPLE}.minimap2.log"

    log "step 2: host removal"
    minimap2 -ax sr -t "$THREADS" --secondary=no "$IDX" "$CUR1" "$CUR2" \
             2>> "$OUTDIR/${SAMPLE}.minimap2.log" \
      | samtools view -f 4 -h - \
      | samtools fastq -NO -1 "$WORK/h_R1.fastq" -2 "$WORK/h_R2.fastq" - \
             2>> "$OUTDIR/${SAMPLE}.minimap2.log"

    AFTER_HOST=$(( $(countfq "$WORK/h_R1.fastq") + $(countfq "$WORK/h_R2.fastq") ))
    printf '2\tremove host\tminimap2+samtools\t%d\t%d\t%d\n' \
        "$PREV" "$AFTER_HOST" "$(( PREV - AFTER_HOST ))" >> "$COUNTS"

    # -f 4 filters mates independently, so R1/R2 are now out of sync.
    log "step 3: repair pairing (via $REPAIR_MODE)"
    if [[ "$REPAIR_MODE" == "bbmap" ]]; then
        repair.sh in="$WORK/h_R1.fastq" in2="$WORK/h_R2.fastq" \
                  out="$WORK/p_R1.fastq" out2="$WORK/p_R2.fastq" \
                  outs="$OUTDIR/${SAMPLE}.singletons.fastq" \
                  repair -Xmx8g 2> "$OUTDIR/${SAMPLE}.repair.log"
    else
        # samtools-only fallback (no bbtools needed).
        # Re-map and keep only pairs where BOTH mates are unmapped:
        #   -f 12 = read unmapped (4) + mate unmapped (8)
        #   -F 256 = drop secondary alignments
        # Sorting by name (-n) guarantees R1/R2 come out in sync.
        minimap2 -ax sr -t "$THREADS" --secondary=no "$IDX" "$CUR1" "$CUR2" \
                 2>> "$OUTDIR/${SAMPLE}.minimap2.log" \
          | samtools view -b -f 12 -F 256 - \
          | samtools sort -n -@ "$THREADS" -o "$WORK/both.bam" - \
                 2>> "$OUTDIR/${SAMPLE}.repair.log"
        samtools fastq -N -1 "$WORK/p_R1.fastq" -2 "$WORK/p_R2.fastq" \
                 -0 /dev/null -s /dev/null "$WORK/both.bam" \
                 2>> "$OUTDIR/${SAMPLE}.repair.log"
        : > "$OUTDIR/${SAMPLE}.singletons.fastq"
    fi
    AFTER_REPAIR=$(( $(countfq "$WORK/p_R1.fastq") * 2 ))
    repair_tool=$([[ "$REPAIR_MODE" == "bbmap" ]] && echo "repair.sh" || echo "samtools (fallback)")
    printf '3\tre-pair orphans\t%s\t%d\t%d\t%d\n' \
        "$repair_tool" "$AFTER_HOST" "$AFTER_REPAIR" "$(( AFTER_HOST - AFTER_REPAIR ))" >> "$COUNTS"
    CUR1="$WORK/p_R1.fastq"; CUR2="$WORK/p_R2.fastq"
else
    log "step 2-3: SKIPPED (no -r host reference given)"
    printf '2\tremove host\tminimap2+samtools\t-\t-\tSKIPPED\n' >> "$COUNTS"
fi

# ---- 4. dereplication (optional) ------------------------------------------
if (( DEREP )); then
    log "step 4: MMseqs2 dereplication"
    mmseqs easy-linclust "$CUR1" "$OUTDIR/${SAMPLE}_rep" "$WORK/mmtmp" \
        --kmer-per-seq-scale 0.3 -c 0.8 --cov-mode 1 \
        --min-seq-id 0.97 --alignment-mode 3 \
        --threads "$THREADS" > "$OUTDIR/${SAMPLE}.mmseqs.log" 2>&1
    REPS=$(grep -c '^>' "$OUTDIR/${SAMPLE}_rep_rep_seq.fasta" 2>/dev/null || echo 0)
    printf '4\tdereplicate 97%%\tMMseqs2\t-\t%s representatives\t0 (counts preserved)\n' \
        "$REPS" >> "$COUNTS"
fi

# ---- output ---------------------------------------------------------------
gzip -c "$CUR1" > "$OUTDIR/${SAMPLE}_R1.clean.fastq.gz"
gzip -c "$CUR2" > "$OUTDIR/${SAMPLE}_R2.clean.fastq.gz"

FINAL=$(( $(countfq "$OUTDIR/${SAMPLE}_R1.clean.fastq.gz") * 2 ))
printf 'total\tretained\t-\t%d\t%d\t%d\n' \
    "$RAW" "$FINAL" "$(( RAW - FINAL ))" >> "$COUNTS"

log "done: $RAW -> $FINAL reads ($(awk -v a="$FINAL" -v b="$RAW" 'BEGIN{printf "%.2f", b?100*a/b:0}')% retained)"
log "counts: $COUNTS"
column -t -s$'\t' "$COUNTS" >&2
