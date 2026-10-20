#!/usr/bin/env bash
# =============================================================================
# Step 2 -- run 01_qc.sh across a whole cohort and build the attrition table.
#
# Looks for  READS_DIR/<sample>_R1.fastq.gz  +  _R2.fastq.gz
# and produces, in OUTDIR:
#   <sample>_R{1,2}.clean.fastq.gz   per-sample cleaned reads
#   <sample>.readcounts.tsv          per-sample step counts
#   qc_attrition.tsv                 one row per sample  <- the per-sample table
#   qc_steps_total.tsv               cohort totals       <- the step table
#
# These are the tutorial equivalents of results/qc/qc_attrition.tsv and
# results/qc/qc_steps_vlp.tsv in this repo.
#
# Usage:
#   ./02_qc_cohort.sh -i READS_DIR -o OUTDIR [-r host.fa] [-p phix.fa]
#                     [-t THREADS] [-j PARALLEL_JOBS]
#
# Serial by default. -j N runs N samples at once (each using -t threads, so
# keep j*t at or below your core count).
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QC="$HERE/01_qc.sh"

READS=""; OUTDIR=""; HOST=""; PHIX=""; THREADS=4; JOBS=1
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) READS="$2"; shift 2 ;;
        -o) OUTDIR="$2"; shift 2 ;;
        -r) HOST="$2"; shift 2 ;;
        -p) PHIX="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -j) JOBS="$2"; shift 2 ;;
        -h|--help) sed -n '2,22p' "$0" | sed 's/^# \?//'; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$READS" && -n "$OUTDIR" ]] || { echo "ERROR: -i and -o required." >&2; exit 1; }
[[ -d "$READS" ]] || { echo "ERROR: no such directory: $READS" >&2; exit 1; }
[[ -x "$QC" ]] || { echo "ERROR: cannot execute $QC" >&2; exit 1; }

mkdir -p "$OUTDIR"

# ---- discover samples ------------------------------------------------------
# bash 3.2 (the macOS default) has no `mapfile`, so read into the array the
# portable way.
SAMPLES=()
while IFS= read -r line; do
    [[ -n "$line" ]] && SAMPLES+=("$line")
done < <(
    find "$READS" -maxdepth 1 \( -name '*_R1.fastq.gz' -o -name '*_R1.fq.gz' \) \
    | sed -E 's|.*/||; s/_R1\.(fastq|fq)\.gz$//' | sort -u
)
(( ${#SAMPLES[@]} )) || { echo "ERROR: no *_R1.fastq.gz found in $READS" >&2; exit 1; }
echo "Found ${#SAMPLES[@]} samples; ${JOBS} at a time, ${THREADS} threads each." >&2

args=()
[[ -n "$HOST" ]] && args+=(-r "$HOST")
[[ -n "$PHIX" ]] && args+=(-p "$PHIX")

run_one() {
    local s="$1"
    local r1 r2
    r1=$(ls "$READS/${s}_R1."{fastq,fq}".gz" 2>/dev/null | head -1)
    r2=${r1/_R1./_R2.}
    if [[ ! -f "$r2" ]]; then
        echo "  SKIP $s -- no R2 mate" >&2; return 0
    fi
    if [[ -f "$OUTDIR/${s}.readcounts.tsv" ]]; then
        echo "  skip $s (already done)" >&2; return 0
    fi
    "$QC" -1 "$r1" -2 "$r2" -s "$s" -o "$OUTDIR" -t "$THREADS" "${args[@]}" \
        > "$OUTDIR/${s}.stdout.log" 2> "$OUTDIR/${s}.stderr.log" \
        || { echo "  FAILED $s -- see $OUTDIR/${s}.stderr.log" >&2; return 1; }
    echo "  done $s" >&2
}
failed=0
# Simple job pool: works on bash 3.2, no GNU parallel or long xargs lines.
running=0
for s in "${SAMPLES[@]}"; do
    if (( JOBS > 1 )); then
        run_one "$s" || failed=1 &
        running=$(( running + 1 ))
        if (( running >= JOBS )); then
            wait -n 2>/dev/null || wait
            running=$(( running - 1 ))
        fi
    else
        run_one "$s" || failed=1
    fi
done
wait

# ---- per-sample attrition table -------------------------------------------
ATTR="$OUTDIR/qc_attrition.tsv"
printf 'sample\traw\tafter_fastp\tafter_host\tfinal\tpct_retained\n' > "$ATTR"
for s in "${SAMPLES[@]}"; do
    f="$OUTDIR/${s}.readcounts.tsv"
    [[ -f "$f" ]] || continue
    awk -F'\t' -v s="$s" '
        $1=="1"            { raw=$4; af=$5 }
        $1=="2"&&$5!="-"   { ah=$5 }
        $1=="total"        { fin=$5 }
        END {
            if (ah=="") ah=af
            printf "%s\t%d\t%d\t%d\t%d\t%.2f\n", s, raw, af, ah, fin,
                   (raw ? 100*fin/raw : 0)
        }' "$f" >> "$ATTR"
done

# ---- cohort totals ---------------------------------------------------------
TOT="$OUTDIR/qc_steps_total.tsv"
awk -F'\t' '
    NR>1 { raw+=$2; af+=$3; ah+=$4; fin+=$5; n++ }
    END {
        if (!n) { print "no samples" > "/dev/stderr"; exit 1 }
        printf "step\toperation\treads_in\treads_out\treads_removed\tpct_of_raw\n"
        printf "1\tfastp trim+filter+dedup\t%d\t%d\t%d\t%.2f%%\n", raw, af, raw-af, 100*(raw-af)/raw
        printf "2\thost removal\t%d\t%d\t%d\t%.2f%%\n",            af, ah, af-ah,  100*(af-ah)/raw
        printf "total\tretained\t%d\t%d\t%d\t%.2f%%\n",            raw, fin, raw-fin, 100*fin/raw
    }' "$ATTR" > "$TOT"

echo >&2
echo "=== cohort totals ($OUTDIR/qc_steps_total.tsv) ===" >&2
column -t -s$'\t' "$TOT" >&2
echo >&2
echo "Per-sample: $ATTR (${#SAMPLES[@]} samples)" >&2
(( failed == 0 )) || { echo "WARNING: at least one sample failed." >&2; exit 1; }
