#!/usr/bin/env python3
"""Build a per-step read attrition table from a Hecatomb preprocessing run.

Hecatomb v1.1.0 ships rules that count reads at each step
(04_summaries_optional.smk), but the `include:` is commented out in the HTCF
install, so they never run. The numbers still exist in the per-step outputs --
this collects them into one table.

Steps tracked (--preprocess paired):
    raw            read pairs in the input FASTQ      (fastp JSON)
    after_fastp    quality/adapter trim + dedup       (fastp JSON)
    after_host     human host removal, viral-masked   (assembly/*_R1.all.fastq.gz)
    after_cluster  dereplication at 97% id            (results/sampleSeqCounts.tsv)
    n_representatives  unique sequences kept          (results/seqtable.fasta)

Counts are REA PAIRS throughout: fastp reports total_reads (both mates), so it
is halved; the host stage counts R1 only; clustering is already per-pair.

Usage:
    python3 collect_qc_attrition.py QC_DIR METADATA_CSV -o qc_attrition.tsv
"""
import argparse
import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path


def read_metadata(path):
    """Sample -> group, from the blinded clinical sheet."""
    groups = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            antler = row.get("Antler_Present", "").strip()
            groups[row["Sample_ID"]] = "HPV+" if antler == "Yes" else "HPV-"
    return groups


def fastp_counts(stats_dir, sample):
    """(raw_pairs, after_fastp_pairs) from fastp's JSON report."""
    f = stats_dir / "p01" / f"{sample}.s1.stats.json"
    if not f.exists():
        return None, None
    with open(f) as fh:
        d = json.load(fh)
    s = d.get("summary", {})
    before = s.get("before_filtering", {}).get("total_reads")
    after = s.get("after_filtering", {}).get("total_reads")
    # fastp counts both mates; we report pairs
    return (before // 2 if before else None, after // 2 if after else None)


def count_fastq_gz(path):
    """Records in a gzipped FASTQ. Uses zcat|wc -l -- much faster than Python."""
    if not path.exists():
        return None
    try:
        out = subprocess.check_output(
            f"zcat {path} | wc -l", shell=True, stderr=subprocess.DEVNULL
        )
        return int(out.strip()) // 4
    except Exception:
        return None


def host_counts(assembly_dir, sample):
    """Read pairs surviving host removal.

    *_R1.all.fastq.gz is the post-host R1 set (paired + singletons), which is
    what carries forward into clustering.
    """
    return count_fastq_gz(assembly_dir / f"{sample}_R1.all.fastq.gz")


def cluster_counts(results_dir):
    """sample -> reads assigned to representative sequences."""
    f = results_dir / "sampleSeqCounts.tsv"
    counts = {}
    if not f.exists():
        return counts
    with open(f) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    counts[parts[0]] = int(float(parts[1]))
                except ValueError:
                    pass
    return counts


def representative_counts(results_dir):
    """sample -> number of representative sequences, from seqtable.fasta headers.

    Headers look like  >VIR_088:1:1.503e-04:0
    """
    f = results_dir / "seqtable.fasta"
    reps = {}
    if not f.exists():
        return reps
    try:
        out = subprocess.check_output(
            f"grep '^>' {f} | cut -d: -f1 | tr -d '>' | sort | uniq -c",
            shell=True, stderr=subprocess.DEVNULL,
        ).decode()
        for line in out.splitlines():
            parts = line.split()
            if len(parts) == 2:
                reps[parts[1]] = int(parts[0])
    except Exception:
        pass
    return reps


def pct(num, den):
    return f"{100.0 * num / den:.2f}" if (num not in (None, "NA") and den) else "NA"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("qc_dir", help="Hecatomb output directory")
    ap.add_argument("metadata", help="jackalope_clinical_metadata.csv")
    ap.add_argument("-o", "--out", default="qc_attrition.tsv")
    args = ap.parse_args()

    qc = Path(args.qc_dir)
    if not qc.exists():
        sys.exit(f"error: no such QC directory: {qc}")

    stats = qc / "processing" / "stats"
    assembly = qc / "processing" / "assembly"
    results = qc / "results"

    groups = read_metadata(args.metadata)
    clustered = cluster_counts(results)
    reps = representative_counts(results)

    cols = ["sample", "group", "raw", "after_fastp", "after_host", "after_cluster",
            "n_representatives", "pct_fastp", "pct_host", "pct_cluster", "pct_retained"]

    rows = []
    for i, sample in enumerate(sorted(groups), 1):
        raw, after_fastp = fastp_counts(stats, sample)
        after_host = host_counts(assembly, sample)
        after_cluster = clustered.get(sample)

        if raw is None and after_cluster is None:
            continue

        rows.append({
            "sample": sample,
            "group": groups[sample],
            "raw": raw if raw is not None else "NA",
            "after_fastp": after_fastp if after_fastp is not None else "NA",
            "after_host": after_host if after_host is not None else "NA",
            "after_cluster": after_cluster if after_cluster is not None else "NA",
            "n_representatives": reps.get(sample, "NA"),
            "pct_fastp": pct(after_fastp, raw),
            "pct_host": pct(after_host, after_fastp),
            "pct_cluster": pct(after_cluster, after_host),
            "pct_retained": pct(after_cluster, raw),
        })
        if i % 20 == 0:
            print(f"  ...{i}/{len(groups)} samples", file=sys.stderr)

    if not rows:
        sys.exit("error: no per-sample stats found -- has the QC run finished?")

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"\nwrote {args.out}  ({len(rows)} samples)\n")

    # Group means. The host-removal rate is a confound check: both groups were
    # simulated with the same 20% human contamination, so a systematic
    # difference here would mean bias entering before any analysis.
    def mean(sub, key):
        vals = [float(r[key]) for r in sub if r[key] != "NA"]
        return sum(vals) / len(vals) if vals else float("nan")

    def imean(sub, key):
        vals = [int(r[key]) for r in sub if r[key] != "NA"]
        return sum(vals) / len(vals) if vals else float("nan")

    print("Mean % reads retained at each step, by group")
    print(f"  {'group':6} {'n':>4} {'fastp':>8} {'host':>8} {'cluster':>9} {'overall':>9}")
    for g in ("HPV-", "HPV+"):
        sub = [r for r in rows if r["group"] == g]
        if not sub:
            continue
        print(f"  {g:6} {len(sub):>4} {mean(sub,'pct_fastp'):>8.2f} {mean(sub,'pct_host'):>8.2f} "
              f"{mean(sub,'pct_cluster'):>9.2f} {mean(sub,'pct_retained'):>9.2f}")

    print("\nMean read pairs, by group")
    print(f"  {'group':6} {'n':>4} {'raw':>12} {'after_fastp':>12} {'after_host':>12} "
          f"{'after_clust':>12} {'n_reps':>10}")
    for g in ("HPV-", "HPV+"):
        sub = [r for r in rows if r["group"] == g]
        if not sub:
            continue
        print(f"  {g:6} {len(sub):>4} {imean(sub,'raw'):>12,.0f} {imean(sub,'after_fastp'):>12,.0f} "
              f"{imean(sub,'after_host'):>12,.0f} {imean(sub,'after_cluster'):>12,.0f} "
              f"{imean(sub,'n_representatives'):>10,.0f}")

    print("\nNote: depth differs by group BY DESIGN (HPV+ 5-20x, HPV- flat 10x).")
    print("Compare the percentage columns, not the raw counts.")


if __name__ == "__main__":
    main()
