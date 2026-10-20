#!/usr/bin/env python3
"""Summarise a Hecatomb preprocessing run: tools, reads removed, totals.

Produces three tables from a qc_attrition*.tsv file:

  1. STEP TABLE     -- one row per QC step: the software that ran, its version,
                       what it does, and how many reads it removed cohort-wide.
  2. GROUP SUMMARY  -- the same attrition split by HPV status.
  3. TOTAL SUMMARY  -- cohort totals and the overall funnel.

Tool versions are read from the Hecatomb install's conda env specs, so they
describe what actually ran rather than what the docs claim.

Usage:
    python3 qc_summary.py qc_attrition_bulk.tsv -o results/qc/ [--label bulk]
"""
import argparse
import csv
import os
import sys
from pathlib import Path

# The --preprocess paired route. Each entry: (step key, tool, env yaml, what it does)
STEPS = [
    ("raw", "—", None,
     "Input read pairs as generated/sequenced"),
    ("after_fastp", "fastp", "fastp.yaml",
     "Adapter removal, quality trim (Q15, 25bp window), length filter (>=90bp), "
     "poly-X trim, deduplication"),
    ("after_host", "minimap2 + samtools", "minimap2.yaml",
     "Map to viral-masked human reference (GRCh38); keep unmapped reads only "
     "(samtools view -f 4)"),
    ("after_cluster", "MMseqs2 linclust", "mmseqs2.yaml",
     "Dereplicate at 97% identity / 80% coverage; counts carried onto "
     "representatives (lossless for reads)"),
]

# Fallback versions if the env specs cannot be read
FALLBACK_VERSIONS = {
    "fastp.yaml": "fastp 0.23.2",
    "minimap2.yaml": "minimap2 >=2.20, samtools >=1.13",
    "mmseqs2.yaml": "mmseqs2 12.113e3",
}


def tool_versions(env_dir):
    """Parse `dependencies:` out of the Hecatomb conda env yamls."""
    versions = {}
    if not env_dir or not Path(env_dir).is_dir():
        return versions
    for f in Path(env_dir).glob("*.yaml"):
        deps, in_deps = [], False
        try:
            for line in f.read_text().splitlines():
                s = line.strip()
                if s.startswith("dependencies:"):
                    in_deps = True
                    continue
                if in_deps:
                    if s.startswith("- "):
                        deps.append(s[2:].split("=")[0] + (
                            "=" + s[2:].split("=")[1] if "=" in s[2:] else ""))
                    elif s and not s.startswith("#"):
                        break
        except OSError:
            continue
        if deps:
            versions[f.name] = ", ".join(deps)
    return versions


def load(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def icol(rows, key, group=None):
    """Sum an integer column, optionally for one group."""
    total = 0
    for r in rows:
        if group and r["group"] != group:
            continue
        v = r.get(key, "NA")
        if v not in ("NA", ""):
            total += int(v)
    return total


def fmt(n):
    return f"{n:,}"


def pct(a, b):
    return f"{100.0 * a / b:.2f}%" if b else "—"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("attrition", help="qc_attrition*.tsv from collect_qc_attrition.py")
    ap.add_argument("-o", "--outdir", default=".")
    ap.add_argument("--label", default="cohort", help="name for output files")
    ap.add_argument("--envs",
                    default="${HECATOMB_ENV}/"
                            "snakemake/workflow/envs",
                    help="Hecatomb conda env directory (for tool versions)")
    args = ap.parse_args()

    rows = load(args.attrition)
    if not rows:
        sys.exit(f"error: no rows in {args.attrition}")

    groups = sorted({r["group"] for r in rows})
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    versions = tool_versions(args.envs)

    # ---------------------------------------------------------------- 1. steps
    step_path = outdir / f"qc_steps_{args.label}.tsv"
    step_cols = ["step", "software", "version", "purpose", "reads_in", "reads_out",
                 "reads_removed", "pct_removed", "pct_of_raw_remaining"]

    raw_total = icol(rows, "raw")
    step_rows = []
    prev_key = None
    for i, (key, tool, env, purpose) in enumerate(STEPS):
        cur = icol(rows, key)
        if i == 0:
            step_rows.append({
                "step": "0. input", "software": tool, "version": "—",
                "purpose": purpose, "reads_in": "—", "reads_out": fmt(cur),
                "reads_removed": "—", "pct_removed": "—",
                "pct_of_raw_remaining": "100.00%",
            })
        else:
            prev = icol(rows, prev_key)
            removed = prev - cur
            step_rows.append({
                "step": f"{i}. {key.replace('after_', '')}",
                "software": tool,
                "version": versions.get(env, FALLBACK_VERSIONS.get(env, "—")),
                "purpose": purpose,
                "reads_in": fmt(prev), "reads_out": fmt(cur),
                "reads_removed": fmt(removed),
                "pct_removed": pct(removed, prev),
                "pct_of_raw_remaining": pct(cur, raw_total),
            })
        prev_key = key

    with open(step_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=step_cols, delimiter="\t")
        w.writeheader()
        w.writerows(step_rows)

    # ---------------------------------------------------------------- 2. totals
    total_path = outdir / f"qc_totals_{args.label}.tsv"
    total_cols = ["scope", "n_samples", "raw", "after_fastp", "after_host",
                  "after_cluster", "n_representatives",
                  "removed_fastp", "removed_host", "removed_cluster",
                  "removed_total", "pct_retained"]

    def scope_row(name, group=None):
        sub = [r for r in rows if (group is None or r["group"] == group)]
        raw = icol(rows, "raw", group)
        fp = icol(rows, "after_fastp", group)
        ho = icol(rows, "after_host", group)
        cl = icol(rows, "after_cluster", group)
        reps = icol(rows, "n_representatives", group)
        return {
            "scope": name, "n_samples": len(sub),
            "raw": fmt(raw), "after_fastp": fmt(fp), "after_host": fmt(ho),
            "after_cluster": fmt(cl), "n_representatives": fmt(reps),
            "removed_fastp": fmt(raw - fp),
            "removed_host": fmt(fp - ho),
            "removed_cluster": fmt(ho - cl),
            "removed_total": fmt(raw - cl),
            "pct_retained": pct(cl, raw),
        }

    total_rows = [scope_row(g, g) for g in groups]
    total_rows.append(scope_row("TOTAL"))

    with open(total_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=total_cols, delimiter="\t")
        w.writeheader()
        w.writerows(total_rows)

    # ---------------------------------------------------------------- print
    print(f"\n{'=' * 78}\n QC STEP SUMMARY -- {args.label}\n{'=' * 78}")
    print(f"{'step':16} {'software':22} {'reads in':>13} {'removed':>13} {'% removed':>10}")
    print("-" * 78)
    for r in step_rows:
        print(f"{r['step']:16} {r['software']:22} {r['reads_in']:>13} "
              f"{r['reads_removed']:>13} {r['pct_removed']:>10}")

    print(f"\n{'=' * 78}\n TOTALS\n{'=' * 78}")
    print(f"{'scope':10} {'n':>4} {'raw':>14} {'after host':>14} {'final':>14} {'retained':>10}")
    print("-" * 78)
    for r in total_rows:
        print(f"{r['scope']:10} {r['n_samples']:>4} {r['raw']:>14} "
              f"{r['after_host']:>14} {r['after_cluster']:>14} {r['pct_retained']:>10}")

    print(f"\nwrote {step_path}")
    print(f"wrote {total_path}")


if __name__ == "__main__":
    main()
