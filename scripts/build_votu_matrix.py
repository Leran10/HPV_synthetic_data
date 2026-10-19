#!/usr/bin/env python3
"""Assemble the analysis-ready vOTU table for the jackalope VLP cohort.

Joins four things that currently live apart:
  - per-sample coverage      votu_abundance/VIR_*.cov_stats   (bbmap)
  - taxonomy                 phabox_output/.../phagcn_prediction.tsv
  - predicted host           phabox_output/.../cherry_prediction.tsv
  - lifestyle                03_genomic_info/lifestyle_consensus.tsv
  - completeness/quality     03_checkv_final/quality_summary.tsv

Writes:
  votu_abundance_matrix.tsv   vOTU x sample mean coverage  (2921 x 100)
  votu_presence_matrix.tsv    vOTU x sample 0/1 presence
  votu_annotations.tsv        one row per vOTU, all annotation tracks merged
  votu_richness.tsv           per-sample richness + group, for the comparison

PRESENCE needs both a depth and a BREADTH threshold. High coverage over a few
percent of a contig is a conserved-region artefact, not a genome that is there,
so breadth is the stricter and more important of the two.

NOTE ON THE COMPARISON: HPV+ samples were simulated at 5-20x and HPV- at a flat
10x, so HPV+ carries ~1.5x the read depth by construction. Richness rises with
depth, and richness is the expected signal -- so raw richness is confounded.
This script reports depth alongside richness; rarefy or model depth before
drawing conclusions.
"""
import argparse
import csv
import sys
from pathlib import Path

MIN_DEPTH = 1.0      # mean fold coverage
MIN_BREADTH = 50.0   # percent of the vOTU covered


def read_groups(path):
    """Sample -> HPV status, from the blinded clinical sheet."""
    g = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            g[row["Sample_ID"]] = ("HPV+" if row.get("Antler_Present", "").strip() == "Yes"
                                   else "HPV-")
    return g


def read_cov(path):
    """vOTU -> (avg_fold, covered_percent) from one bbmap cov_stats file."""
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            vid = row.get("#ID") or row.get("ID")
            if not vid:
                continue
            try:
                out[vid] = (float(row["Avg_fold"]), float(row["Covered_percent"]))
            except (KeyError, ValueError):
                continue
    return out


def read_keyed(path, key_field, want, prefix):
    """Generic annotation reader: key -> {prefix_col: value} for wanted columns."""
    out = {}
    p = Path(path)
    if not p.exists():
        return out
    delim = "," if p.suffix == ".csv" else "\t"
    with open(p) as f:
        for row in csv.DictReader(f, delimiter=delim):
            k = row.get(key_field)
            if not k:
                continue
            out[k] = {f"{prefix}_{c}": row.get(c, "") for c in want if c in row}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abund-dir", required=True)
    ap.add_argument("--pipeline-out", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("-o", "--outdir", default=".")
    ap.add_argument("--min-depth", type=float, default=MIN_DEPTH)
    ap.add_argument("--min-breadth", type=float, default=MIN_BREADTH)
    args = ap.parse_args()

    O = Path(args.pipeline_out)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    groups = read_groups(args.metadata)

    # ---- per-sample coverage
    cov = {}
    for s in sorted(groups):
        f = Path(args.abund_dir) / f"{s}.cov_stats"
        if f.exists():
            cov[s] = read_cov(f)
    if not cov:
        sys.exit(f"error: no .cov_stats in {args.abund_dir}")
    samples = sorted(cov)
    votus = sorted({v for c in cov.values() for v in c})

    # ---- annotations
    pb = O / "03_genomic_info/phabox_output/final_prediction"
    tax = read_keyed(pb / "phagcn_prediction.tsv", "Accession",
                     ["Lineage", "Genus", "GenusCluster", "PhaGCNScore"], "tax")
    host = read_keyed(pb / "cherry_prediction.tsv", "Accession",
                      ["Host", "CHERRYScore", "Method", "Host_NCBI_lineage"], "host")
    life = read_keyed(O / "03_genomic_info/lifestyle_consensus.tsv", "contig_id",
                      ["lifestyle", "confidence", "source"], "life")
    qual = read_keyed(O / "03_checkv_final/quality_summary.tsv", "contig_id",
                      ["contig_length", "gene_count", "viral_genes", "host_genes",
                       "checkv_quality", "completeness"], "checkv")

    # ---- abundance + presence matrices
    present = {s: set() for s in samples}
    with open(outdir / "votu_abundance_matrix.tsv", "w", newline="") as fa, \
         open(outdir / "votu_presence_matrix.tsv", "w", newline="") as fp:
        wa = csv.writer(fa, delimiter="\t")
        wp = csv.writer(fp, delimiter="\t")
        wa.writerow(["vOTU"] + samples)
        wp.writerow(["vOTU"] + samples)
        for v in votus:
            arow, prow = [v], [v]
            for s in samples:
                d, b = cov[s].get(v, (0.0, 0.0))
                arow.append(f"{d:.4f}")
                hit = 1 if (d >= args.min_depth and b >= args.min_breadth) else 0
                prow.append(hit)
                if hit:
                    present[s].add(v)
            wa.writerow(arow)
            wp.writerow(prow)

    # ---- annotation table
    acols = (["vOTU", "n_samples_present"]
             + sorted({k for d in tax.values() for k in d})
             + sorted({k for d in host.values() for k in d})
             + sorted({k for d in life.values() for k in d})
             + sorted({k for d in qual.values() for k in d}))
    with open(outdir / "votu_annotations.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=acols, delimiter="\t", restval="")
        w.writeheader()
        for v in votus:
            row = {"vOTU": v,
                   "n_samples_present": sum(1 for s in samples if v in present[s])}
            for src in (tax, host, life, qual):
                row.update(src.get(v, {}))
            w.writerow(row)

    # ---- per-sample richness, with depth so the confound stays visible
    with open(outdir / "votu_richness.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "group", "votu_richness", "votus_any_coverage",
                    "total_coverage"])
        for s in samples:
            any_cov = sum(1 for d, _ in cov[s].values() if d > 0)
            tot = sum(d for d, _ in cov[s].values())
            w.writerow([s, groups[s], len(present[s]), any_cov, f"{tot:.2f}"])

    # ---- report
    pos = [s for s in samples if groups[s] == "HPV+"]
    neg = [s for s in samples if groups[s] == "HPV-"]

    def mean(xs):
        return sum(xs) / len(xs) if xs else float("nan")

    print(f"\n{'=' * 66}")
    print(" vOTU TABLE BUILT")
    print("=" * 66)
    print(f"  vOTUs          : {len(votus)}")
    print(f"  samples        : {len(samples)}")
    print(f"  presence rule  : depth >= {args.min_depth}x AND breadth >= {args.min_breadth}%")
    print()
    print(f"  {'group':8} {'n':>4} {'mean richness':>14} {'mean coverage':>14}")
    print("  " + "-" * 44)
    for name, grp in (("HPV+", pos), ("HPV-", neg)):
        print(f"  {name:8} {len(grp):>4} "
              f"{mean([len(present[s]) for s in grp]):>14.1f} "
              f"{mean([sum(d for d, _ in cov[s].values()) for s in grp]):>14.1f}")
    print()
    print("  NOTE: HPV+ was simulated at 5-20x vs HPV- flat 10x, so richness is")
    print("  confounded with depth. Rarefy or model depth before interpreting.")
    print(f"\n  wrote 4 files to {outdir}/")


if __name__ == "__main__":
    main()
