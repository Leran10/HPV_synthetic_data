#!/usr/bin/env python3
"""Phage differential-abundance markers between HPV+ and HPV- in the jackalope cohort.

TWO ANALYSES
  1. Case/control   per-vOTU HPV+ (n=48) vs HPV- (n=52)
  2. Severity       per-vOTU abundance ~ Antler_Size_mm, within HPV+ only (n=48)

A vOTU significant in BOTH is the strongest marker candidate: present/abundant in
cases AND tracking disease severity.

THE DEPTH CONFOUND, AND WHY TWO PASSES
HPV+ samples were simulated at 5-20x coverage and HPV- at a flat 10x, so cases
carry ~1.5x the sequencing depth by construction. Detection rises with depth, so
a naive group test partly measures depth. Measured in this cohort:

    richness ~ coverage   all samples  r = +0.31
                          within HPV+  r = +0.39
                          within HPV-  r = +0.03

Coverage ranges overlap heavily (HPV+ 3001-33507, HPV- 4368-30959), so the
confound is controllable rather than fatal. Every vOTU is therefore tested twice:

    RAREFIED  abundances scaled to a common coverage, then rank-based test
    COVARIATE rank-based partial association, adjusting for log(coverage)

Markers that pass both are reported as robust; those that pass only the raw test
are reported separately and flagged as depth-suspect, not quietly dropped.

The severity analysis needs no such correction -- it is within-group, and uses
Spearman, which is invariant to the monotone depth scaling.

OUTPUTS
  markers_case_control.tsv   every vOTU, both passes, with annotation
  markers_severity.tsv       every vOTU, severity correlation, with annotation
  markers_consensus.tsv      vOTUs significant in both analyses
  analysis_summary.txt       counts, thresholds, caveats

Usage:
  python3 differential_abundance.py --indir results/vlp \
      --metadata jackalope_clinical_metadata.csv --outdir results/markers
"""
import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

# Detection thresholds, matching how the matrices were built.
PREVALENCE_MIN = 0.10   # a vOTU must be present in >=10% of samples to be testable
FDR_ALPHA = 0.05


def read_metadata(path):
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            size = row.get("Antler_Size_mm", "NA")
            out[row["Sample_ID"]] = {
                "group": "HPV+" if row.get("Antler_Present", "").strip() == "Yes" else "HPV-",
                "size": None if size in ("NA", "") else float(size),
                "age": float(row["Age_Years"]) if row.get("Age_Years") else None,
            }
    return out


def read_matrix(path):
    """Return (vOTU list, sample list, 2-D array)."""
    with open(path) as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        samples = header[1:]
        votus, rows = [], []
        for line in r:
            votus.append(line[0])
            rows.append([float(x) for x in line[1:]])
    return votus, samples, np.array(rows)


def read_annotations(path):
    out = {}
    if not Path(path).exists():
        return out
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["vOTU"]] = row
    return out


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p downwards
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def cliffs_delta(a, b):
    """Non-parametric effect size in [-1, 1]; sign says which group is higher."""
    a, b = np.asarray(a), np.asarray(b)
    if a.size == 0 or b.size == 0:
        return float("nan")
    gt = sum((a[:, None] > b[None, :]).sum(axis=1))
    lt = sum((a[:, None] < b[None, :]).sum(axis=1))
    return (gt - lt) / (a.size * b.size)


def partial_spearman(x, y, z):
    """Spearman association of x and y, adjusting for z.

    Rank all three, regress x and y on z separately, correlate the residuals.
    Rank-based so it does not assume the depth relationship is linear in counts.
    """
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    Z = np.column_stack([np.ones_like(rz), rz])

    def resid(v):
        beta, *_ = np.linalg.lstsq(Z, v, rcond=None)
        return v - Z @ beta

    ex, ey = resid(rx), resid(ry)
    if np.std(ex) == 0 or np.std(ey) == 0:
        return float("nan"), float("nan")
    r, _ = stats.pearsonr(ex, ey)
    # t-test on the partial correlation, df = n - 3 (two means + one covariate)
    n = len(x)
    df = n - 3
    if df <= 0 or abs(r) >= 1:
        return r, float("nan")
    t = r * math.sqrt(df / (1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), df)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, help="dir with votu_*_matrix.tsv")
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--outdir", default="results/markers")
    ap.add_argument("--prevalence", type=float, default=PREVALENCE_MIN)
    ap.add_argument("--alpha", type=float, default=FDR_ALPHA)
    args = ap.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = read_metadata(args.metadata)
    votus, samples, abund = read_matrix(indir / "votu_abundance_matrix.tsv")
    _, _, presence = read_matrix(indir / "votu_presence_matrix.tsv")
    annot = read_annotations(indir / "votu_annotations.tsv")

    missing = [s for s in samples if s not in meta]
    if missing:
        sys.exit(f"error: {len(missing)} samples absent from metadata, e.g. {missing[:3]}")

    groups = np.array([meta[s]["group"] for s in samples])
    is_pos = groups == "HPV+"
    is_neg = groups == "HPV-"
    depth = abund.sum(axis=0)                     # total coverage per sample
    log_depth = np.log10(np.maximum(depth, 1))

    # Rarefy by scaling each sample to the minimum total coverage. Proportional
    # scaling rather than random subsampling: these are coverage depths, not
    # integer read counts, so multinomial resampling is not well defined.
    scale = depth.min() / np.maximum(depth, 1e-9)
    rarefied = abund * scale

    # Only test vOTUs seen in enough samples; rare ones have no power and
    # inflate the multiple-testing burden.
    prev = presence.sum(axis=1) / len(samples)
    testable = prev >= args.prevalence

    # ---------------------------------------------------------- 1. case/control
    rows_cc = []
    for i, v in enumerate(votus):
        if not testable[i]:
            continue
        raw_p, raw_r = abund[i][is_pos], abund[i][is_neg]
        rar_p, rar_r = rarefied[i][is_pos], rarefied[i][is_neg]

        # raw and rarefied group tests
        try:
            _, p_raw = stats.mannwhitneyu(raw_p, raw_r, alternative="two-sided")
        except ValueError:
            p_raw = 1.0
        try:
            _, p_rar = stats.mannwhitneyu(rar_p, rar_r, alternative="two-sided")
        except ValueError:
            p_rar = 1.0

        # depth-adjusted: association of abundance with group, controlling depth
        r_cov, p_cov = partial_spearman(abund[i], is_pos.astype(float), log_depth)

        # Presence pattern matters for interpretation. The two arms were built
        # from DISJOINT ViroForge collections (24 for cases, 16 for controls), so
        # a phage in one collection and not the other separates the groups
        # perfectly by construction -- delta = 1.00, p ~ 1e-20, no biology.
        # Those are a positive control that the pipeline detects presence/absence.
        # Only vOTUs present in both arms pose a question whose answer could have
        # gone either way, so that is where the abundance test is informative
        # about realistic signal.
        npos = int(presence[i][is_pos].sum())
        nneg = int(presence[i][is_neg].sum())
        if npos > 0 and nneg == 0:
            pattern = "HPV+_exclusive"
        elif nneg > 0 and npos == 0:
            pattern = "HPV-_exclusive"
        elif npos >= 5 and nneg >= 5:
            pattern = "shared_both"
        else:
            pattern = "shared_rare_in_one"

        rows_cc.append({
            "vOTU": v,
            "presence_pattern": pattern,
            "prevalence": f"{prev[i]:.3f}",
            "n_present_pos": npos,
            "n_present_neg": nneg,
            "mean_abund_pos": f"{raw_p.mean():.4f}",
            "mean_abund_neg": f"{raw_r.mean():.4f}",
            "log2fc": (f"{math.log2((raw_p.mean() + 0.01) / (raw_r.mean() + 0.01)):.4f}"),
            "cliffs_delta": f"{cliffs_delta(raw_p, raw_r):.4f}",
            "p_raw": p_raw,
            "p_rarefied": p_rar,
            "partial_rho_depthadj": f"{r_cov:.4f}" if r_cov == r_cov else "NA",
            "p_depthadj": p_cov if p_cov == p_cov else 1.0,
        })

    for key, newkey in (("p_raw", "fdr_raw"), ("p_rarefied", "fdr_rarefied"),
                        ("p_depthadj", "fdr_depthadj")):
        adj = bh_fdr([r[key] for r in rows_cc])
        for r, a in zip(rows_cc, adj):
            r[newkey] = a

    for r in rows_cc:
        sig_rar = r["fdr_rarefied"] < args.alpha
        sig_cov = r["fdr_depthadj"] < args.alpha
        r["verdict"] = ("robust" if (sig_rar and sig_cov)
                        else "depth_suspect" if (r["fdr_raw"] < args.alpha)
                        else "ns")

    # ---------------------------------------------------------- 2. severity
    pos_idx = np.where(is_pos)[0]
    sizes = np.array([meta[samples[j]]["size"] for j in pos_idx], dtype=float)
    pos_logdepth = log_depth[pos_idx]
    rows_sev = []
    for i, v in enumerate(votus):
        if not testable[i]:
            continue
        vals = abund[i][pos_idx]
        if np.std(vals) == 0:
            rho, p = float("nan"), 1.0
            rho_adj, p_adj = float("nan"), 1.0
        else:
            rho, p = stats.spearmanr(vals, sizes)
            # Antler size DROVE simulated coverage (5x at 8mm -> 20x at 63mm),
            # so a raw size correlation partly measures sequencing depth. Partial
            # out log(coverage) to see what tracks severity independently of it.
            rho_adj, p_adj = partial_spearman(vals, sizes, pos_logdepth)
        rows_sev.append({
            "vOTU": v,
            "n_present_pos": int(presence[i][is_pos].sum()),
            "spearman_rho_size": f"{rho:.4f}" if rho == rho else "NA",
            "p_severity": p if p == p else 1.0,
            "partial_rho_size_depthadj": f"{rho_adj:.4f}" if rho_adj == rho_adj else "NA",
            "p_severity_depthadj": p_adj if p_adj == p_adj else 1.0,
        })
    for key, newkey in (("p_severity", "fdr_severity"),
                        ("p_severity_depthadj", "fdr_severity_depthadj")):
        adj = bh_fdr([r[key] for r in rows_sev])
        for r, a in zip(rows_sev, adj):
            r[newkey] = a

    # ---------------------------------------------------------- write
    ANN = ["tax_Lineage", "tax_Genus", "host_Host", "host_CHERRYScore",
           "life_lifestyle", "checkv_checkv_quality", "checkv_completeness",
           "checkv_contig_length"]

    def with_annot(row):
        a = annot.get(row["vOTU"], {})
        out = dict(row)
        out.update({k: a.get(k, "") for k in ANN})
        return out

    cc_cols = (["vOTU", "verdict", "presence_pattern", "prevalence",
                "n_present_pos", "n_present_neg",
                "mean_abund_pos", "mean_abund_neg", "log2fc", "cliffs_delta",
                "p_raw", "fdr_raw", "p_rarefied", "fdr_rarefied",
                "partial_rho_depthadj", "p_depthadj", "fdr_depthadj"] + ANN)
    rows_cc.sort(key=lambda r: r["fdr_rarefied"])
    with open(outdir / "markers_case_control.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cc_cols, delimiter="\t", restval="")
        w.writeheader()
        for r in rows_cc:
            w.writerow(with_annot(r))

    sev_cols = ["vOTU", "n_present_pos", "spearman_rho_size",
                "p_severity", "fdr_severity",
                "partial_rho_size_depthadj", "p_severity_depthadj",
                "fdr_severity_depthadj"] + ANN
    rows_sev.sort(key=lambda r: r["fdr_severity"])
    with open(outdir / "markers_severity.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sev_cols, delimiter="\t", restval="")
        w.writeheader()
        for r in rows_sev:
            w.writerow(with_annot(r))

    # consensus: robust in case/control AND significant on severity
    sev_by_v = {r["vOTU"]: r for r in rows_sev}
    consensus = [r for r in rows_cc
                 if r["verdict"] == "robust"
                 and sev_by_v.get(r["vOTU"], {}).get("fdr_severity", 1) < args.alpha]
    with open(outdir / "markers_consensus.tsv", "w", newline="") as fh:
        cols = ["vOTU", "cliffs_delta", "log2fc", "fdr_rarefied", "fdr_depthadj",
                "spearman_rho_size", "fdr_severity"] + ANN
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        for r in consensus:
            row = with_annot(r)
            s = sev_by_v[r["vOTU"]]
            row["spearman_rho_size"] = s["spearman_rho_size"]
            row["fdr_severity"] = s["fdr_severity"]
            w.writerow({k: row.get(k, "") for k in cols})

    # ---------------------------------------------------------- summary
    n_robust = sum(1 for r in rows_cc if r["verdict"] == "robust")
    n_suspect = sum(1 for r in rows_cc if r["verdict"] == "depth_suspect")
    n_sev = sum(1 for r in rows_sev if r["fdr_severity"] < args.alpha)
    n_sev_adj = sum(1 for r in rows_sev if r["fdr_severity_depthadj"] < args.alpha)
    from collections import Counter
    pat = Counter(r["presence_pattern"] for r in rows_cc if r["verdict"] == "robust")
    n_shared_rob = pat.get("shared_both", 0)
    up = sum(1 for r in rows_cc if r["verdict"] == "robust" and float(r["cliffs_delta"]) > 0)
    dn = n_robust - up

    lines = [
        "=" * 70,
        " PHAGE DIFFERENTIAL MARKERS -- HPV+ vs HPV-",
        "=" * 70,
        "",
        f"  vOTUs total              : {len(votus)}",
        f"  vOTUs testable           : {int(testable.sum())}  (present in >={args.prevalence:.0%} of samples)",
        f"  samples                  : {len(samples)}  (HPV+ {int(is_pos.sum())}, HPV- {int(is_neg.sum())})",
        f"  FDR threshold            : {args.alpha}",
        "",
        "  ANALYSIS 1 -- case/control",
        f"    robust markers         : {n_robust}   (significant rarefied AND depth-adjusted)",
        f"      enriched in HPV+     : {up}",
        f"      enriched in HPV-     : {dn}",
        f"    depth-suspect          : {n_suspect}   (raw only -- do not trust)",
        "",
        "    robust markers by presence pattern:",
        f"      HPV+ exclusive       : {pat.get('HPV+_exclusive', 0):4}   <- collection artefact",
        f"      HPV- exclusive       : {pat.get('HPV-_exclusive', 0):4}   <- collection artefact",
        f"      shared, rare in one  : {pat.get('shared_rare_in_one', 0):4}",
        f"      SHARED in both arms  : {n_shared_rob:4}   <- the informative subset",
        "",
        "  ANALYSIS 2 -- severity gradient (within HPV+, n=48)",
        f"    vOTUs tracking size    : {n_sev}",
        f"    still significant after adjusting for coverage : {n_sev_adj}",
        "",
        f"  CONSENSUS (both analyses): {len(consensus)}",
        "",
        "-" * 70,
        "  CAVEATS",
        "  * HPV+ was simulated at 5-20x and HPV- at flat 10x, so cases carry",
        "    ~1.5x depth by construction. 'robust' means the marker survived both",
        "    rarefaction and depth adjustment; 'depth_suspect' means it did not.",
        "  * Richness tracks antler size at r=0.95 in this cohort -- far tighter",
        "    than real biology. Severity drove simulated coverage, so treat the",
        "    severity analysis as a positive control that the pipeline recovers a",
        "    planted signal, not as a biological dose-response.",
        "  * MMseqs2 nr taxonomy was omitted (2.1 TB index on a 7-10 MB/s mount).",
        "    Taxonomy here is PhaGCN; hosts are CHERRY. Both are single-caller.",
        "  * Synthetic data: absolute effect sizes are not transferable to real",
        "    cervicovaginal samples.",
        "-" * 70,
    ]
    text = "\n".join(lines)
    (outdir / "analysis_summary.txt").write_text(text + "\n")
    print(text)
    print(f"\n  wrote 4 files to {outdir}/")


if __name__ == "__main__":
    main()
