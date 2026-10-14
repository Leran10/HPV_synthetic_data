#!/usr/bin/env python3
"""Expand a Hecatomb preprocessing run into its real QC sub-steps.

The `--preprocess paired` route runs ONE fastp call that performs several
distinct filtering operations, then host removal, then dereplication. The
top-level attrition table shows only three rows, which understates what the QC
actually did -- fastp reports each filter separately in its JSON, so those
sub-steps can be recovered.

Steps reported (reads = single reads, i.e. both mates, matching fastp's own
accounting; the top-level attrition table counts PAIRS):

  1a  adapter trimming        reads with adapter removed (trimmed, not discarded)
  1b  poly-X trimming         reads with poly-X tail removed (trimmed, not discarded)
  1c  duplicate removal       reads dropped as PCR/optical duplicates
  1d  low-quality filter      reads dropped, mean quality below Q15
  1e  too-short filter        reads dropped, <90bp after trimming
  1f  too-many-N filter       reads dropped, >2 ambiguous bases
  2   host removal            pairs dropped, mapped to viral-masked GRCh38
  3   dereplication           reads collapsed onto representatives (lossless)

Usage:
    python3 qc_steps_detailed.py QC_DIR METADATA_CSV -o results/qc/ --label bulk
"""
import argparse
import csv
import json
import sys
from pathlib import Path


TRIM_OPS = {"adapter", "polyx"}          # modify reads, do not discard them


def read_groups(path):
    groups = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            groups[row["Sample_ID"]] = "HPV+" if row.get(
                "Antler_Present", "").strip() == "Yes" else "HPV-"
    return groups


def fastp_detail(stats_dir, sample):
    """Per-filter counts from one sample's fastp JSON (counts are single reads)."""
    f = stats_dir / "p01" / f"{sample}.s1.stats.json"
    if not f.exists():
        return None
    with open(f) as fh:
        d = json.load(fh)

    summ = d.get("summary", {})
    before = summ.get("before_filtering", {}).get("total_reads", 0)
    after = summ.get("after_filtering", {}).get("total_reads", 0)
    filt = d.get("filtering_result", {})
    passed = filt.get("passed_filter_reads", 0)

    # fastp applies dedup separately from the pass/fail filters: anything that
    # passed the filters but is absent from the output was dropped as duplicate.
    dup_removed = max(0, passed - after)

    return {
        "raw": before,
        "adapter_trimmed": d.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0),
        "polyx_trimmed": d.get("polyx_trimming", {}).get("total_polyx_trimmed_reads", 0),
        "dup_removed": dup_removed,
        "low_quality": filt.get("low_quality_reads", 0),
        "too_short": filt.get("too_short_reads", 0),
        "too_many_n": filt.get("too_many_N_reads", 0),
        "after_fastp": after,
        "dup_rate": d.get("duplication", {}).get("rate", 0.0),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("qc_dir")
    ap.add_argument("metadata")
    ap.add_argument("-o", "--outdir", default=".")
    ap.add_argument("--label", default="cohort")
    ap.add_argument("--attrition", help="qc_attrition*.tsv, for the host/cluster rows")
    args = ap.parse_args()

    qc = Path(args.qc_dir)
    stats = qc / "processing" / "stats"
    if not stats.exists():
        sys.exit(f"error: no stats directory under {qc}")

    groups = read_groups(args.metadata)

    # ---- per-sample fastp detail
    per_sample = {}
    for s in sorted(groups):
        d = fastp_detail(stats, s)
        if d:
            per_sample[s] = d
    if not per_sample:
        sys.exit("error: no fastp JSON found")

    def total(key, group=None):
        return sum(v[key] for s, v in per_sample.items()
                   if group is None or groups[s] == group)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- host + cluster rows come from the pair-level attrition table
    host_pairs_in = host_pairs_out = clust_out = reps = None
    if args.attrition and Path(args.attrition).exists():
        with open(args.attrition) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))

        def isum(k):
            return sum(int(r[k]) for r in rows if r.get(k, "NA") not in ("NA", ""))
        host_pairs_in = isum("after_fastp")
        host_pairs_out = isum("after_host")
        clust_out = isum("after_cluster")
        reps = isum("n_representatives")

    # ---- build the step table
    #
    # Every row carries BOTH units so the funnel reads top to bottom without
    # mental arithmetic. fastp counts single reads (both mates), host removal
    # and clustering count pairs; `native_unit` records which one the tool
    # actually reported, the other is derived (pairs = reads / 2).
    raw = total("raw")
    step_cols = ["step", "operation", "software", "version", "native_unit",
                 "reads_in", "pairs_in",
                 "reads_affected", "pairs_affected",
                 "reads_discarded", "pairs_discarded",
                 "pct_discarded", "sequences_out", "compression", "type"]
    S = []

    def _both(n, unit):
        """Return (reads, pairs) given a count in `unit`."""
        if not isinstance(n, int):
            return n, n
        return (n, n // 2) if unit == "reads" else (n * 2, n)

    def row(step, op, sw, ver, unit, rin, affected, discarded, typ,
            seqs_out=None):
        r_in, p_in = _both(rin, unit)
        r_af, p_af = _both(affected, unit)
        r_di, p_di = _both(discarded, unit)

        def f(x):
            return f"{x:,}" if isinstance(x, int) else x

        S.append({
            "step": step, "operation": op, "software": sw, "version": ver,
            "native_unit": unit,
            "reads_in": f(r_in), "pairs_in": f(p_in),
            "reads_affected": f(r_af), "pairs_affected": f(p_af),
            "reads_discarded": f(r_di), "pairs_discarded": f(p_di),
            "pct_discarded": (f"{100.0*discarded/rin:.3f}%"
                              if isinstance(discarded, int) and isinstance(rin, int)
                              and rin else "—"),
            # Only meaningful for the collapse step: how many distinct sequences
            # come out, and how much the search space shrank. Reads discarded is
            # the wrong metric there -- dereplication is lossless by design.
            "sequences_out": f(seqs_out) if seqs_out is not None else "—",
            "compression": (f"{rin/seqs_out:.1f}x"
                            if seqs_out and isinstance(rin, int) else "—"),
            "type": typ,
        })

    FP, FPV = "fastp", "0.23.2"
    after_fp = total("after_fastp")
    fp_discarded = raw - after_fp

    # Stage 1 is ONE fastp call. Its sub-operations are reported separately in
    # fastp's JSON but they are not sequential passes -- trimming happens first,
    # then each read is assigned to exactly one discard reason (the reasons are
    # mutually exclusive and sum to the total). So only the STAGE row carries a
    # reads_in/reads_out flow; the sub-rows are a breakdown of why reads went.
    row("0", "input reads", "—", "—", "reads", raw, "—", 0, "input")
    row("1", "fastp: trim + filter + deduplicate (single pass)", FP, FPV,
        "reads", raw, fp_discarded, fp_discarded, "stage")
    row("1 .1", "  └ trimmed: adapter (--detect_adapter_for_pe)", FP, FPV,
        "reads", "—", total("adapter_trimmed"), 0, "trim")
    row("1 .2", "  └ trimmed: poly-X tail (--trim_poly_x)", FP, FPV,
        "reads", "—", total("polyx_trimmed"), 0, "trim")
    row("1 .3", "  └ discarded: low quality (Q<15, cut_tail w=25)", FP, FPV,
        "reads", "—", total("low_quality"), total("low_quality"), "reason")
    row("1 .4", "  └ discarded: too short (<90bp)", FP, FPV,
        "reads", "—", total("too_short"), total("too_short"), "reason")
    row("1 .5", "  └ discarded: too many N (>2)", FP, FPV,
        "reads", "—", total("too_many_n"), total("too_many_n"), "reason")
    row("1 .6", "  └ discarded: PCR/optical duplicate (--dedup acc=4)", FP, FPV,
        "reads", "—", total("dup_removed"), total("dup_removed"), "reason")

    if host_pairs_in is not None:
        removed = host_pairs_in - host_pairs_out
        row("2", "host removal: map to viral-masked GRCh38, keep unmapped",
            "minimap2 + samtools", ">=2.20 / >=1.13", "pairs",
            host_pairs_in, removed, removed, "discard")
        # Dereplication runs on R1 only, so its input is reads-as-pairs and its
        # output is a count of distinct sequences -- NOT a reads->reads flow.
        # reads_in is blanked to stop it being read as "61M shrank to 5M".
        row("3", "dereplication at 97% id (R1 only; counts preserved)",
            "MMseqs2 linclust", "12.113e3", "pairs",
            host_pairs_out, clust_out, 0, "collapse", seqs_out=reps)
        S[-1]["reads_in"] = "—"
        S[-1]["reads_affected"] = "—"
        S[-1]["reads_discarded"] = "0"

    step_path = outdir / f"qc_steps_detailed_{args.label}.tsv"
    with open(step_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=step_cols, delimiter="\t")
        w.writeheader()
        w.writerows(S)

    # ---- per-group totals for the fastp sub-steps
    grp_path = outdir / f"qc_steps_by_group_{args.label}.tsv"
    gcols = ["group", "n_samples", "raw_reads", "raw_pairs",
             "adapter_trimmed", "polyx_trimmed",
             "low_quality_removed", "too_short_removed", "too_many_n_removed",
             "duplicates_removed", "mean_duplication_rate",
             "after_fastp_reads", "after_fastp_pairs"]
    with open(grp_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=gcols, delimiter="\t")
        w.writeheader()
        for g in sorted({v for v in groups.values()}) + [None]:
            names = [s for s in per_sample if g is None or groups[s] == g]
            if not names:
                continue
            dr = [per_sample[s]["dup_rate"] for s in names]
            raw_g = total("raw", g)
            fp_g = total("after_fastp", g)
            w.writerow({
                "group": g or "TOTAL", "n_samples": len(names),
                "raw_reads": f"{raw_g:,}", "raw_pairs": f"{raw_g // 2:,}",
                "adapter_trimmed": f"{total('adapter_trimmed', g):,}",
                "polyx_trimmed": f"{total('polyx_trimmed', g):,}",
                "low_quality_removed": f"{total('low_quality', g):,}",
                "too_short_removed": f"{total('too_short', g):,}",
                "too_many_n_removed": f"{total('too_many_n', g):,}",
                "duplicates_removed": f"{total('dup_removed', g):,}",
                "mean_duplication_rate": f"{sum(dr)/len(dr):.4f}",
                "after_fastp_reads": f"{fp_g:,}",
                "after_fastp_pairs": f"{fp_g // 2:,}",
            })

    # ---- print
    W = 118
    print(f"\n{'=' * W}\n QC STEPS IN DETAIL -- {args.label}\n{'=' * W}")
    print(f"{'step':6} {'operation':50} {'reads in':>14} "
          f"{'discarded':>14} {'reads out':>14}")
    print("-" * W)
    for r in S:
        # Only real stages have a flow; sub-rows show their contribution only.
        if r["type"] in ("stage", "input", "discard", "collapse"):
            ri = r["reads_in"]
            out = ("—" if r["reads_in"] == "—" else
                   f"{int(r['reads_in'].replace(',','')) - int(r['reads_discarded'].replace(',','')):,}")
        else:
            ri, out = "—", "—"
        print(f"{r['step']:6} {r['operation'][:49]:50} {ri:>14} "
              f"{r['reads_discarded']:>14} {out:>14}")
    print("-" * W)
    print("  Stage 1 is a SINGLE fastp call: the indented rows break down why reads")
    print("  went, not sequential passes. fastp assigns each discarded read exactly")
    print("  one reason, so the reasons sum to the stage total.")
    reasons = sum(int(r["reads_discarded"].replace(",", ""))
                  for r in S if r["type"] == "reason")
    stage1 = next((int(r["reads_discarded"].replace(",", ""))
                   for r in S if r["type"] == "stage"), 0)
    print(f"    check: reasons sum to {reasons:,} vs stage 1 total {stage1:,} "
          f"-> {'OK' if reasons == stage1 else 'MISMATCH'}")
    coll = [r for r in S if r["type"] == "collapse" and r["sequences_out"] != "—"]
    for r in coll:
        # Two different objects are being counted here, so spell it out:
        # reads are not reduced, the number of DISTINCT SEQUENCES stored is.
        print(f"  Stage {r['step']} is LOSSLESS. It reduces how many sequences you "
              f"search,")
        print(f"  not how many reads you counted:")
        print(f"    {r['sequences_out']} representative sequences cover all "
              f"{r['pairs_in']} reads")
        print(f"    (~{r['compression'].rstrip('x')} reads per representative; "
              f"every read's count is carried onto its representative, 0 lost)")
    print("  fastp reports single reads (both mates); host removal and clustering "
          "report pairs.")

    print(f"\nwrote {step_path}")
    print(f"wrote {grp_path}")


if __name__ == "__main__":
    main()
