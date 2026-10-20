#!/usr/bin/env python3
"""
Regenerate the five figures in docs/results_summary.html as standalone SVGs.

The report embeds these as inline SVG; this script is the generator of record,
so the figures can be rebuilt from the tracked result tables rather than edited
by hand.

  fig1_composition.svg   marker classification -- collection artefact vs informative
  fig2_volcano.svg       effect size vs depth-adjusted significance
  fig3_lactobacillus.svg Lactobacillus-host vOTUs, per-vOTU log2FC
  fig4_host_genus.svg    effect size by predicted host genus (strip plot)
  fig5_richness.svg      vOTU richness vs sequencing depth, by group

Inputs (all tracked in this repo):
  results/markers/markers_case_control.tsv   full 596-row test set
  results/vlp/votu_richness.tsv              per-sample richness and depth

Note: use markers_case_control.tsv, not markers_consensus.tsv. The consensus
file is already filtered (402 rows) and lacks the `verdict` and
`presence_pattern` columns the classification depends on.

Usage:
  python scripts/make_figures.py                       # writes to docs/figures/
  python scripts/make_figures.py --outdir /tmp/figs
  python scripts/make_figures.py --markers X --richness Y --outdir Z

No third-party plotting library: SVG is emitted directly so the output matches
the report byte-for-byte in style and needs no matplotlib install.
"""
from __future__ import annotations

import argparse
import csv
import html
import math
from pathlib import Path

# ---------------------------------------------------------------- palette ----
# Diverging pair + neutrals, validated for CVD separation (see dataviz skill).
DEPLETED = "#2a78d6"   # cool pole -- lower in HPV+
ENRICHED = "#e34948"   # warm pole -- higher in HPV+
ARTEFACT = "#b9b6ae"   # grey -- artefact, not a finding
NEUTRAL = "#f0efec"
INK = "#0b0b0b"
INK2 = "#52514e"
INK3 = "#7a7873"
RULE = "#e2e0da"
SANS = ("-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, "
        "'Helvetica Neue', Arial, sans-serif")
MONO = "ui-monospace, 'SF Mono', SFMono-Regular, Menlo, Consolas, monospace"


def esc(s) -> str:
    return html.escape(str(s), quote=True)


class SVG:
    """Minimal SVG builder. Theme-aware via a CSS block using currentColor."""

    def __init__(self, width: int, height: int, title: str):
        self.w, self.h = width, height
        self.parts: list[str] = []
        self.title = title

    def add(self, s: str) -> None:
        self.parts.append(s)

    def text(self, x, y, s, size=11, fill=INK2, anchor="start",
             family=SANS, weight="normal", opacity=1.0):
        self.add(
            f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" '
            f'text-anchor="{anchor}" font-family="{family}" '
            f'font-weight="{weight}" opacity="{opacity:g}">{esc(s)}</text>'
        )

    def rect(self, x, y, w, h, fill, rx=0, opacity=1.0, stroke=None):
        sk = f' stroke="{stroke}" stroke-width="1"' if stroke else ""
        self.add(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w,0):.1f}" '
            f'height="{max(h,0):.1f}" fill="{fill}" rx="{rx}" '
            f'opacity="{opacity:g}"{sk}/>'
        )

    def line(self, x1, y1, x2, y2, stroke=RULE, width=1, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="{width}"{d}/>'
        )

    def circle(self, cx, cy, r, fill, opacity=1.0, stroke=None, sw=1.5):
        sk = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
        self.add(
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{fill}" '
            f'opacity="{opacity:g}"{sk}/>'
        )

    def tip(self, text: str) -> None:
        """Native tooltip on the preceding element group."""
        self.add(f"<title>{esc(text)}</title>")

    def render(self) -> str:
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w} '
            f'{self.h}" width="100%" role="img" aria-label="{esc(self.title)}">'
            f"<title>{esc(self.title)}</title>"
            + "".join(self.parts)
            + "</svg>"
        )

    def write(self, path: Path) -> None:
        path.write_text(self.render(), encoding="utf-8")
        print(f"  wrote {path}  ({path.stat().st_size:,} bytes)")


# ------------------------------------------------------------------- data ----
def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(v, default=float("nan")) -> float:
    try:
        f = float(v)
        return f if math.isfinite(f) else default
    except (TypeError, ValueError):
        return default


def classify(r: dict) -> str:
    """Three-way split, straight from the analysis columns.

    `verdict` is 'robust' only when a vOTU clears FDR in BOTH the rarefied and
    the depth-adjusted test (set by differential_abundance.py). But a robust
    marker present in only ONE arm is a collection artefact: the two ViroForge
    collections do not share those genomes, so absence is by construction, not
    biology. Only robust markers present in BOTH arms are informative.
    """
    verdict = (r.get("verdict") or "").strip()
    pattern = (r.get("presence_pattern") or "").strip()
    if verdict != "robust":
        return "ns"
    if pattern in {"HPV+_exclusive", "HPV-_exclusive"}:
        return "artefact"
    if pattern == "shared_both":
        return "informative"
    return "borderline"       # shared_rare_in_one


def genus_of(r: dict) -> str:
    g = (r.get("host_Host") or "").strip()
    if not g or g.lower() in {"na", "none", "-", "unassigned"}:
        return "unassigned"
    return g.split(";")[-1].replace("g__", "").strip() or "unassigned"


# ------------------------------------------------------------------ fig 1 ----
def fig1(rows, out: Path):
    counts = {"informative": 0, "borderline": 0, "artefact": 0, "ns": 0}
    for r in rows:
        counts[classify(r)] += 1
    total = sum(counts.values()) or 1

    W, H = 760, 190
    s = SVG(W, H, "Marker classification after depth control")
    s.text(0, 16, "Markers surviving depth control", 13, INK, weight="600")
    s.text(0, 34, f"{total:,} vOTUs tested. Informative = clears both depth "
                  f"tests AND is present in both arms.", 11, INK3)

    bar_y, bar_h, bar_w = 56, 34, W
    order = [
        ("artefact", ARTEFACT, "Collection artefact"),
        ("borderline", "#d8d5cc", "Rare in one arm"),
        ("informative", ENRICHED, "Informative"),
        ("ns", NEUTRAL, "Not significant"),
    ]
    x = 0.0
    for key, colour, label in order:
        w = bar_w * counts[key] / total
        s.add("<g>")
        s.rect(x, bar_y, max(w - 2, 0), bar_h, colour, rx=2)
        s.tip(f"{label}: {counts[key]:,} vOTUs ({100*counts[key]/total:.1f}%)")
        s.add("</g>")
        if w > 44:
            s.text(x + 8, bar_y + 22, f"{counts[key]:,}", 12,
                   "#ffffff" if key in {"artefact", "informative"} else INK,
                   weight="600", family=MONO)
        x += w

    ly = bar_y + bar_h + 28
    lx = 0.0
    for key, colour, label in order:
        s.rect(lx, ly - 9, 11, 11, colour, rx=2, stroke=RULE if key == "ns" else None)
        s.text(lx + 17, ly, f"{label} — {counts[key]:,}", 11, INK2)
        lx += 185
    s.text(0, H - 8,
           "A marker present in only one arm cannot be biology: the two source "
           "collections do not share those genomes, so absence is by design.",
           10, INK3)
    s.write(out / "fig1_composition.svg")
    return counts


# ------------------------------------------------------------------ fig 2 ----
def fig2(rows, out: Path):
    pts = []
    for r in rows:
        d = fnum(r.get("cliffs_delta"))
        q = fnum(r.get("fdr_depthadj"))
        if not (math.isfinite(d) and math.isfinite(q)):
            continue
        y = -math.log10(max(q, 1e-300))
        pts.append((d, y, classify(r), r.get("vOTU", ""), genus_of(r)))
    if not pts:
        return

    W, H = 760, 420
    L, R, T, B = 54, 14, 54, 46
    pw, ph = W - L - R, H - T - B
    ymax = max(p[1] for p in pts) * 1.06 or 1

    s = SVG(W, H, "Effect size versus depth-adjusted significance")
    s.text(0, 16, "Effect size vs depth-adjusted significance", 13, INK,
           weight="600")
    s.text(0, 34, "Coloured points are informative markers; grey are "
                  "collection artefacts or non-significant", 11, INK3)

    def X(d):
        return L + (d + 1) / 2 * pw

    def Y(v):
        return T + ph - (v / ymax) * ph

    for gv in (-1, -0.5, 0, 0.5, 1):
        s.line(X(gv), T, X(gv), T + ph, RULE, 1,
               dash=None if gv == 0 else "2,3")
        s.text(X(gv), T + ph + 18, f"{gv:g}", 10, INK3, anchor="middle",
               family=MONO)
    nsteps = 4
    for i in range(nsteps + 1):
        v = ymax * i / nsteps
        s.line(L, Y(v), L + pw, Y(v), RULE, 1, dash="2,3")
        s.text(L - 8, Y(v) + 3, f"{v:.0f}", 10, INK3, anchor="end", family=MONO)

    cut = -math.log10(0.05)
    if cut < ymax:
        s.line(L, Y(cut), L + pw, Y(cut), ENRICHED, 1, dash="4,3")
        s.text(L + pw, Y(cut) - 6, "FDR 0.05", 9, ENRICHED, anchor="end",
               family=MONO)

    for d, y, cls, name, genus in sorted(pts, key=lambda p: p[2] == "informative"):
        if cls == "informative":
            col, rad, op = (ENRICHED if d > 0 else DEPLETED), 4.6, 0.95
        elif cls in {"artefact", "borderline"}:
            col, rad, op = ARTEFACT, 3.4, 0.6
        else:
            col, rad, op = ARTEFACT, 2.6, 0.32
        s.add("<g>")
        s.circle(X(d), Y(y), rad, col, op,
                 stroke="#ffffff" if cls == "informative" else None, sw=1.2)
        s.tip(f"{name}\nhost: {genus}\nCliff's delta {d:+.2f}, "
              f"FDR {10**-y:.2e}")
        s.add("</g>")

    s.text(L + pw / 2, H - 20, "← lower in HPV+  ·  Cliff's delta  ·  "
                               "higher in HPV+ →", 10, INK3, anchor="middle")
    s.add(f'<text transform="translate(13,{T+ph/2}) rotate(-90)" '
          f'font-size="10" fill="{INK3}" text-anchor="middle" '
          f'font-family="{SANS}">−log10 FDR (depth-adjusted)</text>')
    s.write(out / "fig2_volcano.svg")


# ------------------------------------------------------------------ fig 3 ----
def fig3(rows, out: Path, genus_key="Lactobacillus"):
    # Informative markers only. Including the collection artefacts here would
    # add ~150 rows that are absent-by-design and bury the actual pattern.
    sel = [r for r in rows
           if genus_key.lower() in genus_of(r).lower()
           and classify(r) == "informative"
           and math.isfinite(fnum(r.get("log2fc")))]
    sel.sort(key=lambda r: fnum(r["log2fc"]))
    if not sel:
        print(f"  (no informative {genus_key} vOTUs; skipping fig3)")
        return
    n_all = sum(1 for r in rows if genus_key.lower() in genus_of(r).lower())

    rowh = 15
    W = 760
    T, B, L = 62, 40, 150
    H = T + rowh * len(sel) + B
    vmax = max(abs(fnum(r["log2fc"])) for r in sel) * 1.1 or 1
    pw = W - L - 20

    down = sum(1 for r in sel if fnum(r["log2fc"]) < 0)
    s = SVG(W, H, f"{genus_key}-host vOTUs")
    s.text(0, 16, f"{genus_key}-host phages: all shift one way", 13, INK,
           weight="600")
    s.text(0, 34, f"{down} of {len(sel)} informative {genus_key} vOTUs are "
                  f"depleted in HPV+ (of {n_all} total; the rest are "
                  f"collection artefacts)", 11, INK3)

    def X(v):
        return L + pw / 2 + (v / vmax) * (pw / 2)

    s.line(X(0), T - 6, X(0), H - B + 4, INK3, 1)
    for gv in (-vmax / 2, vmax / 2):
        s.line(X(gv), T - 6, X(gv), H - B + 4, RULE, 1, dash="2,3")

    for i, r in enumerate(sel):
        v = fnum(r["log2fc"])
        y = T + i * rowh
        col = ENRICHED if v > 0 else DEPLETED
        op = 0.9
        x0, x1 = (X(0), X(v)) if v > 0 else (X(v), X(0))
        s.add("<g>")
        s.rect(x0, y + 2, x1 - x0, rowh - 5, col, rx=2, opacity=op)
        s.tip(f"{r['vOTU']}\nlog2FC {v:+.2f}\n"
              f"Cliff's delta {fnum(r.get('cliffs_delta')):+.2f}\n"
              f"host: {genus_of(r)}")
        s.add("</g>")
        name = r["vOTU"]
        s.text(L - 8, y + rowh - 4, name[-22:], 8.5, INK3, anchor="end",
               family=MONO)

    s.text(L + pw / 2, H - 14, "← depleted in HPV+  ·  log2 fold change  ·  "
                               "enriched in HPV+ →", 10, INK3, anchor="middle")
    s.write(out / "fig3_lactobacillus.svg")


# ------------------------------------------------------------------ fig 4 ----
def fig4(rows, out: Path, min_n=3):
    by: dict[str, list[float]] = {}
    for r in rows:
        d = fnum(r.get("cliffs_delta"))
        if math.isfinite(d):
            by.setdefault(genus_of(r), []).append(d)
    keep = {g: v for g, v in by.items() if len(v) >= min_n and g != "unassigned"}
    if not keep:
        print("  (no host genera with enough vOTUs; skipping fig4)")
        return

    def med(v):
        v = sorted(v)
        n = len(v)
        return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2

    order = sorted(keep, key=lambda g: med(keep[g]))
    rowh = 26
    W, L, T, B = 760, 150, 60, 42
    H = T + rowh * len(order) + B
    pw = W - L - 20

    s = SVG(W, H, "Effect size by predicted host genus")
    s.text(0, 16, "Effect size by predicted host genus", 13, INK, weight="600")
    s.text(0, 34, f"Genera with ≥{min_n} vOTUs; one dot per vOTU, bar marks the "
                  f"median. Host calls cover only part of the set — see caveat.",
           11, INK3)

    def X(d):
        return L + (d + 1) / 2 * pw

    s.line(X(0), T - 6, X(0), H - B + 2, INK3, 1)
    for gv in (-1, -0.5, 0.5, 1):
        s.line(X(gv), T - 6, X(gv), H - B + 2, RULE, 1, dash="2,3")
        s.text(X(gv), H - B + 18, f"{gv:g}", 10, INK3, anchor="middle",
               family=MONO)

    for i, g in enumerate(order):
        vals = keep[g]
        y = T + i * rowh + rowh / 2
        if i % 2 == 0:
            s.rect(L, y - rowh / 2, pw, rowh, NEUTRAL, opacity=0.5)
        s.text(L - 10, y + 4, g, 11, INK2, anchor="end")
        s.text(L - 10, y + 15, f"n={len(vals)}", 8.5, INK3, anchor="end",
               family=MONO)
        m = med(vals)
        for d in vals:
            s.add("<g>")
            s.circle(X(d), y, 4.2, ENRICHED if d > 0 else DEPLETED, 0.55,
                     stroke="#ffffff", sw=1.0)
            s.tip(f"{g}: Cliff's delta {d:+.2f}")
            s.add("</g>")
        s.line(X(m), y - 9, X(m), y + 9, INK, 2.2)

    s.text(L + pw / 2, H - 10, "← lower in HPV+  ·  Cliff's delta  ·  "
                               "higher in HPV+ →", 10, INK3, anchor="middle")
    s.write(out / "fig4_host_genus.svg")


# ------------------------------------------------------------------ fig 5 ----
def fig5(rich, out: Path):
    pts = []
    for r in rich:
        x = fnum(r.get("total_coverage"))
        y = fnum(r.get("votu_richness"))
        if math.isfinite(x) and math.isfinite(y) and x > 0:
            pts.append((x, y, (r.get("group") or "").strip(),
                        r.get("sample", "")))
    if not pts:
        return

    W, H = 760, 400
    L, R, T, B = 56, 16, 58, 52
    pw, ph = W - L - R, H - T - B
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    x0, x1 = min(xs) * 0.9, max(xs) * 1.05
    y0, y1 = min(ys) * 0.92, max(ys) * 1.04

    s = SVG(W, H, "vOTU richness versus sequencing depth")
    s.text(0, 16, "Richness tracks sequencing depth", 13, INK, weight="600")
    s.text(0, 34, "Each point is a sample. The groups differ in depth, so "
                  "richness differences are confounded by design.", 11, INK3)

    def X(v):
        return L + (math.log10(v) - math.log10(x0)) / (
            math.log10(x1) - math.log10(x0)) * pw

    def Y(v):
        return T + ph - (v - y0) / (y1 - y0) * ph

    dec = int(math.floor(math.log10(x0)))
    while 10 ** dec <= x1:
        for m in (1, 2, 5):
            v = m * 10 ** dec
            if x0 <= v <= x1:
                s.line(X(v), T, X(v), T + ph, RULE, 1, dash="2,3")
                lab = f"{v/1000:g}k" if v >= 1000 else f"{v:g}"
                s.text(X(v), T + ph + 18, lab, 10, INK3, anchor="middle",
                       family=MONO)
        dec += 1
    for i in range(5):
        v = y0 + (y1 - y0) * i / 4
        s.line(L, Y(v), L + pw, Y(v), RULE, 1, dash="2,3")
        s.text(L - 8, Y(v) + 3, f"{v:.0f}", 10, INK3, anchor="end", family=MONO)

    groups = sorted({p[2] for p in pts})
    colour = {g: (ENRICHED if "+" in g else DEPLETED) for g in groups}
    for x, y, g, name in pts:
        s.add("<g>")
        s.circle(X(x), Y(y), 4.6, colour.get(g, ARTEFACT), 0.75,
                 stroke="#ffffff", sw=1.1)
        s.tip(f"{name} ({g})\ndepth {x:,.0f}×\nrichness {y:,.0f} vOTUs")
        s.add("</g>")

    lx = L
    for g in groups:
        n = sum(1 for p in pts if p[2] == g)
        s.circle(lx + 5, T - 14, 4.6, colour.get(g, ARTEFACT), 0.85)
        s.text(lx + 15, T - 10, f"{g} (n={n})", 11, INK2)
        lx += 150

    s.text(L + pw / 2, H - 20, "total coverage (log scale)", 10, INK3,
           anchor="middle")
    s.add(f'<text transform="translate(14,{T+ph/2}) rotate(-90)" '
          f'font-size="10" fill="{INK3}" text-anchor="middle" '
          f'font-family="{SANS}">vOTU richness</text>')
    s.write(out / "fig5_richness.svg")


# ------------------------------------------------------------------- main ----
def main() -> None:
    root = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--markers", type=Path,
                    default=root / "results/markers/markers_case_control.tsv")
    ap.add_argument("--richness", type=Path,
                    default=root / "results/vlp/votu_richness.tsv")
    ap.add_argument("--outdir", type=Path, default=root / "docs/figures")
    ap.add_argument("--host-genus", default="Lactobacillus",
                    help="genus to highlight in figure 3")
    args = ap.parse_args()

    for p in (args.markers, args.richness):
        if not p.exists():
            raise SystemExit(f"error: missing input {p}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(args.markers)
    rich = read_tsv(args.richness)
    print(f"Read {len(rows):,} markers, {len(rich):,} samples")

    counts = fig1(rows, args.outdir)
    fig2(rows, args.outdir)
    fig3(rows, args.outdir, args.host_genus)
    fig4(rows, args.outdir)
    fig5(rich, args.outdir)

    print(f"\nClassification: {counts['informative']:,} informative, "
          f"{counts['artefact']:,} collection artefact, "
          f"{counts['borderline']:,} rare-in-one-arm, {counts['ns']:,} ns")
    print(f"Figures in {args.outdir}")


if __name__ == "__main__":
    main()
