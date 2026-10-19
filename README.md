# Jackalope cohort — phage community analysis by HPV status

**Question:** does the bacteriophage community differ between HPV-positive and
HPV-negative samples? HPV status is the *grouping variable*, not the readout —
this is not HPV detection.

## Headline result

**29 Lactobacillus-host phages are depleted in HPV+ samples** (median log2FC
−1.16), with no exceptions among the 29. The signal survived rarefaction and
depth adjustment, has moderate effect sizes (|Cliff's δ| ~0.54 rather than the
artefactual 1.00), and the host assignments are corroborated by two independent
callers (CHERRY and iPHoP). It is consistent with the *Lactobacillus*-depletion
pattern reported in the HPV literature.

| | |
|---|---|
| vOTUs recovered | 2,921 |
| testable (≥10% prevalence) | 596 |
| **informative markers** (present in both arms) | **74** |
| of those, *Lactobacillus* phages depleted in HPV+ | **29 / 29** |
| severity-tracking vOTUs after depth adjustment | 337 |

## The two cohorts

100 simulated vaginal virome samples generated with
[ViroForge](https://github.com/Leran10/viroforge) 0.20.0, on HTCF:

| | VLP (`jackalope_cohort`) | bulk (`jackalope_bulk`) |
|---|---|---|
| viral fraction | ~98% | ~1% |
| host DNA | 0.67% | ~40% |
| assembly N50 | ~30 kb | ~3 kb |
| longest contig | 193 kb (82% of a CMV genome) | 10 kb (fragment ceiling) |
| phages assemble? | **yes** | no |
| teaches host removal? | no (1.4% removed) | **yes (44% removed)** |

**The VLP cohort carries the analysis.** Bulk mode puts viruses at ~1% of reads
split across 25–54 genomes per sample; because ViroForge allocates abundance by
read share rather than genome size, a 50 kb phage received ~1× coverage and could
not assemble. Only 211 viral contigs were recovered from 1,287 genomes actually
present. VLP lifts the same phage to ~110×.

The bulk cohort is still the better QC teaching case — host removal discards 44%
of reads there versus 1.4% in VLP.

## Groups

`jackalope_clinical_metadata.csv` is a *blinded* clinical sheet: `Antler_Present`
is the case/control label (Yes = HPV+, n=48; No = HPV−, n=52) and
`Antler_Size_mm` (8–63) is a severity score, NA for all controls. One subject per
sample, no batches, no repeated measures. **There is no gender or other clinical
field.**

## ⚠ Two confounds that shape every result

**1. Depth is coupled to severity.** `jackalope_cohort.py` scales coverage with
antler size (5× at 8 mm → 20× at 63 mm) while all controls got a flat 10×. So
HPV+ averages ~1.5× the depth of HPV− by construction, and richness rises with
depth. Every marker is therefore tested three ways — raw, rarefied, and with
log(coverage) partialled out — and labelled `robust` / `depth_suspect` / `ns`.
Richness tracks antler size at r = 0.95, far tighter than real biology.

**2. The two arms were built from disjoint genome sets.** Cases came from
ViroForge collection 24, controls from collection 16. A phage in one collection
and not the other separates the groups *perfectly* (δ = 1.00, p ≈ 1e-20) with no
biology involved. **421 of 575 "robust" markers are this artefact**, now labelled
in the `presence_pattern` column. Only the 74 present in both arms pose a
question whose answer could have gone either way.

## Layout

```
config/                      project + Hecatomb QC configs
scripts/
  generate_bulk_cohort.sh    ViroForge, SLURM array, seeds from metadata
  assemble_vlp.sh            MEGAHIT per sample (strips spaces from headers)
  stage_reads*.sh            gzip + rename to VIR_###_R{1,2}.fastq.gz
  run_qc*.sh                 Hecatomb preprocessing
  collect_qc_attrition.py    per-step read attrition
  qc_summary.py              tool/version/reads-removed tables
  qc_steps_detailed.py       fastp sub-steps broken out
  run_prophage.sh            prophage-pipeline (abandoned, see below)
  run_phage_analysis_vlp.sh  phage-analysis, all resources pinned on the CLI
  votu_abundance_vlp.sh      bbmap each sample vs vOTU reps  <- the missing piece
  build_votu_matrix.py       joins abundance + taxonomy + host + lifestyle
  differential_abundance.py  the marker analysis
results/qc/                  QC attrition and step tables, both cohorts
results/vlp/                 vOTU × sample matrices + annotations
results/markers/             differential abundance output
docs/                        ViroForge fragment-length issue report
```

Data and pipeline outputs stay on HTCF under `/scratch/sahlab/leran/` and are
gitignored. A copy of the scripts and tables also lives in
`/scratch/sahlab/leran/hpv_project_files/`.

## Pipeline notes, hard-won

**prophage-pipeline is structurally incompatible with ViroForge output.** It bins
bacterial genomes to find prophages inside them, but ViroForge supplies bacterial
background as random 10 kb fragments — realism ballast, not whole genomes. DAS
Tool returned *"No bins with bin-score >0.5"* on 230–306 kb assemblies. No
parameter fixes this. See `docs/viroforge_fragment_bias_issue.md`.

**phage-analysis needed six fixes to run clean** (job 49642184, 435/435 steps,
3h51m, after 15 failed attempts):

| failure | cause | fix |
|---|---|---|
| 4 CPUs not 24 | profile `set-resources` silently ignored | pin on the CLI |
| TIMEOUT at 24 h | runtime written as `"24:00:00"`, not integer minutes | pin in minutes |
| iPHoP OOM | memory not applied from profile | 150 GB on the CLI |
| `cluster_phages` died | **our bug** — MEGAHIT headers contain spaces, breaking `sort`/`join` | strip at assembly |
| PHOLD 0 rows from 834 splits | conda env missing `scipy` + `scikit-learn` | installed |
| `taxonomic_consensus` type error | vConTACT3 name filter matched no MEGAHIT contigs → empty frame typed `<logical>` | widened filter |

**MMseqs2 `nr` taxonomy is skipped by design.** The shared index is 2,125 GB
across 9 split files on a `/ref` mount that reads at 7–10 MB/s — one pass would
take ~81 days. Three runs (27 h, 14 h, 10 h) were lost before measuring that.
`skip_mmseqs: true` routes the rule to its empty-output branch. Taxonomy comes
from PhaGCN and vConTACT3 instead; hosts from CHERRY and iPHoP.

## Caveats on the result

- ~60% of vOTUs have no iPHoP host call, so the *Lactobacillus* pattern is drawn
  from the predictable subset.
- 71 of 74 markers were called `virulent` by BACPHLIP — lopsided enough to
  suspect it behaves poorly on short contigs.
- Synthetic data: effect sizes are not transferable to real cervicovaginal
  samples.

## Related

The Hecatomb pipeline itself is documented at
<https://leran10.github.io/virome-pipeline-guide/>.
