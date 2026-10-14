# Jackalope cohort — phage community analysis by HPV status

**Question:** does the bacteriophage community differ between HPV-positive and
HPV-negative samples? HPV status is the *grouping variable*, not the readout —
we are not doing HPV detection.

## The data

100 simulated vaginal virome samples generated with
[ViroForge](https://github.com/Leran10/viroforge) (SLURM job 47308113,
14–15 Sep 2026), living on HTCF:

```
/scratch/sahlab/leran/jackalope_cohort/
├── VIR_001 … VIR_100/       # 43 GB: fastq/, fasta/, metadata/
├── READS/                   # staged: VIR_###_R{1,2}.fastq.gz  <- pipeline input
├── jackalope_clinical_metadata.csv
└── batch_report.json        # 100/100 successful
```

| | HPV− | HPV+ |
|---|---|---|
| n | 52 | 48 |
| ViroForge collection | 16 | 24 |
| Coverage | flat 10× | **5–20×**, scales with `Antler_Size_mm` |
| Read pairs | 283k – 1.36M per sample | |

`jackalope_clinical_metadata.csv` is a *blinded* clinical sheet: the case/control
label is `Antler_Present` (Yes = HPV+), with `Antler_Size_mm` as a graded severity
measure. It is the input that generated the cohort — `jackalope_cohort.py` reads
it and writes the ViroForge batch YAML — so antler size is physically encoded in
sequencing depth.

## Expected signal (from ground truth)

Every sample has a `metadata/*_composition.tsv` giving the true genome content, so
the answer is known in advance and the pipeline can be benchmarked against it:

| | HPV− | HPV+ |
|---|---|---|
| Phage genomes present | 7–10 | **18–23** |
| Total phage abundance | 0.04–0.21 | 0.11–0.22 |
| *Lactobacillus* phage | ~0.0024 | **~0.0050** |

**Richness is the strong signal** (~2–3× more phage genomes in HPV+, no overlap);
total abundance overlaps and will not discriminate.

## ⚠ The depth confound

HPV+ samples were simulated at 5–20× and HPV− at a flat 10×. Richness is highly
depth-sensitive, and richness is the leading signal — so the two are entangled by
construction. Every comparison must **rarefy to even depth or model depth as a
covariate**. Worth checking whether richness tracks coverage *within* the HPV+
group; if it does, that is the confound showing itself.

## Layout

```
config/project.yaml              paths, groups, confounds, pipeline settings
scripts/stage_reads.sh           gzip + rename FASTQs to VIR_###_R{1,2}.fastq.gz
scripts/collect_qc_attrition.py  per-step read attrition table from a Hecatomb run
pipelines/prophage-pipeline/     reads -> assembly -> bins -> prophages
pipelines/phage-analysis/        assembly -> vOTUs -> taxonomy/lifestyle/host
results/qc/                      QC outputs
```

The source FASTQs share filenames across samples (every HPV− sample is
`vaginal_virome_healthy_R1.fastq`), which breaks any pipeline that infers sample
identity from the filename. `stage_reads.sh` gzips them into `READS/` under their
`VIR_###` id, which also joins cleanly to the metadata.

## Plan

1. **Setup** — metadata, staging, pipeline clones, database checks ✅
2. **QC** — Hecatomb v1.1.0, `preprocessing` stage only, `--preprocess paired`
3. **Phage analysis** — `prophage-pipeline` from reads, then `phage-analysis`
   on the assemblies
4. **Stats** — richness and composition by group, controlling for depth,
   validated against ground truth

## Tools

QC uses the Hecatomb v1.1.0 install at
`/ref/sahlab/software/anaconda3/envs/hecatomb_v1.1.0` (73 GB databases, 14 masked
host genomes). Use `--preprocess paired`, **not** `roundAB` — these are simulated
Illumina reads without Round A/B primer artefacts, and the roundAB path in that
shared install has hand-edited NEBNext primer files.

Hecatomb's human reference is viral-masked and the mask includes 21,029
papillomavirus fragments, so HPV reads survive host removal.

All 14 databases needed by both phage pipelines are already present under
`/ref/sahlab/data/` — no downloads required. Both pipeline configs ship with HTCF
paths as their defaults.

## Related

The Hecatomb pipeline itself is documented at
<https://leran10.github.io/virome-pipeline-guide/>.
