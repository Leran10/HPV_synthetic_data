# Where everything lives

Three locations, each with a different job:

| | holds | backed up? |
|---|---|---|
| **GitHub** `Leran10/HPV_synthetic_data` (private) | code, configs, small result tables | yes, this is the backup |
| **Local** `~/Documents/github_projects/HPV_synthetic_data` | working copy of the above | yes, via GitHub |
| **HTCF** `${PROJECT_DIR}/` | ~230 GB of reads, assemblies, pipeline output | **no — scratch only** |

The large data is gitignored and exists in one place. It is reproducible from the
scripts, but the weekly purge-protection cron is what keeps it alive day to day
(see the bottom of this file).

---

## HTCF — the data

### VLP cohort — this is the one the analysis used

```
${PROJECT_DIR}/data/vlp/                       91G
├── VIR_001 … VIR_100/                 ViroForge output per sample
│   ├── fastq/                           raw R1/R2 (uncompressed)
│   ├── fasta/                           source genomes that went in
│   └── metadata/                        GROUND TRUTH: *_composition.tsv,
│                                        *_metadata.json (pre/post-VLP fractions)
├── READS/                             15G   staged VIR_###_R{1,2}.fastq.gz
├── assemblies/                        253M  VIR_###.contigs_1000bp.fasta
├── qc_hecatomb/                        27G  Hecatomb preprocessing
│   ├── results/seqtable.fasta                5.1M representative sequences
│   └── processing/stats/p01/*.json           per-sample fastp reports
├── phage_analysis/                    6.0G
│   ├── all_samples_contigs.fasta             17,046 pooled contigs
│   └── out/                                  see the pipeline section below
├── votu_abundance/                     35M  VIR_###.cov_stats (bbmap, 100 files)
├── results/                           3.7M  vOTU matrices (copied into the repo)
├── logs/                                    SLURM driver logs
└── jackalope_clinical_metadata.csv          the blinded sample sheet
```

### Bulk cohort — superseded for phage work, kept for the QC lesson

```
${PROJECT_DIR}/data/bulk/                        122G
├── VIR_001 … VIR_100/                 ~40% host DNA version
├── READS/                                   staged reads
├── qc_hecatomb/                        18G  QC where host removal does real work
├── prophage_out/                       45G  prophage-pipeline, died at binning;
│                                            its MEGAHIT assemblies are still here
├── phage_analysis/                    1.5G  earlier run, 210 vOTUs
├── votu_abundance/                          per-sample abundance
└── qc_*.tsv                                 QC attrition and step tables
```

### Original September batch — superseded

```
${PROJECT_DIR}/data/sept_batch/     19G
├── healthy/      healthy_S101 … S125
└── hpv_infected/ hpv_S201 … S225
```

### Backup of the small files

```
${PROJECT_DIR}/backup/                      98K
├── config/  docs/  scripts/  results/
```

---

## HTCF — pipeline output in detail

`data/vlp/phage_analysis/out/` — job 49642184, 435/435 steps, 3h51m:

| path | contents |
|---|---|
| `final_contig_summary.tsv` | 2,921 rows — the pipeline's own master table |
| `02_clustering/vOTU_repSeqs.fasta` | the 2,921 vOTU representative sequences |
| `02_clustering/clusters.tsv` | which contig belongs to which vOTU |
| `03_checkv_final/quality_summary.tsv` | CheckV completeness/quality per vOTU |
| `03_genomic_info/consensus_taxonomy.tsv` | 2,688 — merged PhaGCN + vConTACT3 |
| `03_genomic_info/lifestyle_consensus.tsv` | 2,921 — BACPHLIP temperate/virulent |
| `03_genomic_info/phabox_output/final_prediction/phagcn_prediction.tsv` | taxonomy |
| `03_genomic_info/phabox_output/final_prediction/cherry_prediction.tsv` | host (CHERRY) |
| `03_genomic_info/vc3_output/exports/final_assignments.csv` | vConTACT3 taxonomy |
| `03_iphop_results/iphop_predictions_compiled.tsv` | **9,545 iPHoP host calls** |
| `01_phold_output/phold_per_cds_predictions.tsv` | 127,352 PHOLD gene annotations |
| `01_phage_predictions/phageContigs.fasta` | 4,601 contigs called phage |
| `Pipeline_Summary_Report.html` | the pipeline's interactive report |
| `pipeline_summaries/*.json` | per-stage counts |

`mmseqs_taxonomy.tsv` and `mmseqs_output_*` are **intentionally empty** — see the
README on why the `nr` search is skipped.

---

## HTCF — tools and databases (read-only, shared)

| what | where |
|---|---|
| Hecatomb v1.1.0 + 73 GB databases | `${HECATOMB_ENV}` |
| Snakemake 8.20.3 + SLURM plugin | `${CONDA_BASE}/envs/snakemake_v8.20.3` |
| phage-analysis conda envs (24 GB, prebuilt) | `${CONDA_BASE}/envs/smk_envs_phage_pred` |
| prophage-pipeline conda envs (18 GB) | `${CONDA_BASE}/envs/smk_envs_prophage_pipeline` |
| ViroForge 0.20.0 env | `${VIROFORGE_ENV}` |
| ViroForge source | `${VIROFORGE_DIR}` |
| phage-analysis clone (patched) | `${PROJECT_DIR}/pipelines/phage-analysis` |
| prophage-pipeline clone | `${PROJECT_DIR}/pipelines/prophage-pipeline` |
| reference databases (geNomad, CheckV, iPHoP, PHABOX, vConTACT3, GTDB, nr) | `${DB_DIR}/` |

**Local patches applied to the phage-analysis clone** — these are not upstream, so
a fresh clone will fail the same way:

- `workflow/rules/03_analysis.smk` — `skip_mmseqs` flag; `touch` the third output
  in both empty-output branches. Backups: `.orig`, `.bak`
- `workflow/scripts/01_phagePrediction.R` — skip absent/empty mmseqs evidence
- `workflow/scripts/taxonomic_consensus.R` — vConTACT3 name filter widened to
  match `__`-prefixed contigs. Backup: `.bak2`
- `profile/slurm/config.yaml` — `cpus_per_task` in default-resources. Backup: `.orig`
- conda env `a78b593cb45d769cd0d022d77beef7b3_` — scipy + scikit-learn installed
  (this fixed PHOLD; it affects anyone else using that shared env)

---

## In this repo

```
config/
  project.yaml                 paths, groups, confounds
  hecatomb.qc.yaml             QC config, VLP
  hecatomb.qc.bulk.yaml        QC config, bulk
docs/
  FILE_LOCATIONS.md            this file
  viroforge_fragment_bias_issue.md   issue report for the ViroForge authors
scripts/                       16 scripts, generation through analysis
results/
  qc/                          attrition + step tables, both cohorts
  vlp/votu_abundance_matrix.tsv      2,921 × 100 coverage
  vlp/votu_presence_matrix.tsv       2,921 × 100 presence/absence
  vlp/votu_annotations.tsv           taxonomy + host + lifestyle + quality
  vlp/votu_richness.tsv              per-sample richness, group, depth
  markers/markers_case_control.tsv   596 testable vOTUs, 3 statistical passes
  markers/markers_severity.tsv       severity correlation, raw + depth-adjusted
  markers/markers_consensus.tsv      significant in both analyses
  markers/analysis_summary.txt       counts, thresholds, caveats
jackalope_clinical_metadata.csv      the blinded sample sheet (also on HTCF)
```

---

## Purge protection

HTCF deletes `/scratch` files whose mtime exceeds **28 days**. A weekly cron
(Sundays 03:00) touches every file in the project directories, so the clock never
passes 7 days — roughly 21 days of margin.

```
crontab -l                      # on HTCF
0 3 * * 0 find ${PROJECT_DIR} -type f -exec touch {} +
```

One entry now covers the whole project, replacing four separate ones from before
the directories were consolidated.

Backup of the pre-existing crontab: `~/crontab.bak.20261014`.

**Worth verifying occasionally** that the cron actually fires — check an mtime on
a `jackalope_cohort` file on a Monday. Silent cron failure is the realistic way
this data would be lost.

## Key SLURM jobs

| job | what | outcome |
|---|---|---|
| 47308113 | ViroForge, VLP cohort (Sept 14–15) | 100/100 |
| 48885443 | ViroForge, bulk cohort | 100/100 |
| 48885700 | Hecatomb QC, bulk | COMPLETED 8m |
| 49480692 | MEGAHIT assembly, VLP | 100/100 |
| 48989998 | prophage-pipeline | FAILED at 25% — no bins possible |
| **49642184** | **phage-analysis, VLP** | **COMPLETED 435/435, 3h51m** |
| **49642535** | **vOTU abundance, 100 samples** | **COMPLETED** |
