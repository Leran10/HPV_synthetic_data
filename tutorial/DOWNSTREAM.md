# Assembly, annotation and analysis — why there are no scripts here

[`01_qc.sh`](01_qc.sh) and [`02_qc_cohort.sh`](02_qc_cohort.sh) are runnable on a
laptop. The stages after QC are not, and wrapping them in a standalone script
would promise something it cannot deliver. This page says what they need, so you
can decide where to run them — and documents the analysis half, which *is*
runnable.

## Why QC stops at the laptop boundary

QC needs one human genome (~3 GB) and tools that stream reads. Everything after
assembly needs reference databases measured in hundreds of gigabytes. Sizes
below are `du -sh` from the installation that produced our results, not vendor
estimates:

| Stage | Tool | Database | Size |
|---|---|---|---|
| Host prediction | iPHoP | `Aug_2023_pub_rw` | **318 GB** |
| Protein search | MMseqs2 | NCBI `nr` (multi-file index) | **2.5 TB** |
| Structure annotation | PHOLD | ProstT5 + Foldseek | **50 GB** |
| Completeness | CheckV | `checkv-db-v1.5` | 6.4 GB |
| Network taxonomy | vConTACT3 | `DB223` | 6.3 GB |
| Taxonomy + host | PHABOX 2 | `phabox_db_v2` | 1.7 GB |
| Viral identification | geNomad | `genomad_db` | 1.4 GB |
| | | **Total** | **~2.9 TB** (~385 GB without `nr`) |

Add the compute: PHOLD runs a protein language model over every predicted
protein; iPHoP runs four host-prediction methods and combines them; vConTACT3
builds a gene-sharing network across the whole set. Our cohort needed a GPU-less
cluster run of roughly a day per expensive stage, with jobs pinned at 24+ CPUs.

So for assembly onward, use the pipeline that already orchestrates this:
**[phage-analysis](https://github.com/megjohnson1999/phage-analysis)** by
megjohnson1999. Our exact invocation, with every resource override, is in
[`scripts/run_phage_analysis_vlp.sh`](../scripts/run_phage_analysis_vlp.sh), and
each tool's command and mechanism is documented in
[`docs/TOOLS_WALKTHROUGH.md`](../docs/TOOLS_WALKTHROUGH.md) steps 6–17.

One stage we could not run at all: **MMseqs2 taxonomy against `nr`**. The index
is 2.5 TB, and our `/ref` filesystem read it at 7–10 MB/s — a projected ~81 days
for a single pass. We skipped it deliberately and used PhaGCN for taxonomy
instead. If you hit the same wall, that is an infrastructure problem, not a
configuration one: measure your read throughput against the index size before
committing to the stage.

Excluding `nr`, the remaining databases total ~385 GB, which is within reach of
a well-provisioned workstation. `nr` is the one that forces a cluster.

## What you can assemble on a workstation

MEGAHIT is the one post-QC step with no database at all:

```bash
megahit -1 clean_R1.fastq.gz -2 clean_R2.fastq.gz \
        -o assembly_SAMPLE --k-min 21 --k-max 141 --k-step 20 \
        -t 16 -m 0.9
```

Per-sample assembly of a VLP library is minutes, not hours. Prefix the contig
names with the sample before pooling, or downstream joins will silently
mismatch:

```bash
sed -e "s/^>/>${SAMPLE}__/" -e "s/^\(>[^ ]*\).*/\1/" \
    assembly_SAMPLE/final.contigs.fa > SAMPLE.contigs.fa
```

That second expression strips everything after the first space. MEGAHIT writes
`>k141_213 flag=1 multi=2.0 len=5000`, and the space breaks any tool that
assumes the sequence ID is the whole header — it cost us a silent failure in
clustering. See the pitfalls in the walkthrough.

## The analysis half — this you can run

Once you have a vOTU abundance matrix, the statistics and figures need nothing
but Python. Both scripts are parameterised and run on the tracked result tables
in this repo:

```bash
pip install numpy scipy        # that is the whole dependency list

# 1. differential abundance: three passes per vOTU, FDR-corrected
python scripts/differential_abundance.py \
    --indir results/vlp \
    --metadata jackalope_clinical_metadata.csv \
    --outdir results/markers

# 2. regenerate every figure in the report
python scripts/make_figures.py --outdir docs/figures
```

### What `differential_abundance.py` does

Three tests per vOTU, because any one of them alone would mislead on this
cohort:

| Pass | Test | Guards against |
|---|---|---|
| raw | Mann-Whitney U | nothing — reported for reference only |
| rarefied | Mann-Whitney U after scaling every sample to the minimum depth | depth differences between groups |
| depth-adjusted | partial Spearman, log(coverage) partialled out | residual depth effect rarefaction misses |

A marker is `robust` only if it clears FDR 0.05 in **both** the rarefied and the
depth-adjusted pass. Benjamini-Hochberg throughout.

Then a second filter that matters more than the statistics. Each vOTU is labelled
by presence pattern — `HPV+_exclusive`, `HPV-_exclusive`, `shared_both`, or
`shared_rare_in_one`. On our cohort:

```
596 vOTUs tested
  421 robust but present in only ONE arm   <- collection artefact
   80 robust, rare in one arm              <- borderline
   74 robust and present in BOTH arms      <- the informative set
   21 not significant
```

The 421 are not findings. The two ViroForge source collections do not share
those genomes, so their absence from one arm is a property of the simulation,
not of HPV status. **A test that only asked "is this significant?" would have
reported 575 markers and been wrong about 421 of them.** This is the single most
transferable idea in the analysis: significance tells you a difference is real,
not that it is biological.

### What `make_figures.py` does

Regenerates all five report figures as standalone SVG, from
`results/markers/markers_case_control.tsv` and `results/vlp/votu_richness.tsv`:

| File | Shows |
|---|---|
| `fig1_composition.svg` | the 421 / 80 / 74 / 21 split above |
| `fig2_volcano.svg` | effect size vs depth-adjusted significance |
| `fig3_lactobacillus.svg` | 29 of 29 informative *Lactobacillus* vOTUs depleted in HPV+ |
| `fig4_host_genus.svg` | effect size by predicted host genus |
| `fig5_richness.svg` | richness vs depth, showing the confound directly |

No matplotlib — SVG is emitted directly, so the output matches the report's
style and the only dependency is the standard library. Each mark carries a
`<title>`, so hovering in a browser gives the vOTU name and numbers.

Useful flags:

```bash
python scripts/make_figures.py --host-genus Atopobium    # retarget figure 3
python scripts/make_figures.py --markers mine.tsv --richness mine.tsv --outdir figs/
```

Pass `--markers` your own table and the figures rebuild against it, provided it
carries `verdict`, `presence_pattern`, `log2fc`, `cliffs_delta`, `fdr_depthadj`
and `host_Host`.

> **Use `markers_case_control.tsv`, not `markers_consensus.tsv`.** The consensus
> file is already filtered to 402 rows and drops the `verdict` and
> `presence_pattern` columns the classification needs. Pointing the figure script
> at it silently reports every marker as informative — which it did, on my first
> run, until the output disagreed with `analysis_summary.txt`.

## Summary

| Stage | Runnable here? | Why |
|---|---|---|
| QC | **yes** — [`01_qc.sh`](01_qc.sh), [`02_qc_cohort.sh`](02_qc_cohort.sh) | one host genome, streaming tools |
| Assembly | yes, with the snippet above | MEGAHIT needs no database |
| Identification & annotation | **no** | ~385 GB of databases (plus 2.5 TB if you add `nr`), cluster-scale compute |
| Abundance | partly | bbmap is light, but needs the vOTU set from above |
| Statistics & figures | **yes** — `differential_abundance.py`, `make_figures.py` | numpy + scipy, or stdlib alone |
