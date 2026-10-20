# Runnable QC scripts

The scripts in [`../scripts/`](../scripts) are the ones we actually ran: SLURM
batch jobs for one specific cluster, driving Hecatomb through a site Snakemake
profile. They are here as a provenance record, not as something you can run.

These two scripts are the tutorial version — **no cluster, no Snakemake, no
Hecatomb.** Just the tools, with the same parameters, on your own data.

| Script | Scope |
|---|---|
| `01_qc.sh` | one sample, end to end, prints a per-step read-count table |
| `02_qc_cohort.sh` | runs `01_qc.sh` over a directory, builds the attrition tables |

They implement the operations in the [QC tool
list](../docs/TOOLS_WALKTHROUGH.md#qc-tool-list--every-operation-and-the-reads-it-removed).

## Install

```bash
conda create -n virome-qc -c bioconda -c conda-forge \
    fastp minimap2 samtools bbmap mmseqs2 jq
conda activate virome-qc
```

Minimum is `fastp`. `minimap2` + `samtools` add host removal; `bbmap` provides
`bbduk.sh` (PhiX) and `repair.sh` (re-pairing); `mmseqs2` is only needed for
`--derep`; `jq` enriches the per-reason breakdown. The scripts check what is on
your PATH and tell you exactly what is missing before doing any work.

If `repair.sh` is unavailable, `01_qc.sh` falls back to a samtools-only
re-pairing route automatically.

## One sample

```bash
./01_qc.sh -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
           -s SAMPLE01 -o qc_out \
           -r GRCh38.fa.gz \
           -p phix174_ill.ref.fa.gz \
           -t 8
```

| Flag | Meaning |
|---|---|
| `-1` / `-2` | input read pairs (required) |
| `-s` | sample name, used for output filenames (required) |
| `-o` | output directory (required) |
| `-r` | host genome to remove. **Omit to skip host removal** |
| `-p` | PhiX/vector reference. **Omit to skip** — only safe on simulated data |
| `-t` | threads (default 8) |
| `--derep` | also dereplicate with MMseqs2 (skip unless a big DB search follows) |

Outputs `SAMPLE01_R{1,2}.clean.fastq.gz`, the fastp HTML/JSON report, and
`SAMPLE01.readcounts.tsv`:

```
step   operation                     tool               reads_in  reads_out  reads_removed
1      trim+filter+dedup             fastp              680       502        178
1.x    discarded: low_quality_reads  fastp              -         -          0
1.x    discarded: too_short_reads    fastp              -         -          120
1b     remove PhiX+vector            bbduk              -         -          SKIPPED
2      remove host                   minimap2+samtools  502       422        80
3      re-pair orphans               samtools           422       422        0
total  retained                      -                  680       422        258
```

## A whole cohort

Expects `READS_DIR/<sample>_R1.fastq.gz` and `_R2.fastq.gz`:

```bash
./02_qc_cohort.sh -i reads/ -o qc_out/ \
                  -r GRCh38.fa.gz -t 4 -j 6
```

`-j 6` runs six samples at once, each on `-t 4` threads — keep `j × t` at or
under your core count. Re-running skips samples that already have a
`.readcounts.tsv`, so an interrupted run resumes.

Produces, besides the per-sample files:

- **`qc_attrition.tsv`** — one row per sample. The tutorial equivalent of
  [`results/qc/qc_attrition.tsv`](../results/qc/qc_attrition.tsv).
- **`qc_steps_total.tsv`** — cohort totals, equivalent to
  [`results/qc/qc_steps_vlp.tsv`](../results/qc/qc_steps_vlp.tsv).

## Three things worth knowing before you run this

**1. Set `-p` on real data.** We did not, because our reads are
InSilicoSeq-simulated and contain no PhiX. Illumina runs are routinely
PhiX-spiked and not every demultiplexer strips it. Without `-p`, the step is
skipped and the table says `SKIPPED` rather than silently reporting zero.

**2. Re-pairing is not optional, and it is the step most hand-written QC
scripts omit.** `samtools view -f 4` filters reads independently, so a pair with
one human mate leaves an orphan; R1 and R2 end up different lengths and out of
order, and MEGAHIT or any aligner will either error or silently misread them.
Step 3 exists only to fix this.

**3. Read your attrition table as a result, not a receipt.** On our VLP cohort
host removal took 0.86% of reads — that is enrichment working as intended. The
bulk cohort, same tools, lost 34.72% to host. If your VLP library comes back 40%
human, the bench prep failed and no software choice rescues it. Equally, if
duplicates are most of your loss (they were 25.4% of our raw reads, ~93% of all
fastp removal), that is a library-complexity signal, not a QC problem to tune
away.

## Verifying your own run

A quick check that host removal and pairing both did their job:

```bash
# should report 0
minimap2 -ax sr -t 4 GRCh38.fa.gz qc_out/SAMPLE01_R1.clean.fastq.gz 2>/dev/null \
  | samtools view -c -F 4 -

# should be equal
gzip -cd qc_out/SAMPLE01_R1.clean.fastq.gz | wc -l
gzip -cd qc_out/SAMPLE01_R2.clean.fastq.gz | wc -l
```

## Scope

These cover QC only — steps 1–5 of the walkthrough. Assembly onward (MEGAHIT,
Jaeger, geNomad, CheckV, PHOLD, vclust, PHABOX, vConTACT3, BACPHLIP, iPHoP,
bbmap, the statistics) is documented with commands in
[`../docs/TOOLS_WALKTHROUGH.md`](../docs/TOOLS_WALKTHROUGH.md) but is not
wrapped into standalone scripts: those stages need large reference databases
(iPHoP ~900 GB, geNomad ~1.6 GB, CheckV, PHABOX) and realistically a cluster.
Use the [phage-analysis](https://github.com/megjohnson1999/phage-analysis)
pipeline for that half.

## Tested

Both scripts were run end to end on a synthetic fixture with planted defects —
adapter read-through, exact duplicates, low-quality tails, short reads,
ambiguous bases, and reads drawn from a known host genome. Host removal
recovered exactly the planted host reads, zero host reads survived re-screening,
and R1/R2 stayed in sync. Verified on bash 3.2 (macOS default) and in parallel
mode.
