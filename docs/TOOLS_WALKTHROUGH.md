# Every tool in this analysis, step by step

A walkthrough of the 18 steps — 16 distinct tools plus two scripted steps — that take us from simulated reads to a ranked list of
phage markers. For each: **what it does**, **the mechanism** (how it actually
decides), **the command we ran**, and **what it produced on this cohort**.

Commands are reproduced from source as executed, with site paths shown as
`${PLACEHOLDERS}` (see the README table). Provenance, so you know what to trust:

- **Tools 1, 2, 6, 17, 18** — our own scripts in [`scripts/`](../scripts), run directly
- **Tools 3, 4, 5** — Hecatomb v1.1.0's internal rules (`01_preprocessing.smk`)
- **Tools 7–16** — the [phage-analysis](https://github.com/megjohnson1999/phage-analysis)
  pipeline's rules (`01_prediction.smk`, `02_clustering.smk`, `03_analysis.smk`)

Every command below was checked against the rule file that ran it. Snakemake
wildcards are resolved to readable paths; flags and values are verbatim.

**Running order:**

```
 generate ──▶ stage ──▶ QC ──▶ assemble ──▶ identify ──▶ cluster ──▶ annotate ──▶ quantify ──▶ test
 ViroForge   gzip/     fastp   MEGAHIT     Jaeger       vclust      PHOLD        bbmap        Python
             rename    minimap2            geNomad                  PHABOX       (coverage)   (stats)
                       MMseqs2             CheckV                   vConTACT3
                                                                    BACPHLIP
                                                                    iPHoP
```

---

## Stage 0 — Building the dataset

### 1. ViroForge 0.20.0 — synthetic virome generator

**What it does.** Builds FASTQ files with known ground truth: you choose a viral
collection, a sequencing platform, a coverage depth, and a contamination profile,
and it writes reads plus a manifest of exactly which genomes went in at what
abundance.

**Mechanism.** Picks genomes from a curated collection, assigns relative
abundances, optionally applies a VLP-enrichment model (size filtration + nuclease),
adds host/bacterial/fungal/rRNA contamination, then calls **InSilicoSeq** to
simulate platform-specific reads with realistic error profiles. Abundance is
allocated as a *share of reads*, not scaled by genome length — which matters, see
the pitfall below.

```bash
viroforge generate \
    --collection-id 24 \          # 24 = HPV-infected vaginal; 16 = healthy
    --platform novaseq \
    --coverage 12.4 \             # per-sample, from Antler_Size_mm
    --seed 1041 \
    --output ${PROJECT_DIR}/data/vlp/VIR_041
```

**Produced.** 100 samples. VLP cohort ≈98% viral; a second bulk cohort
(`--no-vlp --contamination-level heavy`) at ~40% host DNA.

> **Pitfall that cost us the first analysis.** In bulk mode viruses are ~1% of
> reads, split across 25–54 genomes per sample. Because abundance is a read share
> rather than per-base coverage, a 50 kb phage received ~1× and could not
> assemble, while a 7 kb HPV genome got ~25×. Only 211 viral contigs came back
> from 1,287 genomes actually present. The VLP cohort lifts the same phage to
> ~110×. **If you want phages to assemble, enrich.**

---

## Stage 1 — Read preparation

### 2. gzip + rename — sample identity

**Why it needs a step.** ViroForge names output by collection, so every HPV−
sample was `vaginal_virome_healthy_R1.fastq`. Identical filenames break any
pipeline that infers the sample from the path.

```bash
gzip -c "$src" > ${PROJECT_DIR}/data/vlp/READS/VIR_041_R1.fastq.gz
```

**Produced.** 200 files as `VIR_###_R{1,2}.fastq.gz`, joinable to the metadata.

### 3. fastp 0.23.2 — read QC

**What it does.** One pass that trims adapters, quality-trims, length-filters,
removes poly-X tails and deduplicates.

**Mechanism.** Adapters are detected by overlap analysis of read pairs rather than
a reference list. Quality trimming uses a sliding window from the 3′ end.
Deduplication hashes read pairs and drops repeats.

```bash
fastp -i VIR_041_R1.fastq.gz -I VIR_041_R2.fastq.gz \
      -o trimmed_R1.fastq -O trimmed_R2.fastq \
      --qualified_quality_phred 15 --length_required 90 \
      --detect_adapter_for_pe \
      --cut_tail --cut_tail_window_size 25 --cut_tail_mean_quality 15 \
      --dedup --dup_calc_accuracy 4 --trim_poly_x --thread 12
```

**Produced (bulk cohort, 100 samples).** 134,958,156 reads in → 108,329,338 out.
Deduplication removed 23.5M reads (17.4%); quality and length filters together
under 2.3%.

### 4. minimap2 ≥2.20 + samtools ≥1.13 — host removal

**What it does.** Discards reads of human origin.

**Mechanism.** minimap2 indexes the host genome by minimizers (a sampled subset of
k-mers) and seeds alignments from them, so it maps short reads against a 3 Gb
genome quickly. `samtools view -f 4` then keeps only records flagged *unmapped*.
The reference is **viral-masked** — virus-like regions of the human genome are
blanked, so genuine viral reads are not thrown away with the host.

```bash
minimap2 -ax sr -t 16 --secondary=no ${HOST_INDEX} R1.fastq R2.fastq \
  | samtools view -f 4 -h \
  | samtools fastq -NO -1 unmapped_R1.fastq -2 unmapped_R2.fastq
```

**Produced.** Bulk cohort: 43.3% of reads removed. VLP cohort: 1.4% — the VLP
model had already stripped host DNA to 0.67%, so there was little left to remove.

### 5. MMseqs2 12.113e3 — dereplication (inside Hecatomb QC)

**What it does.** Collapses near-identical reads to representatives, carrying
counts forward.

**Mechanism.** `linclust` is linear-time clustering: each sequence contributes a
few k-mers, sequences sharing k-mers are compared, and a representative is chosen
per cluster. Crucially **lossless for abundance** — every read's count transfers
to its representative.

```bash
mmseqs easy-linclust reads.fastq out_prefix tmp_dir \
    --kmer-per-seq-scale 0.3 -c 0.8 --cov-mode 1 \
    --min-seq-id 0.97 --alignment-mode 3 --threads 8
```

**Produced.** 30,733,965 read pairs → 5,108,133 representatives, a 6× smaller
search space with zero reads lost.

> **Deduplication vs dereplication — different operations.** Deduplication
> (fastp) *deletes* PCR duplicates as artefacts. Dereplication (MMseqs2) *groups*
> distinct molecules and keeps their counts. One loses reads on purpose; the other
> must not.

---

## Stage 2 — Assembly

### 6. MEGAHIT 1.2.9 — metagenome assembler

**What it does.** Reconstructs contigs from short reads.

**Mechanism.** Builds a **succinct de Bruijn graph** — a memory-compressed index
of k-mer overlaps — and iterates over increasing k. Small k resolves low-coverage
regions; large k resolves repeats. Assembly time scales with graph size, not read
count, which is why a small high-coverage viral assembly finishes in under a
minute.

```bash
megahit -1 VIR_041_R1.fastq.gz -2 VIR_041_R2.fastq.gz \
    -o ${PROJECT_DIR}/data/vlp/assemblies/VIR_041 --out-prefix VIR_041 \
    --k-min 21 --k-max 141 --k-step 12 \
    --min-contig-len 1000 -t 16

# sample-prefix the headers AND strip the description after the first space
sed -e "s/^>/>VIR_041__/" -e "s/^\(>[^ ]*\).*/\1/" VIR_041.contigs.fa > out.fasta
```

**Produced.** 17,046 contigs pooled across 100 samples, N50 ≈ 30 kb, longest
193 kb (82% of a CMV genome). ~47 s per sample.

> **Pitfall.** MEGAHIT writes `>k141_30 flag=1 multi=7.0 len=1244`. Those embedded
> spaces silently broke vOTU clustering downstream, because `seqkit fx2tab --name`
> keeps the whole description as field 1 and the later `sort`/`join` split on
> whitespace. Strip everything after the first space at assembly time.

---

## Stage 3 — Which contigs are viral?

Four independent callers, combined by consensus. No single tool is trusted alone.

### 7. Jaeger 1.1.26 — deep-learning viral classifier

**What it does.** Scores each contig on how virus-like its sequence composition is.

**Mechanism.** A convolutional neural network over one-hot encoded DNA, trained on
known viral and non-viral genomes. It learns compositional signatures — k-mer
usage, codon bias — rather than matching a database, so it can flag novel viruses
with no close relative in any reference.

```bash
Jaeger -i pooled_contigs.fasta -o ${PROJECT_DIR}/.../01_jaeger_output \
    -s 2.5 --fsize 1000 --stride 1000
```

**Produced.** All 17,046 contigs scored.

### 8. geNomad 1.8.0 — hybrid gene-based classifier

**What it does.** Identifies viral and plasmid sequences, and marks their topology.

**Mechanism.** Two branches combined. A **marker branch** calls genes and matches
them to curated viral/chromosomal/plasmid profiles. A **neural branch** classifies
raw sequence. It also detects terminal repeats (DTR/ITR) and provirus boundaries —
structural evidence of a viral genome, independent of composition.

```bash
genomad end-to-end --min-score 0.6 --cleanup --threads 16 \
    pooled_contigs.fasta ${PROJECT_DIR}/.../01_genomad_output ${DB_DIR}/genomad_db
```

**Produced.** 5,997 contigs called viral.

### 9. CheckV 1.0+ — completeness and contamination

**What it does.** Estimates how complete each viral genome is and trims host
sequence left on the ends.

**Mechanism.** Compares each contig against a database of complete viral genomes
to estimate expected length (completeness = observed / expected). Separately scans
for a host-to-viral gene transition along the contig, which marks where an
integrated provirus ends.

```bash
checkv end_to_end pooled_contigs.fasta ${PROJECT_DIR}/.../01_checkv_output \
    -d ${DB_DIR}/checkv-db-v1.5 -t 16
```

**Produced.** Quality for all 17,046 contigs, and again on the 2,921 vOTUs:
6 high-quality, 7 medium, 195 low (in the earlier bulk run).

### 10. PHOLD 0.2.0 — structure-aware protein annotation

**What it does.** Annotates phage genes by predicted protein *structure*, and
supplies the functional diversity score the consensus depends on.

**Mechanism.** Calls genes with Pyrodigal, predicts each protein's 3D structure
with ProstT5 (a protein language model), then searches structures with
**Foldseek** against a phage protein database. Structure is conserved long after
sequence has diverged, so this annotates genes that sequence search misses.

```bash
phold run -i contig_chunk.fasta -o phold_out --cpu -t 8
```

**Produced.** 127,352 per-CDS annotations.

> **Pitfall that silenced the whole analysis.** PHOLD's conda env was missing
> `scipy` and `scikit-learn`, so `transformers` failed to import and PHOLD crashed
> on all 834 chunks — each time writing an *empty* output file rather than
> erroring. Because the consensus rule derives `Functionaldiversity` from PHOLD,
> every contig scored 0 and four of the five routes to calling a phage became
> unreachable. **Phage calls collapsed from 4,601 to 214.** Always check that a
> tool produced rows, not just that the job exited 0.

### 11. Consensus integration (R) — combining the four callers

**Mechanism.** A contig is called a phage on any of:

| rule | condition |
|---|---|
| high functional diversity | `Functionaldiversity >= 3` |
| special topology | DTR / ITR / Provirus from geNomad |
| moderate diversity + one tool | diversity 1–3 **and** (Jaeger ≥0.6 **or** geNomad score) |
| dual tool agreement | Jaeger ≥0.6 **and** geNomad score |

**Produced.** 4,601 phage contigs of 17,046 viral.

---

## Stage 4 — Collapsing to species-level units

### 12. vclust 1.2.9 — vOTU clustering

**What it does.** Groups contigs into viral Operational Taxonomic Units — roughly
"species" — so the same phage found in 40 samples counts once.

**Mechanism.** Three steps. `prefilter` finds candidate pairs by k-mer sketching
(cheap, avoids all-vs-all). `align` computes average nucleotide identity for
survivors. `cluster` applies the community standard — **95% ANI over 85%
coverage** — using Leiden community detection on the similarity graph.

```bash
vclust prefilter -i phageContigs.fasta -o vclust_fltr.txt
vclust align    -i phageContigs.fasta -o vclust_ani.tsv --filter vclust_fltr.txt
vclust cluster  -i vclust_ani.tsv -o vclust_clusters.tsv \
                --algorithm leiden --metric ani --ani 0.95 --qcov 0.85 --rcov 0.85
```

**Produced.** 4,601 contigs → **2,921 vOTUs**.

---

## Stage 5 — What is each vOTU, and what does it infect?

### 13. PHABOX 2 — taxonomy and host, four sub-tools

**Mechanism.** Four graph/network models over a protein-sharing network:
**PhaMer** (is it a phage), **PhaGCN** (taxonomy by graph convolution over shared
proteins), **CHERRY** (host prediction), **PhaTYP** (lifestyle).

```bash
phabox2 --task end_to_end --dbdir ${DB_DIR}/phabox_db_v2 \
    --contigs vOTU_repSeqs.fasta --outpth phabox_output --threads 16
```

**Produced.** 2,921 taxonomy calls (PhaGCN) and 2,921 host calls (CHERRY).

### 14. vConTACT3 — network taxonomy

**What it does.** A second, independent taxonomy using gene sharing.

**Mechanism.** Calls proteins, builds a bipartite gene-sharing network against
reference viruses, and clusters it. vOTUs landing in a cluster with known viruses
inherit that taxonomy; isolated clusters are novel genera.

```bash
vcontact3 run --nucleotide vOTU_repSeqs.fasta --output vc3_output \
    --db-path ${DB_DIR}/vcontact3_DB223
```

**Produced.** 8,226 assignment rows covering all 2,921 vOTUs.

> **Pitfall.** Its contig filter only matched names containing `contig`, `edge`, or
> `virus_comp`. MEGAHIT names are `k141_*`, so all 8,226 rows were dropped,
> leaving an empty data frame whose ID column typed as `<logical>` and crashed the
> downstream join. Check that a filter matches *your* naming.

### 15. BACPHLIP 0.9.6 — lifestyle prediction

**What it does.** Classifies each phage as **temperate** (can integrate and lie
dormant) or **virulent** (always lyses).

**Mechanism.** HMMER searches for ~200 protein domains diagnostic of lysogeny —
integrases, excisionases, repressors — then a random-forest classifier on domain
presence/absence.

```bash
bacphlip -i vOTU_repSeqs.fasta --multi_fasta -f
```

**Produced.** 2,921 lifestyle calls. **Treat with caution:** 71 of our 74 markers
came back virulent, lopsided enough to suspect poor performance on short contigs
rather than real biology.

### 16. iPHoP 1.3.0 — host prediction, independent of CHERRY

**What it does.** Predicts the bacterial host, integrating several lines of
evidence.

**Mechanism.** Runs four methods and combines their scores into a calibrated
confidence: **CRISPR spacer matches** (a bacterium that recorded this phage),
**k-mer composition** similarity to host genomes, **BLAST** homology, and
**WIsH** (Markov models of host genomes).

```bash
iphop predict --fa_file vOTU_chunk.fasta --out_dir iphop_out \
    --db_dir ${DB_DIR}/iphop_DBs/Aug_2023_pub_rw --num_threads 16
```

**Produced.** 9,545 host predictions covering 1,175 of 2,921 vOTUs (40%).

> **Why only 40%?** Host prediction matches against *known* bacteria, so it
> succeeds for well-studied hosts and fails for obscure ones. The missingness is
> systematic, not random — and in our data the assigned and unassigned markers
> behaved differently (33/41 vs 16/33 depleted). Any claim must be scoped to the
> assignable subset.

---

## Stage 6 — Quantification

### 17. bbmap 38.90 — per-sample abundance

**What it does.** Measures how much of each vOTU is present in each sample.

**Mechanism.** Maps every sample's reads back to the vOTU representatives and
reports per-reference depth *and* **breadth** (percent of the reference covered).
Breadth matters as much as depth: high depth over 5% of a genome is a
conserved-region artefact, not a genome that is there.

```bash
bbmap.sh ref=vOTU_repSeqs.fasta \
    in=VIR_041_R1.fastq.gz in2=VIR_041_R2.fastq.gz \
    covstats=VIR_041.cov_stats rpkm=VIR_041.rpkm \
    minid=90 physcov=t ambiguous=random threads=8 -Xmx20g
```

**Produced.** A 2,921 × 100 coverage matrix. ~1 min per sample.

> **This step was missing from our first attempt.** The pipeline's own CoverM rule
> never ran, and the cluster table listed only representatives — so we had 2,921
> vOTUs and *no per-sample counts*, which means no comparison was possible at all.
> A catalogue is not an analysis.

---

## Stage 7 — Statistics

### 18. Python (scipy/numpy) — differential abundance

**What it does.** Tests each vOTU for a difference between HPV+ and HPV−, and for
a trend against severity.

**Mechanism.** Three passes per vOTU, because sequencing depth is confounded with
group in this cohort:

| pass | test | purpose |
|---|---|---|
| raw | Mann-Whitney U | the naive comparison |
| rarefied | Mann-Whitney on depth-scaled values | removes depth differences |
| depth-adjusted | partial Spearman, log(coverage) partialled out | models depth instead |

A marker counts as `robust` only if it clears **both** the rarefied and
depth-adjusted tests. Benjamini-Hochberg FDR across all tested vOTUs. Effect size
is Cliff's delta (non-parametric, −1 to +1).

```bash
python3 scripts/differential_abundance.py \
    --indir results/vlp \
    --metadata jackalope_clinical_metadata.csv \
    --outdir results/markers
```

**Produced.** 596 testable vOTUs, 575 robust markers — of which **421 are a
design artefact** (the two arms came from disjoint genome collections, so those
separate perfectly by construction) and **74 are informative**.

---

## What the whole chain yielded

| stage | output |
|---|---|
| reads | 100 samples, ~98% viral (VLP) |
| assembly | 17,046 contigs, N50 ≈ 30 kb |
| viral identification | 4,601 phage contigs |
| clustering | 2,921 vOTUs |
| annotation | taxonomy, host, lifestyle, completeness for all 2,921 |
| abundance | 2,921 × 100 matrix |
| statistics | 74 informative markers |

Results and figures: [results_summary.html](results_summary.html).

---

## The five lessons worth carrying

1. **Enrichment decides what you can analyse.** ~1% viral content means large
   phages cannot assemble, no matter how much you sequence.
2. **A tool exiting 0 is not success.** PHOLD "succeeded" 834 times while writing
   empty files, and quietly cut phage calls by 95%.
3. **Check that filters match your naming.** Two separate failures traced to
   MEGAHIT's header format — spaces breaking a join, and a hardcoded name pattern
   matching nothing.
4. **A catalogue is not an analysis.** vOTUs without per-sample counts cannot
   answer a comparative question.
5. **Control confounds explicitly and report what survives.** Labelling markers
   `robust` / `depth_suspect` / artefact is more honest than one p-value — and it
   was what separated 74 real candidates from 421 artefacts.
