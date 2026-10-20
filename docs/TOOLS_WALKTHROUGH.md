# Every tool in this analysis, step by step

A walkthrough of the 20 steps — 16 distinct tools, two workflow wrappers, and two scripted steps — that take us from simulated reads to a ranked list of
phage markers. For each: **what it does**, **the mechanism** (how it actually
decides), **the command we ran**, and **what it produced on this cohort**.

Commands are reproduced from source as executed, with site paths shown as
`${PLACEHOLDERS}` (see the README table). Provenance, so you know what to trust:

- **Tools 1, 2, 6, 17, 18** — our own scripts in [`scripts/`](../scripts), run directly
- **Tools 3, 4, 5** — run BY Hecatomb v1.1.0 (commands from its `01_preprocessing.smk`);
  standalone equivalents given, since the wrapper is optional
- **Tools 7–16** — the [phage-analysis](https://github.com/megjohnson1999/phage-analysis)
  pipeline's rules (`01_prediction.smk`, `02_clustering.smk`, `03_analysis.smk`)

Every command below was checked against the rule file that ran it. Snakemake
wildcards are resolved to readable paths; flags and values are verbatim.

**Running order:**

```
 generate    stage      QC          assemble   identify    cluster   annotate    quantify  test
 ViroForge   gzip/      fastp       MEGAHIT    Jaeger      vclust    PHOLD       bbmap     Python
             rename     minimap2               geNomad               PHABOX      coverage  scipy
                        MMseqs2                CheckV                vConTACT3
                                                                     BACPHLIP
                                                                     iPHoP
 ─ours─      ─ours─     └── Hecatomb ──┘  ─ours─  └──── phage-analysis ─────┘   ─ours─   ─ours─
```

Two of those columns are **workflow wrappers**, not tools: Hecatomb ran the three
QC tools, phage-analysis ran the ten identification/annotation tools. Each is one
command, and each gets its own section below so you can see what we actually
typed — then the tool sections explain what happened inside.

**If you are here to learn how to do this on your own data,** the wrapper
sections also lay out your alternatives: [how QC was done and what your options
are](#how-qc-was-done-here--and-what-your-options-are) covers running the three
QC tools individually, choosing a host-removal tool, and which steps you can skip
outright. Nothing in this pipeline requires our specific software stack — but the
read counts we report do, so both routes are spelled out.

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

### How QC was done here — and what your options are

Three questions a reader needs answered before the tool pages below make sense:
what did *we* run, do *you* need the same thing, and where are the real choices.

#### What we ran

One command. All three QC tools below were orchestrated by **Hecatomb**
([shandley/hecatomb](https://github.com/shandley/hecatomb), Roach et al.) — a
published Snakemake virome workflow:

```bash
source ${CONDA_BASE}/bin/activate hecatomb_v1.1.0

hecatomb run \
    --reads   ${PROJECT_DIR}/data/vlp/READS \
    --output  ${PROJECT_DIR}/data/vlp/qc_hecatomb \
    --preprocess paired \        # one fastp pass; `roundAB` = six bbduk passes instead
    --host    human \            # selects Hecatomb's viral-masked GRCh38
    --search  fast \
    --threads 16 \
    --profile slurm \
    preprocessing                 # stage target: stop after QC
```

The wrapper resolves the dependency graph, submits each rule as its own SLURM
job, keeps a conda environment per tool, and resumes after a failure. The
`preprocessing` target stops after QC instead of running Hecatomb's full
taxonomic cascade, which would have taken days and duplicated what the phage
pipeline does better.

**This is the only QC the cohort received.** Worth saying because the next
wrapper also takes a `reads_dir`: **phage-analysis contains no QC rules** — no
trimming, no host removal — it uses reads only for coverage in Stage 6. (A third
pipeline we tried, `prophage-pipeline`, does bundle fastp, but that run was
abandoned at binning and never fed anything.)

#### Where are adapter removal, PhiX removal, quality trimming?

If you have run a QC pipeline before, you expect a list like *adapter removal →
primer removal → vector/PhiX removal → quality trimming → length filter → host
removal*. Only one of those appears as its own row in our tables, which looks
like steps are missing. Two different reasons, and one is a genuine gap.

**Most of them are inside fastp.** fastp is a single-pass tool that does what
older pipelines spread across several programs. Mapping the classic list onto
what ran:

| Classic QC step | Did it happen? | Where |
|---|---|---|
| Adapter removal | **Yes** | fastp `--detect_adapter_for_pe` — 3,740,339 reads trimmed |
| Quality trimming | **Yes** | fastp `--cut_tail` w=25, Q<15 — 418,402 reads discarded |
| Length filter | **Yes** | fastp `--length_required 90` — 2,144,966 discarded |
| Poly-G/X tail | **Yes** | fastp `--trim_poly_x` — 786,119 trimmed |
| Ambiguous-base filter | **Yes** | fastp, >2 N — 3,346 discarded |
| Duplicate removal | **Yes** | fastp `--dedup` — 34,359,370 discarded |
| Host removal | **Yes** | minimap2 + samtools — 581,127 pairs |
| **PhiX / vector removal** | **NO — skipped** | see below |
| Primer removal (Round A/B) | **NO — not applicable** | see below |

So the single `fastp` row in the attrition table is six classic steps collapsed
into one pass. That is a reporting artefact of how fastp works, not missing QC —
which is exactly why the per-flag breakdown below matters more than the per-step
totals.

**The real gap: PhiX and vector contamination.** Hecatomb *does* have this step,
but only on its other preprocessing route. The two routes are not minor variants:

| | `--preprocess paired` (what we ran) | `--preprocess roundAB` |
|---|---|---|
| Cleaning rules | **1** (`fastp_preprocessing`) | **6** discrete `bbduk.sh` passes |
| 5′ primer removal | — | `remove_5prime_primer` |
| 3′ contaminant | — | `remove_3prime_contaminant` |
| Primer-free adapter | — | `remove_primer_free_adapter` |
| Adapter-free primer | — | `remove_adapter_free_primer` |
| **PhiX + UniVec** | — | `remove_vector_contamination` |
| Quality/entropy | inside fastp | `remove_low_quality` (adds an entropy filter) |
| Everything after | identical | identical |

The `roundAB` route is built for Hecatomb's own Round A/B amplification protocol,
and its primer files in this install are hand-edited for lab-specific primers.
We chose `paired` because our reads are **simulated by InSilicoSeq** — they carry
no Round A/B primers, no PhiX spike-in, and no vector sequence, so five of those
six passes would have found nothing. The vector pass is the one that would be
mandatory on real sequencer output.

```bash
# what roundAB runs and we did not -- PhiX + NCBI UniVec
bbduk.sh in=R1.fastq in2=R2.fastq ref=vector_contaminants.fa \
    out=clean_R1.fastq out2=clean_R2.fastq stats=vector.stats \
    k=31 hammingdistance=1 ordered=t ow=t threads=16

# roundAB's quality pass -- note the entropy filter fastp has no equivalent for
bbduk.sh in=R1.fastq in2=R2.fastq qtrim=r maxns=2 \
    trimq=${QSCORE} minlength=${READ_MINLEN} \
    entropy=${ENTROPY} entropywindow=${ENTROPYWINDOW}
```

> **If you are adapting this to real data, add PhiX removal.** Illumina runs are
> routinely spiked with PhiX as a sequencing control, and it is not removed by
> every demultiplexer. It is cheap to screen for: `bbduk.sh` against the bundled
> `phix174_ill.ref.fa.gz`, or `--preprocess roundAB` if your library used that
> protocol. We skipped it because simulated reads cannot contain it — that
> justification does not transfer to a real cohort.
>
> The **entropy filter** in `remove_low_quality` is also worth knowing about:
> it discards low-complexity reads (homopolymer runs, simple repeats) that pass
> a quality threshold perfectly well. fastp has no direct equivalent.

#### The full rule list — nine steps, three of which remove reads

"Three QC steps" is a simplification, and worth unpacking, because Hecatomb's
`preprocessing` target actually ran **nine rules** per sample. Only three change
the read count; the rest reshape the data. Counted from the run's own benchmark
and log files:

| # | Rule | Tool | Jobs | Effect on read count |
|---|---|---|---|---|
| 1 | `fastp_preprocessing` | fastp | 100 | **removes** — trim, filter, dedup |
| 2 | `host_removal_mapping` | minimap2 + samtools | 300 | **removes** — drops host |
| 3 | `nonhost_read_repair` | bbtools `reformat.sh` | 100 | re-pairs — splits singletons/orphans back into R1/R2 |
| 4 | `nonhost_read_combine` | `cat` | 100 | re-joins paired + singleton streams per sample |
| 5 | `cluster_similar_sequences` | MMseqs2 linclust | 100 | **collapses** — 0 reads lost, counts preserved |
| 6 | `create_individual_seqtables` | seqkit + awk | 100 | reshapes to a per-sample sequence/count table |
| 7 | `merge_seq_table` | `mergeSeqTable.py` | 1 | merges 100 tables into one cohort seqtable |
| 8 | `calculate_gc` / `calculate_tet_freq` | — | 2 | annotates GC% and tetranucleotide frequency |
| 9 | `seq_properties_table` | — | 1 | assembles the final properties table |

Three observations for anyone reading this as a method:

- **`host_removal_mapping` shows 300 jobs for 100 samples** because the rule is
  three chained invocations (minimap2 → `samtools view` → `samtools fastq`), each
  logging separately. One logical step, three processes.
- **Steps 3–4 exist because host removal breaks pairing.** When one mate maps to
  human and the other doesn't, you get an orphan. `reformat.sh` splits those back
  into synchronised R1/R2 files and `cat` recombines them, so downstream
  assembly gets valid pairs. Skip this and MEGAHIT will complain about mismatched
  files. This is the step most hand-rolled QC scripts forget.
- **Steps 6–9 are bookkeeping, not cleaning.** They build Hecatomb's seqtable —
  the input to its taxonomic cascade, which we did not run. If you QC by hand for
  an assembly-based workflow like ours, you stop after step 4 (or 5).

A fourth note on the counts: `fastp_preprocessing` logged 87 jobs and
`old_fastp_preprocessing` 100. Both appear because the workflow retains a legacy
rule name; the 87 reflects a resumed run after a partial failure, not 13 skipped
samples. All 100 samples have complete fastp output, which the attrition table
below confirms.

#### How many reads each step removed

Of those nine rules, three move the number. VLP cohort, 100 samples, read pairs:

| Step | Tool | Pairs in | Pairs out | Removed | % removed | % of raw left |
|---|---|---|---|---|---|---|
| 0. input | — | — | 67,652,053 | — | — | 100.00% |
| 1. trim + filter + dedup | fastp 0.23.2 | 67,652,053 | 49,189,011 | 18,463,042 | **27.29%** | 72.71% |
| 2. host removal | minimap2 + samtools | 49,189,011 | 48,607,884 | 581,127 | **1.18%** | 71.85% |
| 3. dereplication | MMseqs2 linclust | 48,607,884 | 48,607,884 | 0 | **0.00%** | 71.85% |

**Nearly all the loss is one step, and nearly all of that is duplicates.**
Breaking fastp open (reads, not pairs — a pair can be hit by more than one rule):

| Sub-operation | Flag | Reads affected |
|---|---|---|
| Adapter trimmed | `--detect_adapter_for_pe` | 3,740,339 |
| Poly-X tail trimmed | `--trim_poly_x` | 786,119 |
| Discarded — low quality | Q<15, `--cut_tail` w=25 | 418,402 |
| Discarded — too short | `--length_required 90` | 2,144,966 |
| Discarded — too many N | >2 N | 3,346 |
| **Discarded — duplicate** | `--dedup --dup_calc_accuracy 4` | **34,359,370** |

So of the 27.29% that fastp removed, **duplicates account for ~93%**. Quality
filtering — the thing people usually picture when they hear "QC" — discarded
0.3% of reads. Trimming *modifies* reads without discarding them, which is why
the first two rows don't subtract from the total.

Three things a reader should take from this shape:

- **Dereplication removed nothing.** 0.00% is not a bug — MMseqs2 collapses reads
  onto representatives while carrying counts forward, so no read is lost. It
  compressed 48.6 M pairs to 7.1 M representatives (6.8×) for searching. Deduplication
  deletes; dereplication groups. (See the box after step 5.)
- **Host removal at 1.18% is the VLP prep working**, not a weak filter.
- **Quality is rarely where reads go.** If your own QC drops 30%, look at
  duplicates and length before you touch the quality threshold.

**Per-sample spread** (`results/qc/qc_attrition.tsv`, one row per sample) — useful
because a cohort mean hides the sample that failed:

| Stage | Min | Median | Max |
|---|---|---|---|
| Retained after fastp | 67.91% | 73.96% | 77.74% |
| Retained after host removal | 98.37% | 98.75% | 99.12% |
| **Retained overall** | **67.17%** | **72.91%** | **76.73%** |

Tight — no sample collapsed, which is what you want to confirm before analysing.
But note the group split: **HPV+ retained 70.12% vs HPV− 74.36%** (medians). QC
attrition is not independent of the grouping variable here, and that is the
coverage confound documented in the [results
summary](https://leran10.github.io/HPV_synthetic_data/) showing up as early as
the trimming step. Always break attrition down by group, not just by sample.

**All of it is in the repo:** [`results/qc/`](../results/qc) —
`qc_steps_vlp.tsv` (the first table), `qc_steps_detailed_vlp.tsv` (the fastp
breakdown), `qc_attrition.tsv` (per sample), `qc_steps_by_group_vlp.tsv` (by
group), with `_bulk` counterparts for the 80%-host cohort. Rebuilt from the tool
stats files by [`scripts/qc_steps_detailed.py`](../scripts/qc_steps_detailed.py)
and [`scripts/collect_qc_attrition.py`](../scripts/collect_qc_attrition.py).

#### Do you need Hecatomb?

| Your goal | Answer |
|---|---|
| Reproduce our exact read counts | **Yes** — the parameters and the reference both came from it |
| QC your own VLP virome | **No** — the three commands below are the whole of it |
| QC a bulk metagenome | Optional, but its masked reference starts to matter — see below |

Hecatomb's contribution is orchestration plus a pre-masked host reference. The
transferable knowledge is in the tools: how fastp finds adapters, why a host
reference gets masked, what dereplication preserves. That applies whichever
route you take.

#### Running the same QC without the wrapper

Per sample. Step 1b is marked because we skipped it and you probably should not:

```bash
# 1. adapters, quality, length, poly-X, duplicates
fastp -i VIR_041_R1.fastq.gz -I VIR_041_R2.fastq.gz \
      -o t_R1.fastq -O t_R2.fastq -j VIR_041.fastp.json -h VIR_041.fastp.html \
      --qualified_quality_phred 15 --length_required 90 --detect_adapter_for_pe \
      --cut_tail --cut_tail_window_size 25 --cut_tail_mean_quality 15 \
      --dedup --dup_calc_accuracy 4 --trim_poly_x --thread 12

# 1b. PhiX / vector screen -- WE SKIPPED THIS (simulated reads); add it for real data
bbduk.sh in=t_R1.fastq in2=t_R2.fastq \
         ref=phix174_ill.ref.fa.gz,vector_contaminants.fa \
         out=nophix_R1.fastq out2=nophix_R2.fastq stats=phix.stats \
         k=31 hammingdistance=1 ordered=t ow=t threads=16

# 2. host removal -- index once, reuse
minimap2 -t 16 -d host.idx GRCh38.fa.gz                      # once
minimap2 -ax sr -t 16 --secondary=no host.idx nophix_R1.fastq nophix_R2.fastq \
  | samtools view -f 4 -h \
  | samtools fastq -NO -1 clean_R1.fastq -2 clean_R2.fastq

# 3. REPAIR PAIRING -- do not skip this
# host removal orphans reads whose mate mapped to human; re-synchronise R1/R2
repair.sh in=clean_R1.fastq in2=clean_R2.fastq \
          out=fixed_R1.fastq out2=fixed_R2.fastq \
          outs=singletons.fastq repair -Xmx8g

# 4. dereplication -- SKIP unless a huge database search follows (see below)
mmseqs easy-linclust fixed_R1.fastq VIR_041_rep tmp \
    --kmer-per-seq-scale 0.3 -c 0.8 --cov-mode 1 \
    --min-seq-id 0.97 --alignment-mode 3 --threads 8
```

**Step 3 is the one people leave out.** `samtools view -f 4` filters reads
independently, so a pair with one human mate becomes an orphan — the R1 and R2
files end up different lengths and out of order. MEGAHIT and most aligners will
either error or silently misinterpret that. Hecatomb handles it with
`reformat.sh` + `cat` across two rules; `repair.sh` (same bbtools suite) does it
in one. You can also use `seqkit pair` or `fastq_pair`.

**Step 4 removed 0 reads** for us and changed nothing downstream — it exists to
shrink the query set for Hecatomb's 73 GB taxonomic cascade, which we never ran.
Skip it unless you are feeding a large protein database.

#### The one real choice: host removal

**Hecatomb is not a host-removal tool.** That step *is* minimap2 + samtools; what
the wrapper supplies is a GRCh38 with virus-like regions blanked, so genuine
viral reads are not discarded as human.

How much that masking matters depends entirely on your library prep, and our two
cohorts bracket the answer:

| Cohort | Host reads removed | Does masking matter? |
|---|---|---|
| **VLP** (what we analysed) | 581,127 / 49.2 M = **1.18%** | Barely — too little host to get wrong |
| **Bulk**, 80% host | **43.26%** | Yes — a bad mask here costs you real viral reads |

So for VLP data, map against any human reference you can download. For bulk
data, the masking is worth having (`hecatomb add-host` builds one).

Any of these is a defensible choice for the step:

| Tool | Approach | When to reach for it |
|---|---|---|
| **minimap2** + samtools | minimizer-seeded alignment | What we used; fast on short reads vs a 3 Gb genome |
| **bowtie2** / **bwa** + samtools | FM-index alignment | Identical logic, more familiar to most labs |
| **bbduk** / **bbsplit** | k-mer matching | Fastest; no alignment step |
| **Kraken2** (human DB) | k-mer classification | Already in your stack for taxonomy |
| **hostile** | wrapper, ships clean indexes | Purpose-built for host depletion — easiest start |

If I were starting fresh on someone else's VLP data, **hostile** is the shortest
path; **minimap2 + samtools** is what to use when you want to see and control
every flag.

> **Read the 1.18% as a QC result, not a formality.** A VLP prep *should* leave
> almost no host behind — that is the enrichment working. If your own VLP library
> comes back 40% human, the prep failed, and no amount of software fixes it. The
> bulk column is what an unenriched library looks like.

#### Two things specific to this install

`--preprocess paired` runs one fastp call; the alternative `roundAB` route runs
six sequential `bbduk.sh` passes (compared in full
[above](#where-are-adapter-removal-phix-removal-quality-trimming)), including the
PhiX/vector screen we skipped. And Hecatomb's per-step counting rules
(`04_summaries_optional.smk`) are commented out in our version, so the attrition
tables in [`results/qc/`](../results/qc) were rebuilt from the stats files by
[`scripts/collect_qc_attrition.py`](../scripts/collect_qc_attrition.py).

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

**Produced (VLP cohort, the one analysed).** 67,652,053 pairs → 49,189,011
(27.29% removed), of which duplicates were ~93%. Full per-flag breakdown in the
[attrition tables above](#how-many-reads-each-step-removed). Bulk cohort for
comparison: 67,479,078 → 54,164,669 (19.73%).

### 4. minimap2 ≥2.20 + samtools ≥1.13 — host removal

**What it does.** Discards reads of human origin.

**Mechanism.** minimap2 indexes the host genome by minimizers (a sampled subset of
k-mers) and seeds alignments from them, so it maps short reads against a 3 Gb
genome quickly. `samtools view -f 4` then keeps only records flagged *unmapped*.
Our reference was viral-masked; see the choices table above for why that mattered
little here and a lot on bulk data.

```bash
minimap2 -ax sr -t 16 --secondary=no ${HOST_INDEX} R1.fastq R2.fastq \
  | samtools view -f 4 -h \
  | samtools fastq -NO -1 unmapped_R1.fastq -2 unmapped_R2.fastq
```

**Produced.** VLP cohort (analysed): 581,127 of 49,189,011 reads removed =
**1.18%**. Bulk cohort: **43.26%**.

### 5. MMseqs2 12.113e3 — dereplication (optional; skipped effect here)

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

**Produced (VLP).** 48,607,884 pairs → 7,111,211 representatives: a **6.8×**
smaller search space with **zero reads lost**, which is why its row in the
attrition table reads 0.00% removed. (Bulk: 30,733,965 → 5,108,133, 6.0×.)

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

### — phage-analysis — the wrapper that runs steps 7–16

Same caveat as Hecatomb: the ten tools in stages 3–5 were not invoked separately.
They are rules in the [phage-analysis](https://github.com/megjohnson1999/phage-analysis)
Snakemake workflow, run once:

```bash
export PATH=${CONDA_BASE}/bin:${CONDA_BASE}/condabin:$PATH   # conda must be ON PATH
cd ${PROJECT_DIR}/pipelines/phage-analysis

snakemake \
    --snakefile workflow/Snakefile \
    --profile   profile/slurm/ \
    --configfile config/vlp_config.yaml \
    --set-resources \
        mmseqs_taxonomy:cpus_per_task=24  mmseqs_taxonomy:mem_mb=700000 \
        mmseqs_taxonomy:runtime=4320 \
        genomad_prediction:cpus_per_task=16  genomad_prediction:mem_mb=100000 \
        phabox_prediction:mem_mb=100000   vcontact3_taxonomy:mem_mb=200000 \
        iphop_single_prediction:mem_mb=150000  iphop_single_prediction:runtime=2880 \
        # ...16 rules pinned in total; see scripts/run_phage_analysis_vlp.sh
    --nolock
```

The config names the inputs:

```yaml
output_dir:    ${PROJECT_DIR}/data/vlp/phage_analysis/out
assembly_file: ${PROJECT_DIR}/data/vlp/phage_analysis/all_samples_contigs.fasta
reads_dir:     ${PROJECT_DIR}/data/vlp/READS
skip_mmseqs:   true     # see below
```

**One assembly file, not a directory.** The pipeline takes a single FASTA, so all
100 per-sample assemblies were concatenated with sample-prefixed contig names
(`VIR_041__k141_78`). That prefix is what lets abundance be attributed back to
samples in stage 6.

**Why every resource is pinned on the command line.** The profile's
`set-resources:` block is silently ignored by this Snakemake build — partly
because it writes runtimes as strings (`"24:00:00"`) where Snakemake wants integer
minutes. Four separate runs failed on resource defaults before we stopped trusting
the profile: a rule declaring `threads: 24` was handed 4 CPUs, two rules hit a
24-hour ceiling, and iPHoP was given 50 GB when it needed more. CLI
`--set-resources` does take effect, verified with `--dry-run`.

**MMseqs2 `nr` taxonomy is skipped deliberately.** The shared index is 2,125 GB
across 9 split files on a filesystem reading at **7–10 MB/s** — one pass would take
roughly 81 days. Three runs (27 h, 14 h, 10 h) were lost to it before we timed a
raw read and understood why. `skip_mmseqs: true` routes the rule to its
empty-output branch; taxonomy comes from PhaGCN and vConTACT3 instead, hosts from
CHERRY and iPHoP. **Measure your filesystem before blaming a tool.**

**Result.** 435 jobs, 3h51m, zero errors — after fifteen failed attempts. The six
fixes that got there are in the [README](../README.md#pipeline-notes-hard-won).

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
