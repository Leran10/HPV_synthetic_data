# Background fragments are fixed-length and non-contiguous, which biases assembly benchmarking

## Summary

`add_microbial_background()` and the host-DNA contamination path both sample background
sequence as **random fixed-length windows** (`fragment_length: int = 10000`) rather than
whole genomes. Because the windows are drawn from random offsets, a single source genome
is represented by several disconnected 10 kb pieces with gaps between them.

The consequence is that **non-viral content can never assemble beyond 10 kb, no matter
the coverage**, while viral genomes (taken whole from collections) assemble completely.
Any benchmark metric that rewards contiguity or completeness therefore favours viruses
for reasons that have nothing to do with the tool under test.

I hit this running a bulk-metagenome cohort through a prophage-detection pipeline, but
the concern is general: it affects the `benchmarking/assembly.py` numbers any user will
get on a `--no-vlp` cohort.

## Evidence

Cohort: 100 samples, `--no-vlp --contamination-level heavy`, novaseq 150 bp, 5-20x.
Measured read composition: bacterial 47.6% | host 40.1% | rRNA 6.5% | fungal 3.6% |
viral ~1.0%.

**1. Every reference fragment is exactly 10 kb.**

```
bacterial_fragments.fasta   n=120  min=10000  median=10000  max=10000
host_fragments.fasta        n=48   min=10000  median=10000  max=10000
fungal_fragments.fasta      n=15   min=10000  median=10000  max=10000
```

**2. Assemblies hit that ceiling exactly.** MEGAHIT (`--k-min 21 --k-max 141`, no length
cap) on one sample, longest contigs:

```
9999 9998 9998 9997 9997 9997 9997 9997 9996 9996 9996 9994
```

Thirteen contigs within 6 bp of 10,000. That is not coverage-limited fragmentation --
the assembler reconstructed each fragment end to end and stopped because nothing more
existed. Per-sample assemblies were 38-66 contigs totalling 230-306 kb.

**3. Binning is impossible, not merely hard.** DAS Tool: `No bins with bin-score >0.5
found`; MetaBAT produced an empty contigs2bin file. A 2-5 Mb genome represented by a
handful of non-adjacent 10 kb windows has no assembly path to a genome. Extra coverage
does not help.

**4. Viral content is unaffected**, because collection genomes are supplied whole --
e.g. 7 contigs in the 7-8.5 kb band matching HPV (~7.9 kb), one contig at 13.8 kb.

## Why this is a benchmarking problem

`benchmarking/assembly.py` reports per-genome recovery/completeness/identity and
contiguity (N50/L50). On a `--no-vlp` cohort:

- **N50 / L50 are not interpretable.** The distribution is capped at 10 kb for the
  majority of the sample by construction. A real bulk metagenome would have bacterial
  contigs in the 100 kb+ range pulling N50 up.
- **Completeness is structurally bimodal** -- viral genomes can reach 100%, non-viral
  content cannot exceed `fragment_length / genome_length` (~0.2-0.5% for a bacterium).
- **Any viral-vs-non-viral comparison is an artefact.** "Recovered 85% of viral genomes
  and 0% of bacterial genomes" measures the simulator, not the pipeline.
- **Prophage and binning tools cannot be benchmarked at all** on ViroForge output, since
  there are no host genomes to integrate into or bin.

`assembly.py` already scopes recovery to viral genomes (`# per-genome completeness (viral
only for recovery)`), which limits the damage -- but the contiguity metrics are computed
over the whole assembly, and nothing warns the user that the non-viral fraction is
synthetic ballast.

## Suggestions

Roughly in order of effort:

1. **Document the limitation where the numbers appear.** Add a note to the benchmarking
   report output and to `--no-vlp` help text: background sequence is fixed-length
   non-contiguous fragments; assembly contiguity metrics reflect `fragment_length`, not
   assembler performance; binning and prophage tools cannot be benchmarked on this data.

2. **Expose `fragment_length`.** It is currently a function default, not reachable from
   the CLI or batch YAML. Users cannot change it without editing source. Even just
   `--background-fragment-length` would let people match it to their use case.

3. **Report it in the metadata.** The composition TSV and `*_metadata.json` record
   abundances but not the fragment length used, so a cohort's assembly ceiling is not
   recoverable from its outputs. Add `background_fragment_length` to the metadata.

4. **Make contiguity metrics fragment-aware.** Either compute N50 over the viral fraction
   only, or report viral and non-viral contiguity separately, so the 10 kb ceiling does
   not silently depress a whole-assembly N50.

5. **Optional: a whole-genome background mode.** For users benchmarking binning or
   prophage detection, allow `fragments_path` to be supplied as complete genomes with
   `fragment_length=None` meaning "use the full sequence". Expensive in reference-data
   size, but it is the only way those tool classes become testable. Could ship as an
   opt-in download rather than bundled.

Items 1-3 are cheap and would have saved me a failed 5-hour pipeline run; 4-5 are the
real fix if benchmarking non-viral recovery is ever in scope.

## Environment

ViroForge 0.20.0, `--no-vlp --contamination-level heavy`, platform novaseq, collections
16 and 24, 100 samples. `fragment_length` left at its default throughout -- not
overridden anywhere in my config.
