# Important Considerations — DEVIDER Haplotype Reconstruction

This document explains the factors that determine whether DEVIDER recovers multiple
distinct haplotypes or collapses to a single consensus sequence. Understanding these
factors is essential for interpreting DEVIDER output and for tuning the pipeline when
haplotype recovery is unsatisfactory.

---

## How DEVIDER reconstructs haplotypes

DEVIDER is a de Bruijn graph-based haplotype assembler designed for noisy long reads
from viral quasispecies. The core idea:

1. Reads are broken into overlapping k-mers.
2. A de Bruijn graph is built where nodes are k-mers and edges represent k-1 overlaps.
3. Each variant position (SNP) where two alleles are present in the population creates
   a **bubble** in the graph: the reference allele and the alternate allele produce two
   distinct k-mers, diverging at that position and re-converging after it.
4. Haplotypes are **paths through the graph**. A haplotype that carries a specific
   combination of alleles across N variant positions corresponds to a specific path
   that traverses N bubbles in sequence.
5. Paths are weighted by read support. DEVIDER reports paths whose abundance exceeds
   `--min-abund` and whose coverage exceeds `--min-cov`.

---

## Why de Bruijn graphs collapse to a single haplotype

### The combinatorial explosion

With N variant positions, there are up to 2^N possible allele combinations, each
corresponding to a distinct haplotype. As N grows, the graph rapidly becomes a dense
tangle of interconnected bubbles. The number of valid paths through the graph explodes
exponentially, and the read support for any single minority path drops below the
signal-to-noise threshold.

### The key constraint: reads must span multiple variant positions

For DEVIDER to phase two variant positions together (i.e., determine that allele A at
position 1 co-occurs with allele B at position 2 on the same haplotype), at least some
reads must physically span both positions in a single alignment. A read that covers only
one variant position provides depth at that position but contributes **no phasing
information** linking it to adjacent variants.

When variant positions are dense relative to the read length:

- Most reads span only 1–2 bubbles rather than 5–10.
- The graph has many bubbles but few reads bridging them.
- The assembler cannot confidently assign a read to one path vs. another.
- The result is graph collapse: DEVIDER emits only the majority (consensus) path because
  minority paths lack sufficient multi-spanning read support to be resolved.

### Sequencing error compounds the problem

ONT R10.4.1 reads carry approximately 2–5% per-base error. At any given genomic
position, up to 5 reads in every 100 will carry an incorrect base — not because a minor
variant is present, but due to sequencing noise. From the graph's perspective, a real
1% minor allele and a 1% noise allele are indistinguishable without additional statistical
context. When many noisy low-AF variants are included in the VCF passed to DEVIDER, the
graph acquires spurious bubbles that fragment haplotype paths and prevent clean traversal.

### Summary of collapse conditions

| Condition | Effect on graph | Outcome |
|---|---|---|
| Too many variant positions (high variant density relative to read length) | Many bubbles, few bridging reads | Single consensus haplotype emitted |
| Very low-AF variants included (< 5%) | Noise bubbles indistinguishable from signal | Fragmented or collapsed graph |
| Short reads (< genome / variants_per_read threshold) | Bubbles not bridged by individual reads | Haplotypes cannot be phased |
| Very low coverage at a window | Insufficient reads to traverse minority paths | Window skipped (`--min-cov`) |

---

## Read length

**This is the single most important factor for haplotype recovery.**

The HCV genome is approximately 9,286 bp. A read of length L bp will on average span
`L × N / G` variant positions, where N is the total number of variants and G is the
genome length.

| Read length | SNPs spanned (26 SNPs / 9286 bp genome) | SNPs spanned (100 SNPs / 9286 bp genome) |
|---|---|---|
| 1,000 bp | ~2.8 | ~10.8 |
| 2,000 bp | ~5.6 | ~21.5 |
| 4,000 bp | ~11.2 | ~43.1 |
| 9,000 bp | ~25.2 | ~96.7 |

With ~26 SNPs and 4,000 bp reads (the pipeline default), each read spans roughly 11
variant positions on average — sufficient for DEVIDER to phase the major and minor
haplotypes in a typical quasispecies sample.

If variant density is high (many SNPs), or reads are short (low yield from the run,
or a lot of reads below the length threshold), fewer variants are phased per read and
haplotype recovery degrades. The pipeline default of `--devider_min_read_length 4000`
removes short reads that would contribute depth without contributing phasing signal,
improving the quality of the de Bruijn graph at the cost of read yield.

**If DEVIDER recovers only one haplotype**, the first diagnostic step is to check
the number of SNPs in the VCF passed to DEVIDER (see
`variants/<GT>/lofreq.filtered.vcf.gz`) and the read length distribution after the
length filter (check nanoq QC outputs). If SNP count is high (> ~50) or read lengths
are predominantly shorter than expected, the parameters below should be adjusted.

---

## Variant quality and the VCF passed to DEVIDER

DEVIDER does not use the raw BAM reads alone; it also takes a VCF of variant positions
to define the bubbles in its graph. The pipeline filters the LoFreq output to
`--devider_min_af` (default 5%) before passing it to DEVIDER. This is intentionally
higher than the reporting threshold (`--min_report_af`, default 1%) because:

- LoFreq calls variants down to 0.5% AF. Many of these are real quasispecies variants,
  but at 1–4% AF they are close to the ONT error floor.
- Including them creates additional bubbles in the DEVIDER graph that are not supported
  by enough reads to be phased reliably.
- At 5% AF, variants have a clear signal-to-noise advantage over ONT background error
  (~2–5%), so bubbles are well-supported and graph traversal is more confident.

If you observe only one haplotype and expect multiple, **raising `--devider_min_af`**
(e.g., to 0.10 or 0.15) reduces the number of variants in the graph and can allow
DEVIDER to resolve haplotypes it previously collapsed. The tradeoff is that this
excludes lower-frequency variants from phasing, so minor haplotypes differing only at
low-AF positions may not be distinguished.

---

## Tunable parameters and their effects

| Parameter | Default | What it controls | Effect of increasing | Effect of decreasing |
|---|---|---|---|---|
| `--devider_min_read_length` | `4000` | Minimum read length for DEVIDER input (applied before depth cap) | Fewer reads, but each spans more variant positions; cleaner graph | More reads enter, but short reads add depth without phasing; graph may become noisier |
| `--devider_min_af` | `0.05` | Minimum AF for variants in the VCF passed to DEVIDER | Fewer, higher-confidence variant positions; simpler graph; may miss low-AF haplotype differences | More variant positions; graph is more complex; more likely to collapse at low coverage |
| `--devider_min_cov` | `10` | Minimum per-window depth for DEVIDER to attempt reconstruction | Fewer windows attempted; only high-coverage windows produce haplotypes | More windows attempted; very low-coverage windows are included at risk of noise |
| `--devider_min_abund` | `0.25` | Minimum haplotype abundance to report (literal %, so 0.25 = 0.25%) | Fewer haplotypes reported; minor variants suppressed | More haplotypes reported; increases risk of chimeric or noise haplotypes |
| `--devider_max_depth` | `5000` | rasusa depth cap before DEVIDER; excess reads subsampled | Lower cap → smaller graph → faster; may miss low-coverage minority haplotypes | Higher cap → more reads → better minority coverage but higher memory/runtime |

---

## Resuming a Nextflow run with different DEVIDER settings

Nextflow caches completed process outputs by content hash. When you change only
DEVIDER-relevant parameters and re-run with `-resume`, only the affected downstream
steps are re-executed — all upstream steps (read QC, host depletion, Round 1/2 mapping,
LoFreq calling) are retrieved from cache.

### Which steps re-run for each parameter change

| Parameter changed | Steps re-run |
|---|---|
| `--devider_min_read_length` | `SAMTOOLS_LENGTH_FILTER → RASUSA_ALN → LOFREQ_PREPROCESS (DEVIDER) → DEVIDER → FORMAT_HAPLOTYPES → SAMPLE_REPORT` |
| `--devider_min_af` | `FILTER_VCF_FOR_DEVIDER → DEVIDER → FORMAT_HAPLOTYPES → SAMPLE_REPORT` |
| `--devider_min_cov`, `--devider_min_abund` | `DEVIDER → FORMAT_HAPLOTYPES → SAMPLE_REPORT` |
| `--devider_max_depth` | `RASUSA_ALN → LOFREQ_PREPROCESS (DEVIDER) → DEVIDER → FORMAT_HAPLOTYPES → SAMPLE_REPORT` |

### Example: retry with stricter read length and higher AF threshold

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --outdir results/run1 \
    --devider_min_read_length 6000 \
    --devider_min_af 0.10 \
    -resume
```

Nextflow will reuse the Round 2 BAM and LoFreq VCF from the previous run and re-execute
only the DEVIDER-related steps. The `results/run1` output directory will be updated in
place with the new haplotype outputs.

### Example: retry with relaxed coverage floor

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --outdir results/run1 \
    --devider_min_cov 5 \
    -resume
```

### Keeping multiple parameter sets

If you want to compare results without overwriting, use a different `--outdir` for each
attempt. Nextflow's work directory (`./work/`) is shared between runs, so cached process
outputs are reused regardless of `--outdir`.

```bash
# Conservative (default)
nextflow run main.nf -profile docker --input samplesheet.csv \
    --outdir results/conservative -resume

# Relaxed read length, stricter AF
nextflow run main.nf -profile docker --input samplesheet.csv \
    --outdir results/relaxed_length \
    --devider_min_read_length 2000 --devider_min_af 0.10 -resume
```

---

## Practical checklist when DEVIDER recovers only one haplotype

1. **Check variant count.** Look at `variants/<GT>/lofreq.filtered.vcf.gz`. If it
   contains > ~50 variants, try raising `--devider_min_af` to 0.10 to reduce graph
   complexity.

2. **Check read length distribution.** Review the nanoq QC report
   (`qc/posthost/<sample_id>/nanoq.json`). If the N50 read length is well below 4,000 bp
   across the whole sample, the sequencing run may not have yielded reads long enough
   to phase across multiple variant positions. In this case, lowering
   `--devider_min_read_length` may increase read yield but will not necessarily improve
   haplotype resolution.

3. **Check coverage.** Review mosdepth output
   (`mapping/<GT>/qc/coverage/<GT>/mosdepth.summary.txt`). If mean coverage is below
   ~50×, the minority haplotype may genuinely not have enough reads to be distinguished
   from noise. DEVIDER requires sufficient depth on each haplotype independently.

4. **Check `--devider_min_abund`.** The default is 0.25% (very permissive). If the
   expected minor haplotype is at, say, 5% frequency, this threshold is not the
   limiting factor. If no haplotypes at all are reported, the collapse has occurred
   at the graph level (points 1–3 above), not at the abundance filter.

5. **Inspect the DEVIDER output directory directly.**
   `haplotypes/<GT>/devider/` contains DEVIDER's raw output. The presence of a
   `devider.failed` marker file indicates DEVIDER exited non-zero; check
   `.nextflow.log` and the corresponding `work/` subdirectory for the error message.
