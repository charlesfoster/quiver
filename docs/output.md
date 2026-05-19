# Output — QuIVER

All outputs are written under `--outdir` (default `results/`).
Per-sample outputs are nested under `results/<sample_id>/`.
Pipeline execution metadata is written to `results/pipeline_info/`.
Run-level reports are written to `results/reports/`.

For the steps that produce each file see [docs/data_flow.md](data_flow.md).

---

## `results/<sample_id>/qc/`

Quality metrics for raw, filtered, and host-depleted reads.

| Path | Format | Contents | Absent when |
|---|---|---|---|
| `qc/raw/<sample_id>/NanoPlot-report.html` | HTML | NanoPlot read-length and quality visualisations for the raw FASTQ | Input FASTQ is empty (`EMPTY_INPUT` flag) |
| `qc/raw/<sample_id>/nanoq.json` | JSON | nanoq read statistics (count, N50, median quality) for the raw FASTQ | Input FASTQ is empty |
| `qc/posthost/<sample_id>/nanoq.json` | JSON | nanoq statistics after host depletion | Input FASTQ is empty |

---

## `results/<sample_id>/genotyping/`

Round 1 competitive mapping results and genotype classification.

| Path | Format | Contents | Absent when |
|---|---|---|---|
| `<sample_id>.round1.bam` | BAM | Sorted, indexed Round 1 alignment against the 238-sequence reference panel | Never absent if the sample reached this step |
| `<sample_id>.round1.bam.bai` | BAI | Index for the Round 1 BAM | As above |
| `<sample_id>.genotype_assignments.tsv` | TSV | One row per primary alignment: `read_id`, `ref_id`, `subtype`, `genotype`, `AS`, `mapq` | Fewer than `min_round1_mapped` reads mapped (`NO_HCV_DETECTED`) |
| `<sample_id>.genotype_summary.json` | JSON | Aggregated genotype fractions, `is_mixed` flag, `branches_to_run` list. Schema in [docs/configuration.md](configuration.md) | `NO_HCV_DETECTED` |

---

## `results/<sample_id>/consensus/<GT>/`

Per-genotype sample-specific consensus sequence. `<GT>` is the genotype label, e.g., `1a` or `3a`.

| Path | Format | Contents | Absent when |
|---|---|---|---|
| `consensus.fasta` | FASTA | Single-contig consensus named `<sample_id>_<GT>_consensus`; low-coverage positions masked with N | `NO_HCV_DETECTED`; `LOW_COVERAGE_CONSENSUS` flag |
| `consensus.fasta.fai` | FAI | samtools FASTA index | As above |
| `consensus.mmi` | MMI | minimap2 `map-ont` index of the consensus, used as the Round 2 reference | As above |
| `consensus.vcf.gz` | VCF | High-confidence variants (AF ≥ 50%, DP ≥ 10) applied during consensus build | As above |
| `mask.bed` | BED | Positions masked with N (coverage < `min_consensus_cov`) | As above |

---

## `results/<sample_id>/mapping/<GT>/`

Round 2 alignment of per-genotype reads against the sample-specific consensus.

| Path | Format | Contents | Absent when |
|---|---|---|---|
| `round2.bam` | BAM | Sorted, indexed Round 2 alignment | `NO_HCV_DETECTED`; `LOW_COVERAGE_CONSENSUS` |
| `round2.bam.bai` | BAI | Index | As above |
| `round2.flagstat.txt` | TXT | samtools flagstat output; `LOW_MAPPING_RATE` flag if <90% primary-mapped | As above |
| `qc/coverage/<GT>/mosdepth.summary.txt` | TXT | mosdepth per-contig coverage summary | As above |
| `qc/coverage/<GT>/<GT>.r2.regions.bed.gz` | BED.gz | Per-100-bp-window coverage from mosdepth | As above |

---

## `results/<sample_id>/variants/<GT>/`

Variant calling outputs for each genotype branch.

| Path | Format | Contents | Absent when |
|---|---|---|---|
| `lofreq.vcf.gz` | VCF.gz | Raw LoFreq variant calls (AF ≥ `min_call_af`); tabix-indexed | `NO_HCV_DETECTED`; `LOW_COVERAGE_CONSENSUS` |
| `lofreq.vcf.gz.tbi` | TBI | tabix index | As above |
| `lofreq.filtered.vcf.gz` | VCF.gz | LoFreq calls filtered to AF ≥ `min_report_af` and DP ≥ `min_variant_depth` | As above |
| `lofreq.filtered.vcf.gz.tbi` | TBI | tabix index | As above |
| `variants.tsv` | TSV | Tab-separated variant table: `CHROM`, `POS`, `REF`, `ALT`, `AF`, `DP`, `SB` | As above |
| `clair3/merge_output.vcf.gz` | VCF.gz | Optional Clair3 corroboration calls | Absent unless `--run_clair3`; also absent for the above failure modes |
| `clair3_concordance.tsv` | TSV | Concordance between LoFreq and Clair3 calls at overlapping positions | As above |

---

## `results/<sample_id>/haplotypes/<GT>/`

Haplotype reconstruction outputs from DEVIDER and post-hoc stitching.

| Path | Format | Contents | Absent when |
|---|---|---|---|
| `devider/` | directory | Raw DEVIDER output directory; contents vary by DEVIDER version and windowing. Parse by directory listing, not hard-coded filenames | `LOW_COVERAGE` flag (mean < `min_mean_coverage`); DEVIDER build failure |
| `devider.failed` | marker | Empty file indicating DEVIDER exited non-zero. The pipeline continues; downstream steps emit partial results | Absent when DEVIDER succeeded |
| `<sample_id>_<GT>_haplotypes.fasta` | FASTA | Haplotypes sorted by abundance (highest first) with annotated headers: `>ID abund:<pct> depth:<x> length:<bp>` | `LOW_COVERAGE`; DEVIDER failed; no haplotypes produced |
| `<sample_id>_<GT>_haplotype_map.tsv` | TSV | Maps clean sequential haplotype IDs back to DEVIDER's original headers | As above |
| `<sample_id>_<GT>_haplotype_report.json` | JSON | Per-haplotype summary: abundance, depth, length, and fallback flag if DEVIDER failed | As above |

---

## `results/reports/` and `results/<sample_id>/reports/`

HTML reports summarising the run.

| Path | Format | Contents | Absent when |
|---|---|---|---|
| `results/<sample_id>/reports/<sample_id>_summary.html` | HTML | Per-sample Jinja2 report: read counts at each stage, genotype summary, per-branch coverage, variant table, haplotype summary, all flags. NanoPlot PNGs embedded as base64 | Never absent once sample completes input validation |
| `results/<sample_id>/reports/<sample_id>_summary.json` | JSON | Machine-readable version of the per-sample report | As above |
| `results/reports/multiqc_report.html` | HTML | MultiQC aggregation across all samples: NanoPlot, mosdepth, samtools flagstat, bcftools stats | Only absent if no samples complete successfully |
| `results/reports/run_summary.html` | HTML | Tabular run-level summary listing each sample status (PASS / NO_HCV / LOW_COVERAGE / MIXED) with links to per-sample reports | As above |
| `results/reports/run_summary.json` | JSON | Machine-readable run summary | As above |

---

## `results/pipeline_info/`

Nextflow execution metadata, written automatically by Nextflow.

| Path | Format | Contents |
|---|---|---|
| `execution_timeline.html` | HTML | Per-process wall-clock timeline for the run |
| `execution_report.html` | HTML | Resource usage (CPU, memory, I/O) per process |
| `execution_trace.txt` | TSV | Detailed trace of all tasks: status, duration, CPU%, memory used |
| `pipeline_dag.html` | HTML | Directed acyclic graph of the pipeline topology |
