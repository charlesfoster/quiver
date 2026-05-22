# Parameters — QuIVER

All parameters live under `params { }` in `nextflow.config`.
Defaults are documented in [docs/configuration.md](configuration.md).
Profile-specific overrides (e.g., resource ceilings for Katana, Gadi) are described in [docs/configuration.md](configuration.md) and applied in `conf/`.

Pass any parameter at the command line: `--parameter_name value`.
To set many parameters for a site, write a params JSON file and pass it with `-params-file my_params.json`.

---

## Inputs

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--input` | path | required | Path to the samplesheet CSV. Must contain `sample_id`, `fastq`, and optionally `metadata_json` columns. | Always required. |
| `--reference_panel` | path | `assets/hcv_references.fasta` | HCV reference panel FASTA (186 named-subtype sequences). Used for Round 1 competitive mapping and subtype classification. All sequences must follow the `<subtype>_<accession>` naming convention (e.g. `1a_M62321.1`). | Replace with a custom panel if working with non-standard genotypes or if the bundled panel is updated. |
| `--host_reference` | path | null | Path to a local GRCh38 no-alt FASTA or a pre-built minimap2 `.mmi` index. Only required with `--use_minimap2` or `--use_hostile`. Not needed for the default nohuman mode. | Provide when using `--use_minimap2` or `--use_hostile`. A pre-built `.mmi` skips indexing and saves ~10 min per run. |
| `--outdir` | path | `results` | Output directory. Created if it does not exist. | Change to avoid overwriting a previous run. |

---

## Read filtering

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--min_length` | integer | `200` | Minimum read length in bp (chopper `--minlength`). Reads shorter than this are discarded. | Increase to 500 bp if haplotype phasing quality is more important than read yield. |
| `--max_length` | integer | `10000` | Maximum read length in bp (chopper `--maxlength`). HCV is 9,646 bp; reads above this length are overwhelmingly chimeric. | Rarely need changing. |
| `--min_qual` | integer | `8` | Minimum mean Phred quality score (chopper `-q`). | Increase to 10 for higher-confidence basecalls (newer R10.4.1 v5 models). Decrease only if yield is critically low. |

---

## Genotyping and mixed-infection detection

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--min_secondary_fraction` | float | `0.05` | Minimum fraction of reads mapping to a secondary genotype that triggers the mixed-infection flag and per-genotype branching. | Decrease below 0.05 to catch low-level co-infections (increases false-positive risk from cross-contamination). Increase if minor-genotype branches are not scientifically required. |
| `--ambiguous_delta_as` | integer | `20` | Alignment score gap between a read's best and second-best mapping hit below which the read is called ambiguous and excluded from genotype assignment. | Increase if too many reads are being placed in the ambiguous pool. Decrease for stricter assignment. |
| `--min_round1_mapped` | integer | `100` | Minimum number of primary mapped reads in Round 1. Samples below this threshold are flagged `NO_HCV_DETECTED` and excluded from downstream analysis. | Decrease for very low-titre samples (with caution — downstream steps may produce unreliable results). |

---

## Coverage thresholds

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--min_mean_coverage` | integer | `100` | Mean genome coverage below which DEVIDER is skipped and a `LOW_COVERAGE` flag is set. LoFreq still runs. | Decrease to 50 if haplotype reconstruction is not required and low-coverage variant calls are acceptable. |
| `--min_consensus_cov` | integer | `10` | Per-position coverage below which a position is masked with N in the per-genotype consensus. Masked positions are excluded from Round 2 variant calling. | Increase for higher confidence consensus; decrease to recover more consensus bases from low-coverage samples. |
| `--min_consensus_identity` | float | `0.90` | Pairwise nucleotide identity threshold between the polished sample consensus and its dominant panel reference (minimap2 asm5). Below this, a `DIVERGENT_CONSENSUS` flag is set. The sample continues processing normally — the flag is informational. | Decrease to 0.88 for genotype-6 samples where some subtypes are inherently more distant from panel references. Increase to 0.93 for early-warning QC in well-characterised cohorts. Note: if `LOW_COVERAGE_CONSENSUS` also fired (many N positions), identity will be artificially low — consider both flags together. |
| `--min_variant_depth` | integer | `20` | LoFreq DP (total depth) filter applied during variant filtering. Variants at positions with fewer than this many reads are discarded. | Increase for stricter variant calls; decrease only for very low-coverage samples. |

---

## Variant calling

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--min_call_af` | float | `0.005` | LoFreq lower-bound allele frequency for raw variant calls (0.5%). All calls above this threshold are emitted to `lofreq.vcf.gz`. | Decrease to 0.001 only if ultra-low-frequency variants are required; increases false-positive burden. |
| `--min_report_af` | float | `0.01` | Reporting AF threshold (1%). Applied by the variant filter step to produce `lofreq.filtered.vcf.gz`. | Increase to 0.02–0.05 for conservative clinical reporting. |
| `--min_mq` | integer | `20` | Minimum mapping quality for reads contributing to LoFreq variant calls. | Increase to 30 for stricter calls at the cost of read yield. |
| `--min_bq` | integer | `7` | Minimum base quality for LoFreq variant calls. | Rarely need changing for R10.4.1 data with HAC basecalls. |
| `--min_alt_bq` | integer | `7` | Minimum base quality for LoFreq alternate-allele calls. | Rarely need changing for R10.4.1 data with HAC basecalls. |
| `--lofreq_sig` | float | `0.01` | LoFreq strand-bias significance threshold. | Lower to 0.001 if strand-biased false positives are a concern (e.g., known problematic homopolymers). |
| `--lofreq_pp_threads` | integer | `8` | Number of LoFreq `call-parallel` workers. Values of 1 use serial `lofreq call`; the `docker_mac` profile sets this to 1. | Lower to 1 on Apple Silicon Docker or other environments where `call-parallel` is unstable. |
| `--max_sb` | integer | `200` | Maximum strand-bias Phred score (INFO/SB). Variants above this threshold are discarded during filtering. | Raise to be more permissive; lower to remove more strand-biased calls. |

---

## Depth normalisation

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--lofreq_max_depth` | integer | `5000` | rasusa depth cap for the LoFreq input BAM. LoFreq sensitivity plateaus above ~5,000×; higher depth increases false-positive rate and runtime. | Rarely need changing. Increase to 10,000 only if very deep samples show underdetection. |
| `--devider_max_depth` | integer | `5000` | rasusa depth cap for the DEVIDER input BAM. | Decrease if DEVIDER runs out of memory on very deep samples. |
| `--rasusa_seed_lofreq` | integer | `42` | Random seed for rasusa subsampling of the LoFreq input. Fixed seed ensures reproducible subsampling. | Change only if reproducibility tests reveal a seed-specific bias. |
| `--rasusa_seed_devider` | integer | `43` | Random seed for rasusa subsampling of the DEVIDER input. | As above. |

---

## Haplotype reconstruction

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--devider_min_read_length` | integer | `4000` | Minimum read length (bp) retained for DEVIDER input. Applied before rasusa subsampling so the depth cap is drawn from long reads only. | Decrease to 2000 if read length distribution is short and DEVIDER yields few haplotypes. |
| `--devider_min_cov` | integer | `10` | DEVIDER `--min-cov`: minimum per-window depth required for DEVIDER to attempt reconstruction in that window. Windows below this threshold are skipped. | Increase for stricter reconstruction; decrease to attempt haplotypes from very low-coverage regions. |
| `--devider_min_abund` | float | `0.25` | DEVIDER `--min-abund`: minimum haplotype abundance to report. This is a literal percent value: `0.25` means 0.25%, not 25%. | Decrease to recover minor haplotypes at higher noise risk. |
| `--devider_min_af` | float | `0.05` | Minimum allele frequency for variants in the VCF passed to DEVIDER for phasing. Higher than `min_report_af` to prevent low-AF noise from saturating DEVIDER's graph. | Raise if DEVIDER produces fragmented haplotypes due to noisy SNPs; lower to include more variants in phasing. |

---

## Run-mode toggles

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--run_clair3` | boolean | `false` | Enable optional Clair3 corroboration calling (AF ≥ 25%). Adds significant runtime and disk usage. Outputs written to `variants/<GT>/clair3/`. | Enable when independent corroboration of high-AF variants is required. |
| `--allow_conda_fallback` | boolean | `false` | Allow conda/micromamba to resolve environments when a container is unavailable. | Enable when running on systems without container support. Prefer the `conda` profile for a fully conda-native run. |
| `--skip_host_depletion` | boolean | `false` | Bypass host depletion entirely. Not recommended for clinical samples. | Enable only for cell-culture or purely viral samples. |
| `--use_nohuman` | boolean | `false` | Explicitly select nohuman (Kraken2) host depletion. This is the default when no host depletion flag is set. | Rarely needed; the default already uses nohuman. |
| `--use_minimap2` | boolean | `false` | Use minimap2 alignment against GRCh38 for host depletion. Requires `--host_reference` or will auto-download GRCh38. | Use when precise alignment-based depletion is preferred over Kraken2. |
| `--use_hostile` | boolean | `false` | Use the `hostile` tool for host depletion. Requires `--host_reference`. | Use if a specific hostile database is preferred. |

---

## Cluster / HPC parameters

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--gadi_project` | string | null | NCI Gadi project code. Required when using `-profile gadi`. Sets SLURM account and storage flags. | Required for all Gadi runs. |

---

## Resources

These parameters set ceiling values that feed into Nextflow's `resourceLimits` and are applied across all processes via `conf/base.config`. Individual process resource requests (defined by `process_low` / `process_medium` / `process_high` / `process_high_memory` labels) are capped at these values.

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--max_cpus` | integer | `16` | Maximum CPUs any single process may request. | Lower on shared workstations; raise on large HPC nodes. Profile-specific overrides are already set in `conf/`. |
| `--max_memory` | string | `'64.GB'` | Maximum memory any single process may request. Format: `'<N>.GB'`. | Raise if DEVIDER or host-depletion processes OOM. Lower on memory-limited systems. |
| `--max_time` | string | `'24.h'` | Maximum wall-clock time any single process may request. | Raise to `'48.h'` for very deep or multi-genotype samples on HPC (profiles `katana` and `gadi` already set 48 h). |
