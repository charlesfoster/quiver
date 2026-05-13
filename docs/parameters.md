# Parameters — HCV Quasispecies Pipeline

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
| `--reference_panel` | path | `assets/hcv_references.fasta` | HCV reference panel FASTA (238 sequences). Used for Round 1 competitive mapping and genotype classification. | Replace with a custom panel if working with non-standard genotypes or if the bundled panel is updated. |
| `--host_reference` | path | null (required) | Path to GRCh38 no-alt FASTA or a pre-built minimap2 `.mmi` index. Providing a pre-built `.mmi` skips indexing and saves ~10 min per run. | Always required unless the samples contain no human DNA (e.g., cell-culture only). |
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
| `--min_variant_depth` | integer | `20` | LoFreq DP (total depth) filter applied during variant filtering. Variants at positions with fewer than this many reads are discarded. | Increase for stricter variant calls; decrease only for very low-coverage samples. |

---

## Variant calling

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--min_call_af` | float | `0.005` | LoFreq lower-bound allele frequency for raw variant calls (0.5%). All calls above this threshold are emitted to `lofreq.vcf.gz`. | Decrease to 0.001 only if ultra-low-frequency variants are required; increases false-positive burden. |
| `--min_report_af` | float | `0.01` | Reporting AF threshold (1%). Applied by the variant filter step to produce `lofreq.filtered.vcf.gz`. | Increase to 0.02–0.05 for conservative clinical reporting. |
| `--min_mq` | integer | `20` | Minimum mapping quality for reads contributing to variant calls (`lofreq call-parallel --min-mq`). | Increase to 30 for stricter calls at the cost of read yield. |
| `--min_bq` | integer | `7` | Minimum base quality for variant calls (`lofreq call-parallel --min-bq`). | Rarely need changing for R10.4.1 data with HAC basecalls. |
| `--lofreq_sig` | float | `0.01` | LoFreq strand-bias significance threshold. | Lower to 0.001 if strand-biased false positives are a concern (e.g., known problematic homopolymers). |

---

## Depth normalisation

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--lofreq_max_depth` | integer | `5000` | rasusa depth cap for the LoFreq input BAM. LoFreq sensitivity plateaus above ~5,000×; higher depth increases false-positive rate and runtime. | Rarely need changing. Increase to 10,000 only if very deep samples show underdetection. |
| `--devider_max_depth` | integer | `1000` | rasusa depth cap for the DEVIDER input BAM. DEVIDER memory scales super-linearly above 1,000×. | Increase for higher sensitivity in high-diversity or low-coverage samples; decrease if DEVIDER runs out of memory. |
| `--rasusa_seed_lofreq` | integer | `42` | Random seed for rasusa subsampling of the LoFreq input. Fixed seed ensures reproducible subsampling. | Change only if reproducibility tests reveal a seed-specific bias. |
| `--rasusa_seed_devider` | integer | `43` | Random seed for rasusa subsampling of the DEVIDER input. | As above. |

---

## Haplotype reconstruction

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--devider_min_cov` | integer | `50` | DEVIDER `--min-cov`: minimum per-window depth required for DEVIDER to attempt reconstruction in that window. Windows below this threshold are skipped. | Decrease to recover partial haplotypes from low-coverage regions; increases noise. |
| `--devider_min_abund` | float | `0.25` | DEVIDER `--min-abund`: minimum fractional abundance for a reconstructed haplotype to be reported (25%). | Decrease to 0.05 to recover minor haplotypes; increases the risk of chimeric haplotype artefacts. |
| `--stitch_min_reads` | integer | `5` | Minimum number of reads spanning a window junction required for `bin/stitch_haplotypes.py` to link two adjacent DEVIDER windows into a single haplotype. | Increase for stricter stitching; decrease for lower-coverage samples where spanning reads are scarce. |

---

## Run-mode toggles

| Parameter | Type | Default | Description | When to change |
|---|---|---|---|---|
| `--run_clair3` | boolean | `false` | Enable optional Clair3 corroboration calling (AF ≥ 25%). Adds significant runtime and disk usage. Outputs written to `variants/<GT>/clair3/`. | Enable when independent corroboration of high-AF variants is required. |
| `--allow_conda_fallback` | boolean | `false` | Allow conda/micromamba to resolve environments when a container is unavailable. Required if Docker/Singularity is not available on the system. | Enable when running on systems without container support. Prefer the `conda_local` profile for a fully conda-native run. |
| `--use_hostile` | boolean | `false` | Use the `hostile` tool instead of raw minimap2 for host depletion. `hostile` is a thin wrapper around minimap2 with additional host-database options. | Enable if a specific hostile database is preferred over GRCh38. |

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
