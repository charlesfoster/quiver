# Usage — HCV Quasispecies Pipeline

This document explains how to install, configure, and run `hcv-quasi`.
For architecture and design decisions see [CLAUDE.md](../CLAUDE.md).
For parameter descriptions see [docs/parameters.md](parameters.md).

---

## Requirements

- **Nextflow ≥ 24.10.5**

  ```bash
  curl -s https://get.nextflow.io | bash
  # or via conda/micromamba:
  micromamba install -c bioconda nextflow=24.10.5
  ```

- **Container runtime or conda** — choose one depending on your environment:
  - Docker (any recent version) — required for `docker` and `docker_mac` profiles.
  - Singularity / Apptainer — required for `katana` and `gadi` profiles.
  - conda with micromamba ≥ 1.5 — required for `conda` profile.

- **Host reference genome** (GRCh38 no-alt) for host depletion:

  ```bash
  # Download from NCBI (one-time, ~3 GB compressed):
  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/\
  GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
    -O /reference/GRCh38.fa.gz
  ```

  Pass this path via `--host_reference /reference/GRCh38.fa.gz`.
  The pipeline builds a minimap2 index on first use and caches it.
  This is an explicit local path; the pipeline does not currently auto-download GRCh38.

---

## Profiles

Select a profile with `-profile <name>`. Profiles are defined in `conf/` and documented in [docs/configuration.md](configuration.md).

### `docker` — Mac or Linux workstation with Docker

Runs all processes locally using Docker. Suitable for development and small datasets.

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --host_reference /path/to/GRCh38.fa.gz \
    --outdir results/run1
```

### `docker_mac` — Apple Silicon Docker

Runs locally using Docker Desktop on Apple Silicon. This inherits the Docker profile but sets LoFreq to serial calling to avoid `lofreq call-parallel` OOM/SIGKILL failures under amd64 emulation.

```bash
nextflow run main.nf \
    -profile docker_mac \
    --input samplesheet.csv \
    --host_reference /path/to/GRCh38.fa.gz \
    --outdir results/run1
```

### `conda` — Mac or Linux without Docker

Uses micromamba to resolve the conda environments. Preferred when Docker/Singularity is not available.

```bash
nextflow run main.nf \
    -profile conda \
    --input samplesheet.csv \
    --host_reference /path/to/GRCh38.fa.gz \
    --outdir results/run1
```

### `katana` — UNSW Katana (SLURM + Singularity)

Submits jobs to Katana using SLURM. Singularity images are pulled automatically and cached under `/srv/scratch/${USER}/.singularity_cache`.

```bash
nextflow run main.nf \
    -profile katana \
    --input samplesheet.csv \
    --host_reference /srv/scratch/${USER}/ref/GRCh38.fa.gz \
    --outdir results/run1
```

To override the SLURM account, add `--clusterOptions '--account=<your_account>'` or edit `conf/katana.config`.

### `gadi` — NCI Gadi (SLURM + Singularity)

Requires your Gadi project code. Storage flags are set automatically from the project code.

```bash
nextflow run main.nf \
    -profile gadi \
    --gadi_project <PROJECT> \
    --input samplesheet.csv \
    --host_reference /g/data/<PROJECT>/ref/GRCh38.fa.gz \
    --outdir results/run1
```

### `test` — bundled synthetic test data

Runs the full pipeline on tiny synthetic data included in `test_data/`. No external inputs required. Completes in under 5 minutes.

```bash
nextflow run main.nf -profile test --outdir results/test
```

---

## Providing the samplesheet

The samplesheet is a CSV with three columns:

```
sample_id,fastq,metadata_json
P001,/absolute/path/to/P001.fastq.gz,
P002,/absolute/path/to/P002.fastq.gz,{"collection_date":"2026-01-15"}
```

- `sample_id` must match `^[A-Za-z0-9._-]+$` and be unique within the file.
- `fastq` must be an absolute path to a readable gzip-compressed FASTQ file.
- `metadata_json` is optional; omit or leave blank.

Full schema: [docs/configuration.md](configuration.md).

---

## Smoke test on real data

For manual validation after installation:

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --reference_panel assets/hcv_references.fasta \
    --host_reference /path/to/GRCh38.fa.gz \
    --outdir results/smoke \
    -resume
```

Validation criteria: variants at known polymorphic positions called within ±2% of expected AF; primary genotype matches prior typing. See [docs/testing.md](testing.md) Section 9.4 for full manual criteria.

---

## Resuming failed runs

Nextflow caches completed process outputs by content hash. Use `-resume` to skip already-completed steps:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results/run1 -resume
```

The work directory (default `./work/`) must still exist. Delete it to force a full re-run.

---

## Resource customisation

Override resource ceilings at the command line:

```bash
nextflow run main.nf -profile docker \
    --max_cpus 8 \
    --max_memory '32.GB' \
    --max_time '12.h' \
    --input samplesheet.csv \
    --outdir results/run1
```

Per-process resource labels (`process_low`, `process_medium`, `process_high`, `process_high_memory`) are defined in `conf/base.config`. To change resources for a specific process, add a `withName` block to a custom config file and pass it with `-c my.config`.

---

## Troubleshooting

### No HCV reads detected

The pipeline emits a `NO_HCV_DETECTED` flag for any sample where fewer than `params.min_round1_mapped` (default 100) reads map to the reference panel in Round 1. Possible causes: incorrect FASTQ file, wrong reference panel, or very low viral load. The sample is excluded from all downstream steps and a QC-only report is generated.

### Low coverage

If mean coverage after Round 2 mapping is below `params.min_mean_coverage` (default 100×), the pipeline sets a `LOW_COVERAGE` flag, skips DEVIDER haplotype reconstruction, and continues with LoFreq variant calling. LoFreq results from low-coverage samples should be interpreted with caution.

### Mixed infection detected

When the fraction of reads mapping to a secondary genotype meets or exceeds `params.min_secondary_fraction` (default 5%), the pipeline sets `is_mixed = true` and branches into per-genotype sub-workflows. Each genotype produces its own consensus, variants, and haplotypes. The `genotyping/` output directory contains the full classification summary.

Adjust the threshold with `--min_secondary_fraction`. Setting it below 0.05 may increase false-positive mixed calls from cross-contamination. See design decision D7 in [CLAUDE.md](../CLAUDE.md).

### DEVIDER fails for a sample

DEVIDER failure is non-fatal. The pipeline emits a `devider.failed` marker in the haplotypes directory, continues to reporting, and notes the failure in the per-sample HTML report. Check the Nextflow `.nextflow.log` and the `work/` directory for the DEVIDER process for the underlying error.

### Container pull failure

By default, failed container pulls are retried once. Set `params.allow_conda_fallback = true` (or use `-profile conda`) to fall back to conda. On Apple Silicon, if an x86 image is required and Docker is configured with Rosetta, the pull should succeed but will run under emulation.
