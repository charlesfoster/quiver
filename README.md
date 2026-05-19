# QuIVER

Nextflow DSL2 pipeline for reproducible HCV quasispecies analysis from ONT PromethION reads: genome-wide low-frequency variant calling (LoFreq, ≥1% AF) and global haplotype reconstruction (DEVIDER), with automatic detection and per-genotype branching for mixed-genotype infections.

---

## Quick start

```bash
# 1. Clone the repository
git clone https://github.com/charlesfoster/quiver.git
cd quiver

# 2. Print help
nextflow run main.nf --help

# 3. Run on the bundled test dataset
nextflow run main.nf -profile test --outdir results/test
```

---

## Requirements

| Requirement | Minimum version | Notes |
|---|---|---|
| Nextflow | 24.10.5 | Install via `curl -s https://get.nextflow.io | bash` |
| Docker | any recent | Required for `-profile docker` |
| Singularity / Apptainer | any recent | Required for `-profile singularity` (HPC) |
| conda or mamba | any recent | Required for `-profile conda` |

Disk: allow approximately 3× the size of your input FASTQ files for intermediate BAMs plus final outputs.
The DEVIDER container is the only custom-built image; all other tools use public biocontainers images.

---

## Input

Provide a CSV samplesheet with three columns:

```
sample_id,fastq,metadata_json
P001,/absolute/path/to/P001.fastq.gz,
P002,/absolute/path/to/P002.fastq.gz,{"collection_date":"2026-01-15"}
```

Rules: `sample_id` must match `^[A-Za-z0-9._-]+$`, must be unique, FASTQ must be readable.
`metadata_json` is optional; leave the field empty or omit the column.

Full schema details: [docs/configuration.md](docs/configuration.md).

---

## Profiles

Container engine and executor are separate concerns — combine them as needed:

| Profile | Role | Typical use |
|---|---|---|
| `docker` | Docker engine, local executor | Mac development (resource caps for M-series) |
| `docker_mac` | Docker engine, local executor | Apple Silicon Docker with serial LoFreq calling |
| `conda` | conda/mamba, local executor | Any platform without Docker/Singularity |
| `singularity` | Singularity engine only | Combine with an HPC profile (see below) |
| `katana` | SLURM executor, Katana resources | UNSW Katana HPC |
| `gadi` | SLURM executor, Gadi resources | NCI Gadi HPC |
| `test` | Bundled test data | CI and smoke testing |

```bash
# Mac with Docker
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results

# Apple Silicon Docker, with serial LoFreq calling
nextflow run main.nf -profile docker_mac --input samplesheet.csv --outdir results

# Mac or Linux with conda
nextflow run main.nf -profile conda --input samplesheet.csv --outdir results

# Katana HPC
nextflow run main.nf -profile singularity,katana --input samplesheet.csv --outdir results

# Gadi HPC
nextflow run main.nf -profile singularity,gadi --gadi_project <code> --input samplesheet.csv --outdir results
```

### Katana: specifying your account

The `katana` profile defaults to `--account=oz000`. Override it per run on the command line:

```bash
nextflow run main.nf -profile singularity,katana \
    -process.clusterOptions='--account=<your_account>' \
    --input samplesheet.csv --outdir results
```

`-process.clusterOptions` is Nextflow's config-scope CLI override syntax (note the single dash, not `--`). It takes the highest precedence and requires no config file.

For repeated use, add a one-line file (e.g. `~/.nextflow/config` or a local `my_account.config`):

```groovy
process.clusterOptions = '--account=<your_account>'
```

then pass it with `-c my_account.config`. A `~/.nextflow/config` is loaded automatically on every run without any `-c` flag.

Invocation examples and per-profile notes: [docs/usage.md](docs/usage.md).

---

## Key parameters

| Parameter | Default | Purpose |
|---|---|---|
| `--input` | required | Samplesheet CSV |
| `--reference_panel` | `assets/hcv_references.fasta` | 238-sequence HCV reference panel |
| `--host_reference` | null | GRCh38 FASTA or .mmi; only required with `--use_minimap2` or `--use_hostile` |
| `--outdir` | `results` | Output directory |
| `--min_secondary_fraction` | `0.05` | Read fraction threshold for mixed-infection detection |
| `--lofreq_max_depth` | `5000` | rasusa depth cap for LoFreq |
| `--lofreq_pp_threads` | `8` | LoFreq parallel workers; `1` uses serial `lofreq call` |
| `--devider_max_depth` | `5000` | rasusa depth cap for DEVIDER |
| `--run_clair3` | `false` | Enable optional Clair3 corroboration |

Full parameter reference with descriptions and tuning guidance: [docs/parameters.md](docs/parameters.md).

---

## Outputs

```
results/
  <sample_id>/
    qc/
    genotyping/
    consensus/<GT>/
    mapping/<GT>/
    variants/<GT>/
    haplotypes/<GT>/
    reports/
  pipeline_info/
  reports/
```

Description of every output file: [docs/output.md](docs/output.md).

---

## Documentation

| Document | Contents |
|---|---|
| [docs/usage.md](docs/usage.md) | Installation, profiles, running, troubleshooting |
| [docs/parameters.md](docs/parameters.md) | Full parameter reference |
| [docs/important_considerations.md](docs/important_considerations.md) | DEVIDER haplotype recovery: read length, SNP density, error rates, tuning |
| [docs/output.md](docs/output.md) | Output file descriptions |
| [docs/configuration.md](docs/configuration.md) | Config defaults, samplesheet schema, compute profiles |
| [docs/data_flow.md](docs/data_flow.md) | Step-by-step data flow with exact commands |
| [docs/testing.md](docs/testing.md) | Testing and validation strategy |
| [docs/limitations.md](docs/limitations.md) | Known limitations and future work |

---

## Citation

If you use this pipeline, please cite the underlying tools — in particular LoFreq, DEVIDER, minimap2, and Nextflow.
