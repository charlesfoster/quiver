# hcv-quasi

Nextflow DSL2 pipeline for reproducible HCV quasispecies analysis from ONT PromethION reads: genome-wide low-frequency variant calling (LoFreq, ≥1% AF) and global haplotype reconstruction (DEVIDER), with automatic detection and per-genotype branching for mixed-genotype infections.

Architecture, design decisions, and tool version rationale: [CLAUDE.md](CLAUDE.md).

---

## Quick start

```bash
# 1. Clone the repository
git clone https://github.com/charlesfoster/hcv-quasi.git
cd hcv-quasi

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
| Docker | any recent | Required for `local` profile |
| Singularity / Apptainer | any recent | Required for `katana` and `gadi` profiles |
| conda / micromamba | micromamba ≥1.5 | Required for `conda_local` profile |

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

| Profile | Executor | Containers | Typical use |
|---|---|---|---|
| `local` | local | Docker (Rosetta on Apple Silicon) | Development and small runs on a Mac |
| `conda_local` | local | conda/micromamba | Apple Silicon Mac without Docker, or Docker-free environments |
| `katana` | SLURM | Singularity | UNSW Katana HPC |
| `gadi` | SLURM | Singularity | NCI Gadi HPC (requires `--gadi_project`) |
| `test` | local | Docker | CI and smoke testing |

Invocation examples and per-profile notes: [docs/usage.md](docs/usage.md).

---

## Key parameters

| Parameter | Default | Purpose |
|---|---|---|
| `--input` | required | Samplesheet CSV |
| `--reference_panel` | `assets/hcv_references.fasta` | 238-sequence HCV reference panel |
| `--host_reference` | required | GRCh38 FASTA or .mmi for host depletion |
| `--outdir` | `results` | Output directory |
| `--min_secondary_fraction` | `0.05` | Read fraction threshold for mixed-infection detection |
| `--lofreq_max_depth` | `5000` | rasusa depth cap for LoFreq |
| `--devider_max_depth` | `1000` | rasusa depth cap for DEVIDER |
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
| [CLAUDE.md](CLAUDE.md) | Architecture, design decisions, tool inventory |
| [docs/usage.md](docs/usage.md) | Installation, profiles, running, troubleshooting |
| [docs/parameters.md](docs/parameters.md) | Full parameter reference |
| [docs/output.md](docs/output.md) | Output file descriptions |
| [docs/configuration.md](docs/configuration.md) | Config defaults, samplesheet schema, compute profiles |
| [docs/data_flow.md](docs/data_flow.md) | Step-by-step data flow with exact commands |
| [docs/testing.md](docs/testing.md) | Testing and validation strategy |
| [docs/limitations.md](docs/limitations.md) | Known limitations and future work |

---

## Citation

If you use this pipeline, please cite the underlying tools listed in [CLAUDE.md](CLAUDE.md) Section 4, in particular LoFreq, DEVIDER, minimap2, and Nextflow.
