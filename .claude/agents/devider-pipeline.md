---
name: devider-pipeline
description: Implement the DEVIDER haplotype reconstruction pipeline — rasusa subsampling for DEVIDER depth cap, the DEVIDER Nextflow module, the DEVIDER Dockerfile (multi-stage Rust build), and the cargo install script. Covers Prompts 19–20 from docs/implementation_prompts.md.
model: claude-sonnet-4-6
---

You are implementing **Prompts 19–20** of the `hcv-quasi` pipeline: DEVIDER subsampling, the DEVIDER module, its container, and its install script.

## Context
Read before starting:
- `docs/data_flow.md` — Steps 5.18–5.19 (exact DEVIDER command and flag notes)
- `docs/architecture_reasoning.md` — Section 9 (haplotype reconstruction strategy)
- `docs/implementation_prompts.md` — Prompts 19–20 gotchas (read every bullet carefully)
- `CLAUDE.md` — Section 4 (DEVIDER version note), Section 3 D11

## Files to create

### `subworkflows/local/prep_devider_input.nf`
1. Call `RASUSA` with `coverage_cap = params.devider_max_depth`, `seed = params.rasusa_seed_devider`, `suffix = "devider"`
2. Re-map subsampled reads to consensus:
   ```bash
   minimap2 -ax map-ont -t ${task.cpus} -Y --MD --eqx \
       ${consensus_mmi} ${reads} \
     | samtools sort -@ ${task.cpus} -O bam -o devider.bam -
   samtools index devider.bam
   ```
3. Emit `[meta, devider_bam, devider_bam_bai, consensus_fasta, consensus_mmi]`

### `modules/local/devider.nf`
**Input:** `[meta, devider_bam, devider_bam_bai, consensus_fasta, filtered_vcf, filtered_vcf_tbi]`

**Command (exact — verified against DEVIDER v0.0.1 `--help`):**
```bash
devider \
    -b ${devider_bam} \
    -r ${consensus} \
    -v ${filtered_vcf} \
    -o devider_out \
    -O \
    -t ${task.cpus} \
    --preset nanopore-r10 \
    --min-cov ${params.devider_min_cov} \
    --min-abund ${params.devider_min_abund} \
    --output-reads \
    --allele-output
```

**Failure handling — this is critical:**
```bash
devider ... || {
    echo "DEVIDER exited non-zero for ${meta.id} ${meta.genotype}" >&2
    mkdir -p devider_out
    touch devider_out/devider.failed
}
```

Always emit the output directory — whether DEVIDER succeeded or failed. The `devider.failed` marker signals downstream that no haplotypes are available without failing the sample.

**Output:** `path "devider_out"` (the directory). Downstream consumers must `ls devider_out/` and parse what is present.

- Container: built from `containers/devider/Dockerfile` (see below)
- Label: `process_high_memory` (worst case: 64 GB, 4 hours)

**DEVIDER flag notes — do not change without verifying `devider --help`:**
- `--preset nanopore-r10` — correct for R10.4.1 chemistry. `nanopore-r9` is wrong (too lenient). There is **no** preset called `ont`.
- `-O` — overwrites output directory (required for Nextflow work-dir re-runs)
- `--output-reads` — emits haplotype-tagged BAM needed by the stitching agent
- `--allele-output` — writes nucleotide alleles (A/C/G/T), not 0/1 codes
- v0.0.1 has **no `--merge-windows` flag** — do not add it
- DEVIDER requires the filtered VCF as input; it phases over those SNPs but does **not** call variants itself

### `containers/devider/Dockerfile`
Multi-stage build: compile DEVIDER from source (Rust), copy only the binary to a slim runtime image.

```dockerfile
# Stage 1: build
FROM rust:1.78-slim AS builder

RUN apt-get update && apt-get install -y \
    git cmake pkg-config libssl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN cargo install \
    --git https://github.com/bluenote-1577/devider \
    --tag v0.0.1 \
    --root /usr/local

# Stage 2: runtime
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y \
    libssl3 ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/bin/devider /usr/local/bin/devider

ENTRYPOINT ["/usr/local/bin/devider"]
```

Build target: `linux/amd64`. The M5 Max Mac runs this via Rosetta — acceptable for development.

Add a `containers/devider/build.sh` helper:
```bash
#!/usr/bin/env bash
set -euo pipefail
docker build --platform linux/amd64 -t hcv-quasi/devider:v0.0.1 containers/devider/
echo "Built hcv-quasi/devider:v0.0.1"
```

### `bin/install_devider.sh`
For use in the conda/no-Docker profile. Installs DEVIDER from source via cargo.

```bash
#!/usr/bin/env bash
set -euo pipefail

DEVIDER_VERSION="v0.0.1"
INSTALL_DIR="${1:-$HOME/.local/bin}"

if command -v devider &>/dev/null; then
    installed=$(devider --version 2>&1 | head -1)
    echo "devider already installed: $installed"
    exit 0
fi

if ! command -v cargo &>/dev/null; then
    echo "ERROR: cargo not found. Install Rust via https://rustup.rs/" >&2
    exit 1
fi

echo "Installing DEVIDER $DEVIDER_VERSION..."
cargo install \
    --git https://github.com/bluenote-1577/devider \
    --tag $DEVIDER_VERSION \
    --root "$INSTALL_DIR/.."

echo "DEVIDER installed at $(command -v devider)"
devider --version
```

## DEVIDER version note
DEVIDER v0.0.1 is the only tagged release as of 2026-05-13. At implementation time:
1. Check `https://github.com/bluenote-1577/devider/releases` for newer tags
2. If a newer tag exists, test the CLI flags against `devider --help` before updating the pin
3. The flag interface can change between point releases — do not assume the v0.0.1 invocation above works unchanged

## Success criteria
- `docker build --platform linux/amd64 -t hcv-quasi/devider:v0.0.1 containers/devider/` succeeds
- `docker run hcv-quasi/devider:v0.0.1 --help` prints DEVIDER usage
- DEVIDER module runs on test data without process failure (even if DEVIDER exits non-zero — the module itself must not fail)
- `devider.failed` marker is present when DEVIDER fails; absent when it succeeds
- Output directory contains haplotype FASTA files when DEVIDER succeeds
