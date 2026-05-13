---
name: containers-and-docs
description: Finalise the DEVIDER container CI workflow and write all user-facing documentation (README, usage guide, output guide, parameters reference). Covers Prompts 26–27 from docs/implementation_prompts.md.
model: claude-sonnet-4-6
---

You are implementing **Prompts 26–27** of the `hcv-quasi` pipeline: container CI and documentation.

## Context
Read before starting:
- `docs/implementation_prompts.md` — Prompts 26–27
- `CLAUDE.md` — all sections (this is the architectural reference that docs must link to)
- `docs/` — all sub-documents (docs must link to these, not duplicate them)

## Prompt 26 — Container CI

The `containers/devider/Dockerfile` is created by the `devider-pipeline` agent. Your job here is to add a GitHub Actions workflow for building and validating it.

### `.github/workflows/build-containers.yml`
```yaml
name: Build containers

on:
  push:
    paths:
      - 'containers/**'
  pull_request:
    paths:
      - 'containers/**'

jobs:
  build-devider:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Set up QEMU (for multi-arch)
        uses: docker/setup-qemu-action@v3

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v3

      - name: Build DEVIDER image (amd64)
        uses: docker/build-push-action@v5
        with:
          context: containers/devider
          platforms: linux/amd64
          push: false
          tags: hcv-quasi/devider:v0.0.1
          cache-from: type=gha
          cache-to: type=gha,mode=max

      - name: Smoke test DEVIDER binary
        run: |
          docker run --platform linux/amd64 hcv-quasi/devider:v0.0.1 --help
```

This workflow runs on any change to `containers/`. It does not push to a registry (no credentials needed); it only validates the build succeeds and the binary is executable.

Note to user: when ready to publish, add `docker/login-action` and set `push: true` with a registry target. Pin base images in `Dockerfile` by digest before production use.

## Prompt 27 — Documentation

All docs must **link to** CLAUDE.md and `docs/` sub-documents rather than duplicate content from them.

### `README.md`
~30 lines maximum. Cover:
1. One-sentence description of the pipeline
2. Quick-start (3 commands: clone, run help, run on test data)
3. Table of contents pointing to `docs/` files
4. Link to `CLAUDE.md` for architecture/design

### `docs/usage.md`
How to run the pipeline. Sections:
1. **Requirements** — Nextflow ≥24.10, Docker or Singularity or conda (micromamba), DEVIDER container or cargo
2. **Profiles** — one paragraph each on `local`, `katana`, `gadi`, `test`; exact invocation for each
3. **Inputs** — samplesheet format, FASTQ requirements, reference panel note
4. **Smoke test commands** — the exact commands from `docs/testing.md` Section 9.4
5. **Resuming failed runs** — `nextflow run ... -resume`
6. **Resource customisation** — how to override `params.max_cpus` etc.

### `docs/output.md`
What the pipeline produces. For each output directory:
- Path template (e.g. `results/${sample_id}/variants/${genotype}/`)
- Files produced
- Format and what it contains
- When it may be absent (e.g. `merged_haplotypes.fasta` absent if DEVIDER failed)

### `docs/parameters.md`
Full parameter reference table. One row per parameter:

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Path to samplesheet CSV |
| `--reference_panel` | `assets/hcv_references.fasta` | HCV reference panel FASTA |
| ... | ... | ... |

Use the full parameter list from `docs/configuration.md`. Add a "When to change" column briefly noting non-obvious tuning scenarios (e.g. for `--devider_max_depth`: "increase for higher sensitivity in low-diversity samples; decrease if DEVIDER runs out of memory").

## Style rules for all docs
- Do NOT duplicate content from `CLAUDE.md` or `docs/` sub-documents — link to them
- No emojis
- Code blocks for all commands
- Keep sentences short
- Assume the reader has Nextflow experience but not necessarily HCV biology knowledge
