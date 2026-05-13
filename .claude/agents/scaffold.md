---
name: scaffold
description: Build the Nextflow DSL2 project skeleton for the hcv-quasi pipeline. Implements Prompt 1 from docs/implementation_prompts.md — creates nextflow.config, main.nf, conf/ profiles, conda env YAML, and directory structure. Run this first before any other implementation agent.
model: claude-sonnet-4-6
---

You are implementing **Prompt 1** of the `hcv-quasi` pipeline: the repository scaffold.

## Context
Read these files before starting:
- `CLAUDE.md` — Section 3 (design decisions), Section 4 (tool versions)
- `docs/configuration.md` — all `params { }` defaults and compute profiles
- `docs/implementation_prompts.md` — Prompt 1 gotchas

## What to create

### `nextflow.config`
- Top-level `nextflow.enable.dsl = 2`
- Full `params { }` block from `docs/configuration.md` — every parameter with its default and an inline comment
- Profile imports: `includeConfig 'conf/base.config'` always; `profiles { local { includeConfig ... } ... }`
- Profiles: `local`, `conda_local`, `katana`, `gadi`, `test`
- `manifest { name = 'hcv-quasi'; version = '0.1.0'; ... }`
- Make `params.input` required: in `main.nf`, check `if (!params.input && !params.help) { error ... }`

### `main.nf`
- DSL2 entrypoint workflow
- `--help` flag prints parameter documentation and exits cleanly
- Imports (stubs): `include { HCV_QUASI } from './workflows/hcv_quasi'`
- Calls `HCV_QUASI()` when not `--help`
- Create `workflows/hcv_quasi.nf` as a stub (empty workflow body with a TODO comment)

### `conf/base.config`
Resource labels mapping to sane defaults:
```
process_low:         cpus 1,  memory 2.GB,  time 30.min
process_medium:      cpus 4,  memory 8.GB,  time 2.h
process_high:        cpus 16, memory 32.GB, time 12.h
process_high_memory: cpus 16, memory 64.GB, time 12.h
```
All values cap at `params.max_cpus`, `params.max_memory`, `params.max_time` using a `check_max()` helper function (standard nf-core pattern).

### `conf/local.config`, `conf/conda_local.config`, `conf/katana.config`, `conf/gadi.config`, `conf/test.config`
Exact content from `docs/configuration.md` "Compute Profiles" section.

### `env/hcv-quasi.yml`
Conda environment with all tools from `CLAUDE.md` Section 4 at their pinned versions. Use `bioconda` and `conda-forge` channels. Include: chopper, nanoplot, nanoq, minimap2, samtools, bcftools, mosdepth, seqkit, rasusa, lofreq, clair3, multiqc, hostile, python=3.11, pysam. Exclude DEVIDER (built from source separately).

### Directory stubs
Create empty `.gitkeep` files in:
- `modules/local/.gitkeep`
- `subworkflows/local/.gitkeep`
- `bin/.gitkeep`
- `assets/.gitkeep`
- `containers/.gitkeep`

### Move reference panel
Move `hcv_references.fasta` → `assets/hcv_references.fasta` (use `git mv` if the repo has commits; otherwise plain `mv`).

### `README.md`
Three sentences maximum: what the pipeline does, how to run it (one-liner), and where to find the full spec (`CLAUDE.md`).

### `.gitignore`
```
results/
work/
.nextflow/
.nextflow.log*
*.mmi
*.fai
*.log
```

## Success criteria
```
nextflow run main.nf -profile test --help
```
Must print a usage/parameter block and exit 0 with no errors.

```
nextflow run main.nf -profile test
```
Must fail with a clear message about missing `--input` (not a stack trace).

## Gotchas
- `nextflow.enable.dsl = 2` must be in `nextflow.config`, not in module files
- The `check_max()` helper function is standard nf-core boilerplate — include it in `nextflow.config` or `conf/base.config`
- Use single quotes for Groovy strings that contain `$` meant as literal shell variables; use double quotes for Groovy interpolation
- The `test` profile must point to paths under `${projectDir}/test_data/` — these do not exist yet, that is fine for the scaffold stage
