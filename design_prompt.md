---

You are an expert bioinformatician specialising in viral genomics, long-read sequencing, and reproducible pipeline development. Your task is to **reason carefully through the architecture** of a production-grade bioinformatics pipeline, then **write a complete `CLAUDE.md` file** that a separate Claude Sonnet 4.6 instance will use as its sole reference to implement the pipeline from scratch.

---

## Scientific Context

**Organism:** Hepatitis C virus (HCV)
**Goal:** Infer the quasispecies present in each patient sample — specifically:
1. Genome-wide low-frequency variant calling
2. Global haplotype reconstruction where possible

**Sequencing platform:** Oxford Nanopore Technologies (ONT), PromethION
**Basecalling model:** `dna_r10.4.1_e8.2_400bps_hac@v5.2.0`

**Representative read statistics (highest coverage sample):**

| Metric | Value |
|---|---|
| num_seqs | 911,758 |
| sum_len | ~2.84 Gb |
| min_len | 63 bp |
| avg_len | 3,119.8 bp |
| max_len | 143,252 bp |
| Q1 / Q2 / Q3 | 842 / 1,782 / 4,483 bp |
| N50 | 6,130 bp |
| AvgQual | 15.2 |
| Q20 reads | 80% |

**HCV genome size:** ~9,600 bp. Reads substantially longer than this are chimeric artefacts and must be removed early.

**Reference panel:** A curated FASTA containing representative sequences across all known HCV subtypes (hcv_references.fasta).

**Samples:** Multiple patient samples. Some may have >100,000× depth against the HCV genome after host depletion.

---

## Constraints and Preferences

**Variant calling:** LoFreq is the preferred variant caller for its sensitivity to low-frequency variants. You may propose an alternative only if you can make a strong, evidence-based case that it handles ONT error profiles and low-frequency viral variants better. If you retain LoFreq, reason through what preprocessing is required to make it perform optimally on ONT data (e.g., indel quality score insertion via `lofreq indelqual`).

**Haplotype reconstruction:** DEVIDER is the preferred tool for global haplotype reconstruction from long reads. You may propose an alternative only if you can make a strong case for it. Reason through whether global haplotype reconstruction is tractable across the full ~9.6 kb genome or whether a windowed/regional approach is more appropriate given ONT error rates and typical quasispecies complexity in HCV.

**Mixed genotype infections:** The pipeline must detect and handle the possibility (rare but real) that a sample contains reads from more than one HCV genotype. Reason through how and where in the pipeline this should be detected and what should happen downstream.

**Host depletion:** Samples are derived from clinical material and will contain human reads. Removal must be efficient given the read volumes involved. Choose the most computationally appropriate approach for this data scale and the available compute environments.

**Depth normalisation:** Some samples exceed 100,000× coverage against the HCV genome. Reason through whether this is problematic for any tools in the pipeline (particularly LoFreq and DEVIDER), and if so, incorporate `rasusa` or an equivalent for depth-capped subsampling prior to affected steps. Preserve the full-depth BAM for other steps where depth is beneficial.

**Pipeline approach:** Reason explicitly through the following candidate strategies and recommend one (or a hybrid):
- *De novo assembly → mapping → variant calling → quasispecies inference*
- *Competitive mapping to reference panel → variant calling → quasispecies inference*
- *Competitive mapping → variant calling → per-sample consensus → second-round mapping → quasispecies inference* (two-round approach)
- Any other approach you consider superior

Key considerations for this decision: reference panel diversity vs. patient strain divergence, ONT error rate implications for assembly quality, whether a per-sample consensus improves variant calling accuracy enough to justify the added complexity.

**Workflow manager:** Choose between **Nextflow (DSL2)** and **Snakemake** and justify your choice. Requirements: full reproducibility, containerisation (Docker + Singularity/Apptainer for HPC compatibility), resumability, per-sample parallelism.

**Compute environments:**
- *Primary:* MacBook Pro, Apple M5 Max, 64 GB RAM (ARM64, macOS). No GPU requirement.
- *HPC:* UNSW Katana cluster (SLURM-based)
- *HPC:* NCI Gadi (national facility, SLURM-based)
The pipeline must run correctly on all three with environment-specific profiles.

**Reproducibility:** All tools must be version-pinned. Containers must be specified per process/rule. A `conda`/`mamba` environment file should also be provided as a fallback.

---

## Reasoning Instructions

Before writing any output, work through the following questions explicitly in a `<reasoning>` block. Do not skip steps:

1. What is the best overall pipeline strategy for this data type and scientific goal, and why?
2. Does a two-round mapping approach meaningfully improve quasispecies inference accuracy for ONT HCV data, and is the added complexity justified?
3. At what read length should the chimeric read filter be set, and why?
4. What is the most efficient host depletion strategy at this data scale?
5. How should mixed genotype infections be detected, and what should the pipeline do when one is found?
6. What alignment tool and parameters are optimal for ONT reads against a diverse HCV reference panel?
7. What preprocessing is required before LoFreq to make it perform well on ONT data?
8. Is depth normalisation required before LoFreq? Before DEVIDER? What depth ceiling is appropriate?
9. Is global haplotype reconstruction tractable across the full HCV genome with DEVIDER given ONT error rates, or should it be windowed?
10. What QC metrics should be computed and at what stages?
11. Nextflow or Snakemake — which is better suited to this pipeline and these compute environments?
12. What are the key failure modes or edge cases the pipeline must handle gracefully (e.g., very low viral load samples, single-genotype vs. mixed, samples with poor coverage over specific genome regions)?

---

## Output Instructions

After your reasoning block, produce a single, complete `CLAUDE.md` file. This file will be placed in the root of the pipeline repository and will be the **only context** available to Claude Sonnet 4.6 when it implements the pipeline. Write it accordingly — it must be self-contained and unambiguous.

The `CLAUDE.md` must contain the following sections:

### 1. Project Overview
Brief description of the pipeline's purpose, inputs, outputs, and scientific goals.

### 2. Architecture Diagram (ASCII or Mermaid)
A visual representation of the full pipeline DAG, including branching logic for mixed infection detection.

### 3. Design Decisions
For each major design decision (pipeline strategy, tool choices, workflow manager, depth normalisation, etc.), state: the decision, the alternatives considered, and the rationale. This is the record Sonnet must not deviate from.

### 4. Tool Inventory
A table of every tool used, its version, its role in the pipeline, and its container/conda source.

### 5. Data Flow Specification
For each pipeline step, specify exactly:
- Step name
- Input(s): file type, naming convention, expected format
- Output(s): file type, naming convention
- Tool + exact command template with all parameters
- Resource requirements (CPUs, memory, time estimate per sample)
- Pass/fail criteria (what constitutes a failed step)

### 6. Configuration
Document all user-configurable parameters (reference panel path, depth caps, minimum coverage thresholds, frequency cutoffs for variant calling, etc.) and their default values.

### 7. Compute Profiles
Specify the Nextflow/Snakemake profile configuration for: `local` (Mac M5 Max), `katana` (UNSW HPC), `gadi` (NCI HPC).

### 8. Implementation Prompts for Sonnet
A numbered list of self-contained implementation tasks for Sonnet 4.6 to execute in order. Each prompt must be specific enough that Sonnet can implement it without architectural judgement. Include: what to build, what file(s) to create, what the expected output is, and any gotchas or constraints to observe. Tasks should be scoped to be completable in a single focused session.

### 9. Testing and Validation Strategy
How to verify each module works correctly. Include: a small synthetic test dataset specification, expected outputs for each step, and integration test criteria.

### 10. Known Limitations and Future Work
Honest documentation of what this pipeline does not handle, and what would be needed to extend it.

---

Write the `CLAUDE.md` section completely and without truncation. Do not summarise or use placeholders — Sonnet will have no other reference.

---