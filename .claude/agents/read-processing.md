---
name: read-processing
description: Implement the read intake and QC subworkflow — samplesheet validation, raw NanoPlot/nanoq QC, chopper length/quality filtering, and minimap2 host depletion. Covers Prompts 2–5 from docs/implementation_prompts.md.
model: claude-sonnet-4-6
---

You are implementing **Prompts 2–5** of the `hcv-quasi` pipeline: everything from samplesheet parsing through host-depleted reads.

## Context
Read before starting:
- `docs/data_flow.md` — Steps 5.1–5.5 (exact commands and resource specs)
- `docs/implementation_prompts.md` — Prompts 2–5 gotchas
- `CLAUDE.md` — Section 4 (container tags and conda packages)
- `conf/base.config` — resource label definitions

## Files to create

### `bin/check_samplesheet.py`
Python 3.11 script. Reads CSV from `sys.argv[1]`.
- Validates columns: `sample_id`, `fastq` are required; `metadata_json` is optional
- Rejects: duplicate `sample_id`; `sample_id` matching `[^\w.\-]` (whitespace, `/`, `\`); FASTQ path that doesn't exist or isn't readable; empty FASTQ (0 bytes)
- On success: prints a table summary of validated samples to stderr; writes validated tuples as JSON to stdout (consumed by the Nextflow process)
- On failure: prints a clear error message to stderr and exits non-zero; the first violation found stops processing (fail-fast)
- Schema: `[{"id": "P001", "fastq": "/abs/path.fastq.gz", "metadata": {...}}]`

### `modules/local/input_check.nf`
Process wrapping `bin/check_samplesheet.py`. Emits channel of `[meta, file(fastq)]` tuples.
- `meta` is a Groovy map: `[id: row.id, metadata: row.metadata]`
- Container: `python:3.11-slim`; conda: `conda-forge::python=3.11`
- Label: `process_low`

### `modules/local/nanoplot.nf`
- Command: `NanoPlot --fastq ${reads} -o ${meta.id}_nanoplot --tsv_stats --no_static --threads ${task.cpus}`
- Output: directory of NanoPlot HTML + TSVs, published to `${params.outdir}/${meta.id}/qc/raw/`
- Container/conda from `CLAUDE.md` Section 4 (NanoPlot 1.43.0)
- Label: `process_medium`
- Do NOT fail if NanoPlot warns about low read count (common in test data)

### `modules/local/nanoq.nf`
- Command: `nanoq -i ${reads} --json -o ${meta.id}.nanoq.json --report`
- Output: JSON stats file
- Container/conda from `CLAUDE.md` Section 4 (nanoq 0.10.0)
- Label: `process_low`

### `subworkflows/local/raw_qc.nf`
Includes and calls `NANOPLOT` and `NANOQ` in parallel on the same input channel. Emits combined QC outputs for MultiQC collection.

### `modules/local/chopper.nf`
- Command (exact):
  ```bash
  zcat ${reads} | chopper \
      -q ${params.min_qual} \
      --minlength ${params.min_length} \
      --maxlength ${params.max_length} \
      --threads ${task.cpus} \
      2> chopper.log \
    | pigz -p ${task.cpus} > ${meta.id}.filtered.fastq.gz
  ```
- Output: filtered FASTQ + `chopper.log`
- If output FASTQ is empty (0 reads): emit a sentinel file `${meta.id}.ALL_READS_FILTERED` alongside the empty FASTQ; do NOT fail the process
- Container/conda: chopper 0.9.2
- Label: `process_medium`

### `modules/local/host_deplete.nf`
Two process variants controlled by `params.use_hostile`:

**Default (minimap2):**
```bash
minimap2 -ax map-ont -t ${task.cpus} ${host_index} ${reads} \
  | samtools view -@ ${task.cpus} -b -f 4 - \
  | samtools fastq -@ ${task.cpus} - \
  | pigz -p ${task.cpus} > ${meta.id}.hostdep.fastq.gz

# Compute stats
total=$(zcat ${reads} | awk 'NR%4==1' | wc -l)
kept=$(zcat ${meta.id}.hostdep.fastq.gz | awk 'NR%4==1' | wc -l)
python3 -c "
import json, sys
total, kept = int(sys.argv[1]), int(sys.argv[2])
host = total - kept
print(json.dumps({'total': total, 'kept': kept, 'host_removed': host,
                  'host_fraction': round(host/total, 4) if total else 0}))
" $total $kept > ${meta.id}.host_stats.json
```
- Input: filtered FASTQ + pre-built host index (`.mmi` file from a separate INDEX_HOST process)
- Container: minimap2 2.28 + samtools 1.21 (use the samtools biocontainer which includes both, or chain them)
- Label: `process_high_memory`

**Alternative (hostile):** a second process `HOST_DEPLETE_HOSTILE` using `hostile clean` — activated when `params.use_hostile = true`. Same inputs/outputs.

### `bin/get_host_reference.sh`
Shell script: if `params.host_reference` is a URL or the file doesn't exist, download GRCh38 no-alt from NCBI (`https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`) and index it with minimap2. Accept the output path as `$1`. Skip download if the `.mmi` already exists.

## Sentinel handling
The `INPUT_CHECK` module must emit a sentinel `${meta.id}.EMPTY_INPUT` (an empty file) for samples whose FASTQ is empty — downstream processes gate on its absence.

The `CHOPPER` module emits `${meta.id}.ALL_READS_FILTERED` for samples where all reads are filtered out.

Both sentinels are collected and passed to a simple `EMIT_SAMPLE_FLAG` process that writes them to `${params.outdir}/${meta.id}/qc/flags.txt` for the reporter.

## Success criteria
- `check_samplesheet.py` rejects: duplicate IDs, paths with spaces (if unquoted in CSV), missing files
- A FASTQ with reads spanning lengths 100–15000 bp, after chopper with defaults, retains only reads 200–10000 bp
- Host-depleted FASTQ contains no reads that map to `mini_host.fasta`
