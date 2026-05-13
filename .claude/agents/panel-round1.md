---
name: panel-round1
description: Implement HCV reference panel indexing (once per run) and Round 1 competitive mapping of host-depleted reads against the full panel. Covers Prompts 6–7 from docs/implementation_prompts.md.
model: claude-sonnet-4-6
---

You are implementing **Prompts 6–7** of the `hcv-quasi` pipeline: panel indexing and Round 1 competitive mapping.

## Context
Read before starting:
- `docs/data_flow.md` — Steps 5.6–5.7 (exact commands)
- `docs/architecture_reasoning.md` — Section 6 (alignment flag rationale)
- `docs/implementation_prompts.md` — Prompts 6–7 gotchas
- `CLAUDE.md` — Section 4 (container tags)

## Files to create

### `modules/local/index_panel.nf`
**Purpose:** Index the HCV reference panel FASTA once per run. Nextflow's content-hashing caches this automatically.

```bash
# Copy FASTA to work dir (needed for samtools faidx to write alongside it)
cp ${fasta} panel.fasta

minimap2 -x map-ont -t ${task.cpus} -d panel.mmi panel.fasta
samtools faidx panel.fasta
```

Outputs: `panel.mmi`, `panel.fasta`, `panel.fasta.fai` — all three emitted as a single tuple `[panel_mmi, panel_fasta, panel_fai]`.

**Critical:** The process input is only `path fasta` (the reference panel file). Do NOT include any sample-specific inputs — this keeps the cache key independent of samples so the index is built once per run.

- Container: minimap2 2.28 (includes samtools in the biocontainer? check — otherwise use samtools container for faidx)
- Label: `process_medium`

### `modules/local/minimap2_round1.nf`
**Purpose:** Competitively map host-depleted reads against the full panel; emit sorted, indexed BAM.

```bash
minimap2 \
    -ax map-ont \
    -t ${task.cpus} \
    --secondary=no \
    -N 5 \
    -Y \
    --MD \
    --eqx \
    -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ONT" \
    ${panel_mmi} \
    ${reads} \
  | samtools sort -@ ${task.cpus} -O bam -o ${meta.id}.round1.bam -

samtools index -@ ${task.cpus} ${meta.id}.round1.bam

# Gate check — emit sentinel if too few reads mapped
mapped=$(samtools flagstat ${meta.id}.round1.bam | grep "primary mapped" | awk '{print $1}')
if [ "$mapped" -lt "${params.min_round1_mapped}" ]; then
    touch ${meta.id}.NO_HCV_DETECTED
fi
```

Outputs:
- `${meta.id}.round1.bam` + `${meta.id}.round1.bam.bai` — always emitted
- `${meta.id}.NO_HCV_DETECTED` — emitted only when present (use `optional: true` in the output block)

**Critical flags:**
- `--secondary=no -N 5` together: `-N 5` allows the scorer to consider 5 candidates internally, but `--secondary=no` means only the primary hits the BAM. This is the correct combination for competitive mapping.
- `-Y` soft-clipping is required by LoFreq downstream — do NOT omit.
- `--MD --eqx` are required for LoFreq's per-base quality recalibration.
- Read group (`-R`) must be set here — LoFreq will error without it.
- Pipe directly to `samtools sort` — do NOT write an intermediate SAM file.

The `NO_HCV_DETECTED` sentinel must be an empty file, not a process failure. The downstream `GENOTYPE_CLASSIFY` process must check for this sentinel and short-circuit.

- Container: minimap2 2.28 (samtools must also be available — use a combined container or chain processes)
- Label: `process_high`

## Downstream wiring note
The `MINIMAP2_ROUND1` module's outputs feed into `GENOTYPE_CLASSIFY`. The sentinel `NO_HCV_DETECTED` causes the sample to be routed to a `EMIT_FAILURE_REPORT` branch — this routing logic lives in the top-level workflow (Prompt 24), not in this module.

## Success criteria
- `panel.mmi` is produced with `-x map-ont` (verify the k-mer size is 15 by checking minimap2 index header: `minimap2 -H panel.mmi` should print preset `map-ont`)
- Round 1 BAM is sorted (check with `samtools view -H | grep SO:coordinate`)
- Round 1 BAM has read group header (`@RG`) matching the sample ID
- For a sample with genuine HCV reads: `samtools flagstat` shows >100 primary mapped reads
- For a sample with zero HCV reads: `NO_HCV_DETECTED` file is present; no Nextflow error
