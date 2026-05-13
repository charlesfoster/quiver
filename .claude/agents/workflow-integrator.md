---
name: workflow-integrator
description: Wire all subworkflows and modules into the top-level hcv_quasi.nf workflow, implement correct channel routing through the genotype branching expand/collapse, and ensure end-to-end execution on test data. Covers Prompt 24 from docs/implementation_prompts.md. Requires careful reasoning about Nextflow channel topology and the per-genotype fan-out/fan-in pattern.
model: claude-opus-4-7
---

You are implementing **Prompt 24** of the `hcv-quasi` pipeline: top-level workflow wiring.

## Context
Read before starting:
- `docs/data_flow.md` — all steps (you are wiring them together)
- `docs/architecture_reasoning.md` — Section 12 (failure modes — every sentinel must route correctly)
- `docs/implementation_prompts.md` — Prompt 24 gotchas
- `CLAUDE.md` — the ASCII pipeline diagram in Section 2 (your wiring must match this)
- All other agent definitions — understand what each module/subworkflow emits and consumes

## Prerequisite check
Before writing code, verify that all modules and subworkflows referenced exist:
- `modules/local/input_check.nf`
- `subworkflows/local/raw_qc.nf`
- `modules/local/chopper.nf`
- `modules/local/host_deplete.nf`
- `modules/local/index_panel.nf`
- `modules/local/minimap2_round1.nf`
- `modules/local/genotype_classify.nf`
- `subworkflows/local/genotype_branch.nf`
- `subworkflows/local/build_consensus.nf`
- `modules/local/minimap2_round2.nf`
- `modules/local/mosdepth.nf`
- `subworkflows/local/prep_lofreq_input.nf`
- `modules/local/lofreq_call.nf`
- `modules/local/variant_filter.nf`
- `subworkflows/local/prep_devider_input.nf`
- `modules/local/devider.nf`
- `modules/local/stitch_haplotypes.nf`
- `modules/local/sample_report.nf`
- `modules/local/multiqc.nf`

If any are missing, list them as blocking dependencies and stop.

## File to create

### `workflows/hcv_quasi.nf`

This is the primary workflow. High-level structure:

```groovy
include { INPUT_CHECK         } from '../modules/local/input_check'
include { RAW_QC              } from '../subworkflows/local/raw_qc'
include { CHOPPER             } from '../modules/local/chopper'
include { HOST_DEPLETE        } from '../modules/local/host_deplete'
include { INDEX_PANEL         } from '../modules/local/index_panel'
include { MINIMAP2_ROUND1     } from '../modules/local/minimap2_round1'
include { GENOTYPE_CLASSIFY   } from '../modules/local/genotype_classify'
include { GENOTYPE_BRANCH     } from '../subworkflows/local/genotype_branch'
include { BUILD_CONSENSUS     } from '../subworkflows/local/build_consensus'
include { MINIMAP2_ROUND2     } from '../modules/local/minimap2_round2'
include { MOSDEPTH            } from '../modules/local/mosdepth'
include { PREP_LOFREQ_INPUT   } from '../subworkflows/local/prep_lofreq_input'
include { LOFREQ_CALL         } from '../modules/local/lofreq_call'
include { VARIANT_FILTER      } from '../modules/local/variant_filter'
include { PREP_DEVIDER_INPUT  } from '../subworkflows/local/prep_devider_input'
include { DEVIDER_RUN         } from '../modules/local/devider'
include { STITCH_HAPLOTYPES   } from '../modules/local/stitch_haplotypes'
include { SAMPLE_REPORT       } from '../modules/local/sample_report'
include { MULTIQC             } from '../modules/local/multiqc'
```

### Channel topology — the hardest part

The pipeline has a **fan-out** (one sample → N genotype branches) and a **fan-in** (N branches → one per-sample report). You must get this right.

#### Fan-out point
`GENOTYPE_BRANCH` emits a channel where each item carries `branch_meta = [id: "P001", genotype: "1a"]`. From this point forward, every downstream process receives one item per branch (e.g., P001_1a, P001_3a for a mixed sample).

The key: `branch_meta.id` is the SAMPLE key; `branch_meta.genotype` is the BRANCH key. Downstream processes use `${meta.id}.${meta.genotype}` for output naming.

#### Fan-in point
At `SAMPLE_REPORT`, you need all branches for a sample collected together. Use:
```groovy
// Collect all per-branch outputs by sample ID
ch_report_inputs = ch_variants.join(ch_haplotypes, by: [0])  // join on branch_meta
    .map { branch_meta, variants, haplotypes -> 
        [branch_meta.id, branch_meta.genotype, variants, haplotypes] 
    }
    .groupTuple(by: 0)  // group all branches by sample ID
    .map { sample_id, genotypes, variants_list, haplotypes_list ->
        def report_meta = [id: sample_id]
        [report_meta, genotypes, variants_list, haplotypes_list]
    }
```

This pattern — `map` to extract the sample ID key, then `groupTuple(by: 0)` — is the standard approach for collapsing per-branch outputs back to per-sample.

#### Sentinel handling
After `MINIMAP2_ROUND1`, samples with `NO_HCV_DETECTED` must not enter `GENOTYPE_CLASSIFY` or anything downstream. Gate using:
```groovy
ch_round1.branch {
    no_hcv:       it[2] != null  // third element is the sentinel path (optional emit)
    processable:  true
}.set { ch_r1_routed }

// no_hcv branch → emit failure flag to results, then stop
ch_r1_routed.no_hcv.map { meta, bam, bai, sentinel ->
    log.warn "Sample ${meta.id}: NO_HCV_DETECTED — skipping all downstream processing"
    [meta, sentinel]
} | EMIT_FAILURE_FLAG

// processable branch continues normally
ch_r1_routed.processable.map { meta, bam, bai, _ -> [meta, bam, bai] }
    | GENOTYPE_CLASSIFY
```

Similarly gate on:
- `ALL_READS_FILTERED` after `CHOPPER`
- `LOW_COVERAGE` after `MOSDEPTH` — skip DEVIDER but NOT LoFreq
- `LOW_COVERAGE_CONSENSUS` after `BUILD_CONSENSUS` — downstream continues but flagged
- `devider.failed` marker from `DEVIDER_RUN` — `STITCH_HAPLOTYPES` handles this gracefully

#### Collecting QC for MultiQC
MultiQC needs a flat list of all QC directories from all samples and all branches. Collect with:
```groovy
ch_multiqc_inputs = Channel.empty()
    .mix(ch_nanoplot_dirs)
    .mix(ch_nanoq_jsons)
    .mix(ch_mosdepth_dirs)
    .mix(ch_lofreq_stats)
    .collect()
```

#### Optional Clair3
Gate on `params.run_clair3`:
```groovy
if (params.run_clair3) {
    CLAIR3_CORROBORATE(ch_round2_bam, ch_consensus)
    ch_clair3_out = CLAIR3_CORROBORATE.out
} else {
    ch_clair3_out = Channel.empty()
}
```

### `EMIT_FAILURE_FLAG` stub
A minimal process that writes a failure flag to the output directory. Needed for samples that hit `NO_HCV_DETECTED`, `ALL_READS_FILTERED`, etc.

```groovy
process EMIT_FAILURE_FLAG {
    label 'process_low'
    publishDir "${params.outdir}/${meta.id}/", mode: 'copy'

    input:
    tuple val(meta), path(flag_file)

    output:
    tuple val(meta), path("${meta.id}.PIPELINE_FLAG.txt")

    script:
    """
    echo "Pipeline status flag: \$(basename ${flag_file})" > ${meta.id}.PIPELINE_FLAG.txt
    echo "Sample: ${meta.id}" >> ${meta.id}.PIPELINE_FLAG.txt
    echo "Timestamp: \$(date -u)" >> ${meta.id}.PIPELINE_FLAG.txt
    """
}
```

### `main.nf` update
Update to call `HCV_QUASI()`:
```groovy
include { HCV_QUASI } from './workflows/hcv_quasi'

workflow {
    if (params.help) {
        // print usage
        exit 0
    }
    if (!params.input) {
        error "ERROR: --input is required. Run with --help for usage."
    }
    HCV_QUASI()
}
```

## Reasoning guidance

Think through these before writing any code:

1. **The `groupTuple` timing problem:** `groupTuple` waits for ALL items with the same key to arrive before emitting. For a mixed sample with 2 branches, this works fine — both branches will complete and the tuple will be emitted. But if you use `groupTuple` too early (e.g., on the raw sample channel), it might deadlock. Use it only at the fan-in point.

2. **Join key consistency:** The `branch_meta` map must have the same structure throughout. A map `[id: "P001", genotype: "1a"]` is equal to another `[id: "P001", genotype: "1a"]` in Groovy (maps compare by value). Use this as the join key when combining outputs from parallel processes.

3. **The optional emit pattern:** When a process may or may not emit a sentinel file (e.g., `NO_HCV_DETECTED`), use `optional: true` in the output block. The channel item for that process will include `null` for the sentinel path when it's absent, and a file path when it's present. This lets you `branch` on null-vs-non-null.

4. **publishDir and fan-out:** With fan-out, multiple processes will publish to `${params.outdir}/${meta.id}/variants/1a/` and `${params.outdir}/${meta.id}/variants/3a/`. This works because `meta.genotype` differentiates the paths. Verify the publishDir templates in each module use `${meta.id}/${meta.genotype}` correctly.

## Success criteria
- `nextflow run main.nf -profile test` completes end-to-end on both test samples
- `NO_HCV_DETECTED` sample emits a flag file and does not cause any downstream errors
- Mixed sample produces two sets of outputs (one per genotype)
- Single-genotype sample produces one set of outputs
- MultiQC collects data from all samples
- Per-sample HTML reports are generated for all samples that reach the reporting stage
- `-resume` works correctly: re-running with the same inputs does not re-execute cached processes
