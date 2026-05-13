---
name: genotype-brancher
description: Implement the per-genotype read partitioning script and the Nextflow branching subworkflow that routes reads, dominant references, and metadata to per-genotype downstream paths. Covers Prompts 9–10 from docs/implementation_prompts.md. Requires careful reasoning about Nextflow channel semantics and mixed vs. single-genotype routing.
model: claude-opus-4-7
---

You are implementing **Prompts 9–10** of the `hcv-quasi` pipeline: read partitioning and the genotype branching subworkflow.

## Context
Read before starting:
- `docs/data_flow.md` — Steps 5.9 and the routing note in Step 5.8
- `docs/architecture_reasoning.md` — Section 5 (mixed detection), Section 12 (failure modes for NO_HCV/LOW_COVERAGE)
- `docs/implementation_prompts.md` — Prompts 9 and 10 gotchas
- `CLAUDE.md` — Section 3 D12 (why per-genotype branching)
- `docs/configuration.md` — `genotype_summary.json` schema (you must parse this)

## Prompt 9 — Per-genotype read partitioning

### `bin/partition_reads.py`
Python 3.11 + pysam. Only invoked when `is_mixed = true`.

**Invocation:**
```
python partition_reads.py \
    --bam round1.bam \
    --assignments genotype_assignments.tsv \
    --genotypes 1,3 \
    --sample-id P001 \
    --outdir .
```

**Algorithm:**
1. Load the TSV produced by `classify_genotype.py` into a dict: `{read_id: (genotype, ambiguous)}`
2. Open the Round 1 BAM with pysam.
3. Walk primary alignments. For each read:
   - Look up the read's genotype from the dict
   - If `ambiguous = True`: write to `{sample_id}.ambiguous.fastq.gz`
   - Otherwise: write to `{sample_id}.{genotype}.reads.fastq.gz`
4. Emit the read's sequence from `read.query_sequence` (BAM SEQ field) — do NOT reverse-complement. The SEQ field in a BAM is always the sequence as it was input to the aligner; minimap2 stores it in original orientation regardless of the mapping strand. Leave strand handling to minimap2 on re-mapping.
5. Emit quality scores from `read.query_qualities` converted to Phred+33 ASCII.
6. FASTQ record format: `@{read.query_name}\n{seq}\n+\n{qual}`

**Output files:**
- `{sample_id}.{genotype}.reads.fastq.gz` — one per genotype in `--genotypes`
- `{sample_id}.ambiguous.fastq.gz` — ambiguous reads (QC only, not used downstream)

Compress with `gzip` module or pipe through `pigz` if available.

**Edge case:** A read_id in the BAM may not appear in the TSV (e.g. supplementary alignments that were skipped in the classifier). Skip these reads silently.

### `modules/local/partition_reads.nf`
Only invoked when `is_mixed = true`.

```groovy
process PARTITION_READS {
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), path(assignments_tsv)
    val genotypes_list  // e.g. ["1", "3"]

    output:
    tuple val(meta), path("${meta.id}.*.reads.fastq.gz"), emit: partitioned
    path "${meta.id}.ambiguous.fastq.gz", emit: ambiguous

    script:
    def gt_arg = genotypes_list.join(',')
    """
    python ${projectDir}/bin/partition_reads.py \\
        --bam ${bam} \\
        --assignments ${assignments_tsv} \\
        --genotypes ${gt_arg} \\
        --sample-id ${meta.id} \\
        --outdir .
    """
}
```

The output `partitioned` channel emits a list of FASTQ paths — one per genotype. Downstream, this must be `flatMap`-ed and combined with the genotype label.

## Prompt 10 — Genotype branching subworkflow

### `subworkflows/local/genotype_branch.nf`
This is the routing hub of the pipeline. Its job: given a `genotype_summary.json` and the host-depleted FASTQ, emit one channel item per downstream branch.

**Input channels:**
1. `[meta, genotype_summary_json]` — from `GENOTYPE_CLASSIFY`
2. `[meta, hostdep_fastq]` — from `HOST_DEPLETE`
3. `[meta, round1_bam, round1_bai, assignments_tsv]` — from `MINIMAP2_ROUND1` + `GENOTYPE_CLASSIFY`
4. `[panel_fasta, panel_mmi, panel_fai]` — singleton from `INDEX_PANEL`

**What it must emit:**
A channel of tuples: `[branch_meta, reads_fastq, dom_ref_fasta]` where:
- `branch_meta = meta + [genotype: "1a"]` (copy of `meta` with `genotype` field added)
- `reads_fastq` = per-genotype FASTQ (partitioned, or full hostdep for single-genotype)
- `dom_ref_fasta` = the dominant panel reference FASTA for this genotype (extracted from `panel_fasta` by ref_id from the JSON)

**Routing logic:**

```groovy
workflow GENOTYPE_BRANCH {
    take:
    ch_summary    // [meta, summary_json]
    ch_hostdep    // [meta, fastq]
    ch_round1_bam // [meta, bam, bai, assignments_tsv]
    ch_panel      // [panel_fasta, panel_mmi, panel_fai]

    main:
    // 1. Parse summary JSON to determine routing
    ch_parsed = ch_summary.map { meta, json ->
        def summary = new groovy.json.JsonSlurper().parse(json)
        [meta, summary]
    }

    // 2. Route: NO_HCV → failure branch; LOW_COVERAGE → flag; otherwise → process
    ch_parsed.branch {
        no_hcv:       it[1].branches_to_run.isEmpty()
        processable:  true
    }.set { ch_routed }

    // 3. For processable samples: emit one item per branch
    ch_branches = ch_routed.processable.flatMap { meta, summary ->
        summary.branches_to_run.collect { gt ->
            def branch_meta = meta + [genotype: gt]
            [branch_meta, summary]
        }
    }

    // 4. For mixed samples: invoke PARTITION_READS
    // For single-genotype: pass hostdep FASTQ directly

    // 5. Extract dominant reference per branch from panel FASTA
    // (seqkit or samtools faidx using the top_reference from the JSON)

    emit:
    branches = ch_out   // [branch_meta, reads_fastq, dom_ref_fasta]
    no_hcv   = ch_routed.no_hcv.map { meta, _ -> meta }
}
```

**Dominant reference extraction:**
For each branch, use `samtools faidx panel.fasta "${top_reference}"` to extract the single reference sequence as a temporary FASTA. This small process (`EXTRACT_DOM_REF`) wraps:
```bash
samtools faidx ${panel_fasta} "${ref_id}" > ${meta.id}.${meta.genotype}.dom_ref.fasta
```

**Mixed vs. single-genotype routing:**
- If `is_mixed = false`: bypass `PARTITION_READS`; emit `[branch_meta, hostdep_fastq, dom_ref_fasta]`
- If `is_mixed = true`: call `PARTITION_READS` first; then for each output FASTQ emit `[branch_meta, per_gt_fastq, dom_ref_fasta]`

The routing must work correctly for both cases using Nextflow's `join`, `combine`, or `branch` operators. Use `meta.id` as the join key for combining channels.

**The `no_hcv` emit:**
Samples in the `no_hcv` branch are passed to a `EMIT_FAILURE_REPORT` stub (to be wired in Prompt 24). They do not enter any downstream processing.

## Reasoning guidance
Think carefully about Nextflow channel semantics:
- `flatMap` is the right operator for expanding one sample into multiple branches
- The join between `ch_branches` (which now has `branch_meta` including genotype) and the per-genotype FASTQs from `PARTITION_READS` requires that the key includes the genotype label — use `branch_meta.id + "_" + branch_meta.genotype` or a structured key `[branch_meta.id, branch_meta.genotype]`
- Avoid using `.collect()` at this stage — it would block until all samples are processed, preventing streaming parallelism
- The panel FASTA is a singleton; use `.combine(ch_branches)` to pair it with every branch item

## Success criteria
- For a single-genotype sample: exactly one tuple emitted, with `branch_meta.genotype` set correctly
- For a 70%/30% mixed sample: exactly two tuples emitted (one per genotype)
- `NO_HCV_DETECTED` samples do not enter downstream processing
- Each emitted FASTQ contains only reads assigned to that genotype (verified by mapping back to the panel)
- `dom_ref_fasta` contains exactly one FASTA record matching the `top_reference` field from the JSON
