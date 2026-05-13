---
name: genotype-classifier
description: Design and implement the genotype classification script and Nextflow module — walks the Round 1 BAM, assigns reads to genotypes, detects mixed infections, computes ambiguous pools, and emits the genotype_summary.json. Covers Prompt 8 from docs/implementation_prompts.md. This task requires careful reasoning about the classification algorithm.
model: claude-opus-4-7
---

You are implementing **Prompt 8** of the `hcv-quasi` pipeline: genotype classification with mixed-infection detection.

## Context
Read before starting:
- `docs/data_flow.md` — Step 5.8 (full script logic description)
- `docs/architecture_reasoning.md` — Section 5 (mixed genotype detection rationale)
- `docs/configuration.md` — `genotype_summary.json` schema (mandatory — your output must conform exactly)
- `docs/implementation_prompts.md` — Prompt 8 gotchas
- `CLAUDE.md` — Section 3 D6, D7

## Files to create

### `bin/classify_genotype.py`
Python 3.11. Uses `pysam`.

**Invocation:**
```
python classify_genotype.py \
    --bam round1.bam \
    --panel-fasta panel.fasta \
    --sample-id P001 \
    --min-secondary-fraction 0.05 \
    --ambiguous-delta-as 20 \
    --out-tsv genotype_assignments.tsv \
    --out-json genotype_summary.json
```

#### Step 1 — Parse subtype/genotype from FASTA headers

The FASTA header prefix encodes both genotype and subtype. Examples from the HCV panel:
- `1a_M62321.1` → genotype `1`, subtype `1a`
- `2b_D10988.1` → genotype `2`, subtype `2b`
- `6xj_EF589068.1` → genotype `6`, subtype `6xj` (multi-letter subtype)
- `1_AJ238799.1` → genotype `1`, subtype `1` (no letter — pure genotype assignment)

Regex: `^([0-9]+)([a-z]*)_` where group 1 is the genotype number and groups 1+2 together form the subtype.

Build a dict `{ref_id: (genotype, subtype)}` by scanning all FASTA headers in `panel.fasta` with this regex. Any header that does not match should be logged as a warning and skipped.

#### Step 2 — Walk primary alignments

Use `pysam.AlignmentFile(bam, 'rb')`. For each read:
- Skip if `read.is_secondary` or `read.is_supplementary` or `read.is_unmapped`
- The primary alignment's `reference_name` is the best-hit reference (minimap2 with `--secondary=no` ensures this)
- Extract `AS` tag: `read.get_tag('AS')` (alignment score)
- Extract `XS` tag if present: `read.get_tag('XS')` (second-best alignment score; minimap2 emits this when `-N > 1`)

#### Step 3 — Ambiguous pool detection

A read is **ambiguous** if: `XS` tag is present AND `(AS - XS) < ambiguous_delta_as`.

Ambiguous reads are counted separately and excluded from the genotype fraction computation. They are NOT assigned to any genotype. Rationale: reads at the boundary between two subtypes (e.g., 1a vs 1b) have similar scores to both; assigning them to either inflates one genotype's count.

#### Step 4 — Compute genotype fractions

- `total_mapped` = all non-ambiguous primary aligned reads
- For each genotype `g`: `reads_g` = count of non-ambiguous reads whose best-hit reference is in genotype `g`
- `fraction_g` = `reads_g / total_mapped`
- `top_subtype` for each genotype = the subtype within that genotype with the most reads
- `top_reference` for each genotype = the specific reference (ref_id) within that genotype's top subtype with the most reads

#### Step 5 — Mixed infection detection

Sort genotypes by `reads` descending. The genotype with the most reads is `primary_genotype`.

`is_mixed = True` if any non-primary genotype has `fraction >= min_secondary_fraction`.

`secondary_genotypes` = list of non-primary genotype labels where `fraction >= min_secondary_fraction`, sorted by fraction descending.

`branches_to_run` = `[primary_genotype] + secondary_genotypes`

#### Step 6 — TSV output

One row per primary mapped read (excluding unmapped):
```
read_id  ref_id  subtype  genotype  AS  XS  ambiguous
```
`ambiguous` column: `True` if the read is in the ambiguous pool, else `False`.

#### Step 7 — JSON output

Must conform exactly to the schema in `docs/configuration.md`. Validate with `jsonschema` if available; otherwise manually check all required keys before writing.

### Edge cases to handle explicitly

1. **All reads ambiguous:** `total_mapped = 0` → set `is_mixed = false`, `primary_genotype = null`, `branches_to_run = []`. This is unusual and should be logged as a warning.

2. **Single genotype, no XS tags:** minimap2 with `-N 5` should produce XS tags, but if absent (e.g., only one reference in panel matches), skip ambiguous detection gracefully. Log a note.

3. **Genotype prefix not in panel:** if a reference in the panel has a non-standard header format, emit a warning and exclude those reads from classification but include them in the total read count for `total_mapped_reads`.

4. **Subtype-mixed note:** if within the dominant genotype, any non-dominant subtype has ≥ 20% of that genotype's reads, add an informational `subtype_mixed_note` field to the JSON. This does NOT trigger a branch — it is purely informational for the report.

5. **Tie in primary genotype:** if two genotypes have identical read counts, assign primary genotype arbitrarily (e.g., numerically lowest) and log a warning.

### `modules/local/genotype_classify.nf`

```groovy
process GENOTYPE_CLASSIFY {
    label 'process_medium'
    container '...'  // python:3.11-slim with pysam

    input:
    tuple val(meta), path(bam), path(bai)
    path panel_fasta

    output:
    tuple val(meta), path("${meta.id}.genotype_summary.json"), emit: summary
    tuple val(meta), path("${meta.id}.genotype_assignments.tsv"), emit: assignments

    script:
    """
    python ${projectDir}/bin/classify_genotype.py \\
        --bam ${bam} \\
        --panel-fasta ${panel_fasta} \\
        --sample-id ${meta.id} \\
        --min-secondary-fraction ${params.min_secondary_fraction} \\
        --ambiguous-delta-as ${params.ambiguous_delta_as} \\
        --out-tsv ${meta.id}.genotype_assignments.tsv \\
        --out-json ${meta.id}.genotype_summary.json
    """
}
```

Container must include `pysam`. Use the lofreq biocontainer (which ships pysam), or add a separate `python:3.11-slim` container with pysam installed.

## Success criteria
- JSON output for a single-genotype sample: `is_mixed = false`, `branches_to_run = ["1"]` (or correct genotype)
- JSON output for a mixed sample (70% 1a + 30% 3a): `is_mixed = true`, `secondary_genotypes = ["3"]`
- TSV has one row per primary mapped read, no unmapped reads
- `sum(fraction for all genotypes) + ambiguous_fraction ≈ 1.0` (within float rounding)
- Genotype label `"1"` not `"1a"` — genotype is the numeric prefix only; subtype is the full prefix
- The script exits non-zero and prints a clear error if the BAM has no primary alignments at all
