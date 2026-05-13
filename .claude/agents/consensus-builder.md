---
name: consensus-builder
description: Implement the per-genotype consensus building subworkflow — re-map reads to the dominant panel reference, call high-confidence majority variants, mask low-coverage positions, and apply bcftools consensus. Covers Prompt 11 from docs/implementation_prompts.md.
model: claude-sonnet-4-6
---

You are implementing **Prompt 11** of the `hcv-quasi` pipeline: per-genotype sample-specific consensus building.

## Context
Read before starting:
- `docs/data_flow.md` — Step 5.10 (all sub-steps with exact commands)
- `docs/architecture_reasoning.md` — Section 2 (two-round mapping rationale) and Section 1 (why not medaka)
- `docs/implementation_prompts.md` — Prompt 11 gotchas
- `CLAUDE.md` — Section 4 (container tags); Section 3 D8 (consensus decision)

## Inputs
This subworkflow receives a channel of tuples: `[meta, reads_fastq, dominant_ref_fasta]`
- `meta` carries `id` (sample ID) and `genotype` (e.g. `"1a"`)
- `reads_fastq` is the per-genotype reads (or all hostdep reads for single-genotype samples)
- `dominant_ref_fasta` is the single FASTA record for the dominant panel reference for this genotype

## Files to create

### `modules/local/minimap2_dom_ref.nf`
Re-map per-genotype reads to the dominant panel reference (not the full panel).

```bash
minimap2 -ax map-ont -t ${task.cpus} -Y --MD --eqx \
    -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ONT" \
    ${dom_ref} ${reads} \
  | samtools sort -@ ${task.cpus} -O bam -o ${meta.id}.${meta.genotype}.dom.bam -
samtools index ${meta.id}.${meta.genotype}.dom.bam
```

Label: `process_high`

### `modules/local/consensus_variants.nf`
Call high-confidence variants for consensus construction. This is NOT the final variant calling step — it is only used to build the consensus reference.

```bash
bcftools mpileup \
    -f ${dom_ref} \
    -d 8000 \
    -Q 7 \
    -q 20 \
    -a AD,DP \
    ${bam} \
  | bcftools call --ploidy 1 -mv -Oz -o ${meta.id}.${meta.genotype}.consensus.vcf.gz

bcftools index ${meta.id}.${meta.genotype}.consensus.vcf.gz

# Filter: keep only majority-allele variants (DP >= 10, AD[1]/DP >= 0.5)
bcftools view -i 'INFO/DP >= 10 && (INFO/AD[1] / INFO/DP) >= 0.5' \
    -Oz -o ${meta.id}.${meta.genotype}.consensus.filt.vcf.gz \
    ${meta.id}.${meta.genotype}.consensus.vcf.gz
bcftools index ${meta.id}.${meta.genotype}.consensus.filt.vcf.gz
```

Note: `--ploidy 1` is correct — HCV is a haploid RNA virus.

Label: `process_medium`

### `modules/local/build_mask_bed.nf`
Generate a BED file of positions to mask (coverage < `params.min_consensus_cov`).

```bash
mosdepth --quantize 0:1:${params.min_consensus_cov}: \
    --no-per-base \
    ${meta.id}.${meta.genotype}.cov \
    ${bam}

# Extract positions in the 0:1 and 1:{threshold} quantize bins → these are the low-coverage regions
zcat ${meta.id}.${meta.genotype}.cov.quantized.bed.gz \
  | awk -v thresh=${params.min_consensus_cov} '$4 ~ /^0:/ || ($4 ~ /^[0-9]/ && int($4) < thresh)' \
  > mask.bed

# If no low-coverage regions, create an empty mask.bed
touch mask.bed
```

The quantize bins produce labels like `0:1`, `1:10`, `10:500` — extract only bins below the threshold.

Label: `process_low`

### `modules/local/apply_consensus.nf`
Apply variants and masking to produce the final sample-specific consensus FASTA.

```bash
bcftools consensus \
    -f ${dom_ref} \
    -m ${mask_bed} \
    -o ${meta.id}.${meta.genotype}.consensus.fasta \
    ${consensus_vcf}

# Rename header to stable identifier
sed -i "s/>.*/>${meta.id}_${meta.genotype}_consensus/" \
    ${meta.id}.${meta.genotype}.consensus.fasta

# Index and build minimap2 index
samtools faidx ${meta.id}.${meta.genotype}.consensus.fasta
minimap2 -x map-ont -d ${meta.id}.${meta.genotype}.consensus.mmi \
    ${meta.id}.${meta.genotype}.consensus.fasta
```

Outputs: `consensus.fasta`, `.fai`, `.mmi` emitted as tuple `[meta, consensus_fasta, consensus_fai, consensus_mmi]`

Pass/fail check:
```bash
length=$(grep -v "^>" ${meta.id}.${meta.genotype}.consensus.fasta | tr -d '\\n' | wc -c)
n_count=$(grep -v "^>" ${meta.id}.${meta.genotype}.consensus.fasta | tr -d '\\n' | tr -cd 'Nn' | wc -c)
n_frac=$(python3 -c "print($n_count / $length if $length > 0 else 0)")

if [ "$length" -lt 4000 ] || [ "$length" -gt 11000 ]; then
    echo "WARNING: consensus length $length out of expected range [4000, 11000]" >&2
fi
if python3 -c "import sys; sys.exit(0 if $n_frac >= 0.30 else 1)"; then
    touch ${meta.id}.${meta.genotype}.LOW_COVERAGE_CONSENSUS
fi
```

Label: `process_medium`

### `subworkflows/local/build_consensus.nf`
Wire the four modules above in sequence. Input: `[meta, reads, dom_ref_fasta]`. Output: `[meta, consensus_fasta, consensus_fai, consensus_mmi]` plus optional `LOW_COVERAGE_CONSENSUS` sentinel.

## Success criteria
- Consensus FASTA has exactly one record with header `${sample_id}_${genotype}_consensus`
- Length is in [4000, 11000] bp for normal test samples
- `.fai` and `.mmi` files are present alongside the FASTA
- Masked positions (coverage < 10) appear as `N` in the consensus
- `LOW_COVERAGE_CONSENSUS` sentinel is emitted when N-fraction ≥ 30%

## Gotchas
- `bcftools call --ploidy 1` — HCV is haploid. Using ploidy 2 would produce diploid genotypes and corrupt the consensus.
- The mask BED logic: `mosdepth --quantize` bins are open on the right, so `0:10` means [0, 10). Use the `0:1` (zero coverage) and `1:{threshold}` bins.
- After `sed -i` on macOS: `sed -i ''` (BSD sed). In Linux containers: `sed -i` (GNU sed). For portability: `sed -i.bak "..." && rm *.bak` or use Python.
- The `.mmi` index is used downstream by Round 2 mapping — it must be built with `-x map-ont`, same preset as all other minimap2 steps.
