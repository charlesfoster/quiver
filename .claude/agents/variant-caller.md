---
name: variant-caller
description: Implement Round 2 mapping, LoFreq preprocessing (indelqual + alnqual), coverage QC with mosdepth, rasusa subsampling for LoFreq, LoFreq variant calling, variant filtering, and optional Clair3 corroboration. Covers Prompts 12–17 from docs/implementation_prompts.md.
model: claude-sonnet-4-6
---

You are implementing **Prompts 12–17** of the `hcv-quasi` pipeline: Round 2 mapping through variant filtering.

## Context
Read before starting:
- `docs/data_flow.md` — Steps 5.11–5.16 (exact commands)
- `docs/architecture_reasoning.md` — Section 7 (LoFreq ONT preprocessing) and Section 8 (depth normalisation)
- `docs/implementation_prompts.md` — Prompts 12–17 gotchas
- `CLAUDE.md` — Section 3 D9 (LoFreq decision), Section 4 (container tags)

## Files to create

### `modules/local/minimap2_round2.nf`
Map per-genotype reads to the sample-specific consensus (single reference).

```bash
minimap2 \
    -ax map-ont \
    -t ${task.cpus} \
    -Y --MD --eqx \
    -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ONT" \
    ${consensus_mmi} \
    ${reads} \
  | samtools sort -@ ${task.cpus} -O bam -o ${meta.id}.${meta.genotype}.round2.bam -
samtools index -@ ${task.cpus} ${meta.id}.${meta.genotype}.round2.bam

# Emit LOW_MAPPING_RATE sentinel if <90% primary mapped
pct=$(samtools flagstat ${meta.id}.${meta.genotype}.round2.bam \
  | awk '/primary mapped/ {mapped=$1} /primary$/ {total=$1} END {printf "%.0f", 100*mapped/total}')
if [ "$pct" -lt 90 ]; then touch ${meta.id}.${meta.genotype}.LOW_MAPPING_RATE; fi
```

**Key difference from Round 1:** No `-N 5 --secondary=no` — mapping against a single reference.
Label: `process_high`

### `modules/local/mosdepth.nf`
Run on the full-depth Round 2 BAM (NOT subsampled).

```bash
mosdepth -t ${task.cpus} -n --fast-mode --by 100 \
    ${meta.id}.${meta.genotype}.r2 \
    ${bam}

# Parse mean coverage from summary and emit sentinel if too low
mean_cov=$(grep "total_region" ${meta.id}.${meta.genotype}.r2.mosdepth.summary.txt \
           | awk '{print $4}')
if python3 -c "import sys; sys.exit(0 if float('$mean_cov') < ${params.min_mean_coverage} else 1)"; then
    touch ${meta.id}.${meta.genotype}.LOW_COVERAGE
fi
```

Outputs: summary txt, region BED (gzipped), optional `LOW_COVERAGE` sentinel.
Label: `process_medium`

### `modules/local/rasusa.nf`
Generic subsampling module (used for both LoFreq and DEVIDER caps).

```bash
rasusa reads \
    --coverage ${coverage_cap} \
    --genome-size ${params.genome_size ?: 9646} \
    --seed ${seed} \
    -o ${meta.id}.${meta.genotype}.${suffix}.fastq.gz \
    ${reads}
```

Parameterised: `coverage_cap`, `seed`, `suffix` passed via `task.ext.args` or ext.args2. If input coverage is already below the cap, rasusa passes reads through unchanged.
Label: `process_medium`

### `subworkflows/local/prep_lofreq_input.nf`
1. Call `RASUSA` with `coverage_cap = params.lofreq_max_depth`, `seed = params.rasusa_seed_lofreq`
2. Re-map the subsampled reads to consensus using `MINIMAP2_ROUND2` (same module, different inputs)
3. Run `LOFREQ_PREPROCESS` on the subsampled BAM
4. Emit the LoFreq-ready BAM alongside the full-depth BAM (both paths needed downstream)

### `modules/local/lofreq_preprocess.nf`
```bash
# Step 1: indel quality scoring (required for indel calling)
lofreq indelqual --dindel \
    -f ${consensus} \
    -o ${meta.id}.${meta.genotype}.iq.bam \
    ${bam}
samtools index ${meta.id}.${meta.genotype}.iq.bam

# Step 2: alignment quality recalibration (optional, graceful fallback)
lofreq alnqual -b \
    ${meta.id}.${meta.genotype}.iq.bam \
    ${consensus} \
    > ${meta.id}.${meta.genotype}.iq.alnq.bam \
  || {
    echo "WARNING: lofreq alnqual failed, using indelqual-only BAM" >&2
    cp ${meta.id}.${meta.genotype}.iq.bam ${meta.id}.${meta.genotype}.iq.alnq.bam
  }
samtools index ${meta.id}.${meta.genotype}.iq.alnq.bam
```

Output BAM name: `${meta.id}.${meta.genotype}.iq.alnq.bam`
Label: `process_medium`

### `modules/local/lofreq_call.nf`
```bash
lofreq call-parallel \
    --pp-threads ${task.cpus} \
    --call-indels \
    --min-mq ${params.min_mq} \
    --min-bq ${params.min_bq} \
    --min-cov 20 \
    --sig ${params.lofreq_sig} \
    -f ${consensus} \
    -o ${meta.id}.${meta.genotype}.lofreq.vcf \
    ${bam}

bgzip ${meta.id}.${meta.genotype}.lofreq.vcf
tabix -p vcf ${meta.id}.${meta.genotype}.lofreq.vcf.gz
```

Note: `lofreq call-parallel` uses `--pp-threads` for its thread count (not `--threads`). Bind to `task.cpus`.
Label: `process_high`

### `modules/local/variant_filter.nf`
```bash
bcftools view \
    -i "INFO/AF >= ${params.min_report_af} & INFO/DP >= ${params.min_variant_depth}" \
    -Oz -o ${meta.id}.${meta.genotype}.lofreq.filtered.vcf.gz \
    ${vcf}
tabix -p vcf ${meta.id}.${meta.genotype}.lofreq.filtered.vcf.gz

# TSV for reporting
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/DP\t%INFO/SB\n' \
    ${meta.id}.${meta.genotype}.lofreq.filtered.vcf.gz \
    > ${meta.id}.${meta.genotype}.variants.tsv
```

Use `-i` (inclusion filter), not `-e`. LoFreq INFO fields: `AF`, `DP`, `SB`, `DP4`, `INDEL`.
Label: `process_low`

### `modules/local/clair3.nf` (opt-in)
Only invoked when `params.run_clair3 = true`.

```bash
run_clair3.sh \
    --bam_fn=${bam} \
    --ref_fn=${consensus} \
    --threads=${task.cpus} \
    --platform="ont" \
    --model_path="${params.clair3_model_dir ?: '/opt/models/r1041_e82_400bps_hac_v520'}" \
    --output=${meta.id}.${meta.genotype}.clair3 \
    --haploid_sensitive
```

Then run `bin/clair3_concordance.py` to compare Clair3 calls at AF ≥ 0.25 against LoFreq filtered calls. Emit `clair3_concordance.tsv`.

Model must be `r1041_e82_400bps_hac_v520` (R10.4.1 HAC v5.2.0). Accept model path override via `params.clair3_model_dir`.
Label: `process_high`

### `bin/clair3_concordance.py`
Compare two VCFs at positions where LoFreq AF ≥ 0.25. For each such variant, check if Clair3 also calls it. Output TSV: `chrom, pos, ref, alt, lofreq_af, clair3_called (bool), clair3_gt`.

## Success criteria
- Round 2 BAM has read group; is sorted and indexed
- `lofreq indelqual` output BAM has `BI` and `BD` tags (check with `samtools view -H`)
- LoFreq VCF is bgzip-compressed and tabix-indexed
- Filtered VCF contains only variants with AF ≥ `params.min_report_af`
- `bcftools stats` runs without error on both VCFs
- All AF values in VCF are in (0, 1]
