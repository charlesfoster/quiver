---
name: test-data
description: Generate synthetic test data and build the integration test harness for the hcv-quasi pipeline. Covers Prompt 25 from docs/implementation_prompts.md. Run after all pipeline modules are implemented.
model: claude-sonnet-4-6
---

You are implementing **Prompt 25** of the `hcv-quasi` pipeline: synthetic test data and the integration test harness.

## Context
Read before starting:
- `docs/implementation_prompts.md` — Prompt 25 (all details)
- `docs/testing.md` — Section 9.2 (expected outcomes and tolerances)
- `docs/configuration.md` — test profile parameters
- `assets/hcv_references.fasta` — source of reference sequences for simulation

## Files to create

### `test_data/mini_panel.fasta`
Extract three reference sequences from `assets/hcv_references.fasta`:
- One genotype 1a reference (header prefix `1a_`)
- One genotype 2b reference (header prefix `2b_`)
- One genotype 3a reference (header prefix `3a_`)

Use `seqkit grep -r -p "^1a_" assets/hcv_references.fasta | seqkit head -n 1` (and similarly for 2b, 3a). Concatenate into `test_data/mini_panel.fasta`.

### `test_data/mini_host.fasta`
A small human sequence for host depletion testing. Extract 10,000 bp from a known human sequence that will NOT cross-map to HCV. Use any 10 kb window from a simple repeat-free region of human chr1 or use a synthetic sequence. Label it `>chr1_mini`.

### `test_data/single_gt.fastq.gz`
Simulated ONT reads from one 1a reference, with ~1% host contamination.

Use BadRead:
```bash
# HCV reads (~500 reads, ~50× coverage of a 9.6 kb genome)
badread simulate \
    --reference test_data/mini_panel.fasta \
    --quantity 500x \
    --seed 42 \
    --length 6000,3000 \
    --error_model nanopore2023 \
    --qscore_model nanopore2023 \
    | head -n $((HCV_READS * 4)) \
    > single_gt_hcv.fastq

# Host contamination (~1% = 5 reads)
badread simulate \
    --reference test_data/mini_host.fasta \
    --quantity 5x \
    --seed 43 \
    --length 1000,500 \
    --error_model nanopore2023 \
    | head -n 20 \
    > single_gt_host.fastq

cat single_gt_hcv.fastq single_gt_host.fastq | shuf --random-source=<(seq 999) | pigz > test_data/single_gt.fastq.gz
```

If BadRead is not available, use a simplified simulation script `bin/simulate_reads.py` (see below).

### `test_data/mixed_gt.fastq.gz`
70% reads from the 1a reference + 30% from the 3a reference + ~1% host. Same approach as above, proportion the read counts accordingly. Use `--seed 44` and `--seed 45` for reproducibility.

### `bin/simulate_reads.py`
Fallback read simulator when BadRead is not installed. Not ONT-accurate but sufficient for pipeline logic testing.
- Takes: `--reference`, `--count`, `--seed`, `--min-length`, `--max-length`, `--error-rate` (default 0.05)
- Generates reads by: random start position on reference, random length in [min, max], introduces substitutions at `error_rate` frequency, generates quality scores as ASCII 73 (Q40 synthetic)
- Outputs FASTQ to stdout

### `test_data/samplesheet.csv`
```
sample_id,fastq,metadata_json
single_gt,test_data/single_gt.fastq.gz,
mixed_gt,test_data/mixed_gt.fastq.gz,{"note":"synthetic_mixed"}
```

### `test_data/expected_outputs/single_gt_genotype_summary.json`
```json
{
  "primary_genotype": "1",
  "is_mixed": false,
  "secondary_genotypes": []
}
```
(Tolerances applied during comparison: genotype fractions ±5%.)

### `test_data/expected_outputs/mixed_gt_genotype_summary.json`
```json
{
  "primary_genotype": "1",
  "is_mixed": true,
  "secondary_genotypes": ["3"]
}
```

### `tests/test_integration.sh`
```bash
#!/usr/bin/env bash
set -euo pipefail

OUTDIR="test_out_$(date +%Y%m%d_%H%M%S)"

echo "=== Running integration test ==="
nextflow run main.nf \
    -profile test \
    --outdir "$OUTDIR" \
    2>&1 | tee nextflow_test.log

if [ "${PIPESTATUS[0]}" -ne 0 ]; then
    echo "FAIL: Nextflow run exited non-zero" >&2
    tail -50 nextflow_test.log >&2
    exit 1
fi

echo "=== Validating outputs ==="
python tests/check_outputs.py "$OUTDIR" test_data/expected_outputs

echo "PASS: Integration test complete"
```

### `tests/check_outputs.py`
Validates pipeline outputs against expected values within tolerances.

Check for each sample:
1. `genotype_summary.json` exists
2. `primary_genotype` matches exactly
3. `is_mixed` matches exactly (boolean)
4. `secondary_genotypes` is a superset of expected secondary genotypes
5. `consensus/${GT}/consensus.fasta` exists for each expected branch
6. `variants/${GT}/lofreq.filtered.vcf.gz` exists for each expected branch

Tolerances:
- Genotype fractions: ±5 percentage points
- Variant counts: within 20% of expected count (if expected count > 0)

Output format: one line per check (`PASS` or `FAIL: description`). Exit non-zero if any FAIL.

## Runtime target
The full `nextflow run main.nf -profile test` must complete in under 5 minutes on a laptop. If BadRead simulation is too slow, use `bin/simulate_reads.py` instead — simple error-rate simulation is sufficient for pipeline logic testing.

## Seed pinning
All simulation seeds must be fixed constants. Document them in `test_data/README.md`:
- HCV 1a reads: seed 42
- HCV 3a reads (mixed sample): seed 44
- Host reads: seed 43, 45
- rasusa LoFreq subsample: 42 (from `params.rasusa_seed_lofreq`)
- rasusa DEVIDER subsample: 43 (from `params.rasusa_seed_devider`)
