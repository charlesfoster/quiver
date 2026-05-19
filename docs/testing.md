# Testing and Validation Strategy — QuIVER

## 9.1 Unit-level (per module)

Each module should be invokable in isolation via `nextflow run modules/local/<module>.nf` with a tiny test input. Include a `test` workflow within each module file for this purpose.

## 9.2 Synthetic test dataset

The dataset is deliberately small (≤500 reads per sample) so a full integration run completes in <5 min.

Expected outputs and tolerances:
- `single_gt` sample:
  - `is_mixed = false`
  - `primary_genotype = "1"`
  - At least 1 variant called at simulated variant positions
  - DEVIDER produces ≥1 haplotype
- `mixed_gt` sample:
  - `is_mixed = true`
  - `primary_genotype = "1"`, `secondary_genotypes = ["3"]`
  - Two consensus FASTAs produced (one per genotype)
  - Each branch produces its own variants + haplotypes

## 9.3 Integration test harness

`tests/test_integration.sh`:
```bash
#!/usr/bin/env bash
set -euo pipefail
nextflow run main.nf -profile test --outdir test_out
python tests/check_outputs.py test_out test_data/expected_outputs
```

`tests/check_outputs.py` validates JSON outputs against expected values within tolerance (genotype fractions ±5%, variant counts ±20% of expected).

## 9.4 Smoke test on real data

For each compute profile, the smoke test command (documented in `docs/usage.md`):
```
nextflow run main.nf -profile docker --input samplesheet.csv \
    --reference_panel assets/hcv_references.fasta \
    --host_reference /path/to/GRCh38.fa.gz \
    --outdir results/smoke -resume
```

Manual validation criteria for a known sample: variants at known polymorphic positions called within ±2% of expected AF; primary genotype matches prior typing.

## 9.5 Reproducibility test

Run the pipeline twice with identical inputs, same seeds (`rasusa_seed_lofreq`, `rasusa_seed_devider`), same profile. All variant calls and consensus sequences must be byte-identical. DEVIDER haplotype outputs may differ in tie-breaking; document accepted variance.

## 9.6 Failure-mode tests

Specific tests for each failure mode (see `docs/architecture_reasoning.md` Section 12):
- Empty FASTQ → graceful exit with `EMPTY_INPUT`.
- All-host FASTQ → `NO_VIRAL_READS_LIKELY` flag, no crash.
- Sample with no HCV reads → `NO_HCV_DETECTED`, no downstream errors.
- Low-coverage sample (<100× mean) → DEVIDER skipped, LoFreq still runs flagged.
