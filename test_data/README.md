# test_data/ — Synthetic test dataset for hcv-quasi

All files are generated deterministically by `tests/generate_test_data.py`.

## Seed pinning

| File / step                         | Seed |
|-------------------------------------|------|
| HCV 1a reference sequence           |   42 |
| HCV 2b reference (25% mut from 1a)  |   43 |
| HCV 3a reference (25% mut from 1a)  |   44 |
| Human host sequence                 |   45 |
| single_gt reads (1a + host)         |  100 |
| mixed_gt reads (1a + 3a + host)     |  200 |
| rasusa LoFreq subsample             |   42 (`params.rasusa_seed_lofreq`)  |
| rasusa DEVIDER subsample            |   43 (`params.rasusa_seed_devider`) |

## Contents

| File                                          | Description                                  |
|-----------------------------------------------|----------------------------------------------|
| `mini_panel.fasta`                            | 3 synthetic HCV references (1a, 2b, 3a)      |
| `mini_host.fasta`                             | 10 kb synthetic human chr1 window            |
| `single_gt.fastq.gz`                          | 200 reads: ~98% HCV-1a + ~2% host           |
| `mixed_gt.fastq.gz`                           | 300 reads: ~69% HCV-1a + 30% HCV-3a + 1% host |
| `samplesheet.csv`                             | Pipeline input samplesheet                   |
| `expected_outputs/single_gt_genotype_summary.json` | Expected genotype call for single-gt    |
| `expected_outputs/mixed_gt_genotype_summary.json`  | Expected genotype call for mixed-gt     |

## Regenerating

```bash
python3 tests/generate_test_data.py --output-dir test_data/
```
