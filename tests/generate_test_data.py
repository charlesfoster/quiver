#!/usr/bin/env python3
"""
generate_test_data.py — Deterministic synthetic test data generator for QuIVER.

Run once to create all files under test_data/:
    python3 tests/generate_test_data.py --output-dir test_data/

Seeds used:
    42  — HCV 1a reference sequence generation
    43  — HCV 2b reference (mutate 25% of 1a positions)
    44  — HCV 3a reference (mutate 25% of 1a positions with different pattern)
    45  — Human host synthetic sequence (40% GC)
   100  — single_gt sample reads (1a + host)
   200  — mixed_gt sample reads (1a + 3a + host)
"""

import argparse
import gzip
import os
import random
import sys

# ── Constants ──────────────────────────────────────────────────────────────────
HCV_LEN   = 9646   # canonical HCV genome length
HOST_LEN  = 10000  # synthetic human chr1 window

SEED_1A   = 42
SEED_2B   = 43
SEED_3A   = 44
SEED_HOST = 45
SEED_SINGLE = 100
SEED_MIXED  = 200

# ── Sequence generation ────────────────────────────────────────────────────────

def random_seq(length: int, gc_frac: float, seed: int) -> str:
    """Generate a random DNA sequence with approximate GC fraction."""
    rng = random.Random(seed)
    at_frac = 1.0 - gc_frac
    # Build weighted nucleotide pool
    gc_weight = gc_frac / 2.0
    at_weight = at_frac / 2.0
    pool = (
        ['G'] * round(gc_weight * 1000) +
        ['C'] * round(gc_weight * 1000) +
        ['A'] * round(at_weight * 1000) +
        ['T'] * round(at_weight * 1000)
    )
    return ''.join(rng.choice(pool) for _ in range(length))


def mutate_seq(seq: str, mutation_rate: float, seed: int) -> str:
    """Apply point mutations at the given rate, preserving GC neighbourhood."""
    rng = random.Random(seed)
    bases = list(seq)
    alts = {'A': ['T', 'G', 'C'], 'T': ['A', 'G', 'C'],
            'G': ['A', 'T', 'C'], 'C': ['A', 'T', 'G']}
    for i in range(len(bases)):
        if rng.random() < mutation_rate:
            bases[i] = rng.choice(alts.get(bases[i].upper(), ['A', 'T', 'G', 'C']))
    return ''.join(bases)


def write_fasta(path: str, records: list):
    """Write FASTA records [(header, seq), ...] with 60-char line wrapping."""
    with open(path, 'w') as fh:
        for header, seq in records:
            fh.write(f'>{header}\n')
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + '\n')
    print(f'  Written: {path}  ({len(records)} records)')


# ── Read simulation ────────────────────────────────────────────────────────────

COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def sim_read(ref: str, rng: random.Random, min_len: int, max_len: int,
             error_rate: float) -> tuple:
    """
    Simulate one ONT read from *ref*.
    Returns (sequence, quality_string).
    """
    rlen = rng.randint(min_len, max_len)
    rlen = min(rlen, len(ref))
    start = rng.randint(0, len(ref) - rlen)
    read = ref[start:start + rlen]

    # Reverse-complement 50% of reads (ONT reads both strands)
    if rng.random() < 0.5:
        read = revcomp(read)

    # Introduce errors
    alts = {'A': ['T', 'G', 'C'], 'T': ['A', 'G', 'C'],
            'G': ['A', 'T', 'C'], 'C': ['A', 'T', 'G']}
    bases = list(read)
    for i in range(len(bases)):
        if rng.random() < error_rate:
            b = bases[i].upper()
            bases[i] = rng.choice(alts.get(b, ['A', 'T', 'G', 'C']))
    read = ''.join(bases)

    # Quality: Q15–Q25 range, ASCII 48–58 ('0'–':')
    qual = ''.join(chr(rng.randint(48, 58)) for _ in range(len(read)))
    return read, qual


def simulate_reads(
    refs: dict,
    counts: dict,
    seed: int,
    min_len: int = 150,
    max_len: int = 8000,
    error_rate: float = 0.05,
) -> list:
    """
    Return list of FASTQ lines (4 lines per read) in the order read-by-read,
    interleaved by ref group to ensure all refs contribute before mixing.
    """
    rng = random.Random(seed)
    fastq_lines = []
    i = 0
    for ref_name, seq in refs.items():
        n = counts[ref_name]
        for _ in range(n):
            read_seq, qual = sim_read(seq, rng, min_len, max_len, error_rate)
            fastq_lines.append(f'@read_{i} {ref_name}\n')
            fastq_lines.append(read_seq + '\n')
            fastq_lines.append('+\n')
            fastq_lines.append(qual + '\n')
            i += 1
    # Shuffle reads to interleave reference origins
    reads = [fastq_lines[j:j+4] for j in range(0, len(fastq_lines), 4)]
    rng2 = random.Random(seed + 1)
    rng2.shuffle(reads)
    result = []
    for r in reads:
        result.extend(r)
    return result


def write_fastq_gz(path: str, lines: list):
    with gzip.open(path, 'wt', compresslevel=6) as fh:
        fh.writelines(lines)
    n_reads = len(lines) // 4
    print(f'  Written: {path}  ({n_reads} reads)')


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Generate synthetic test data for QuIVER integration test')
    parser.add_argument('--output-dir', default='test_data',
                        help='Output directory (default: test_data)')
    args = parser.parse_args()

    outdir = args.output_dir
    expected_dir = os.path.join(outdir, 'expected_outputs')
    os.makedirs(expected_dir, exist_ok=True)

    print('Generating synthetic test data...')

    # ── 1. HCV reference sequences ─────────────────────────────────────────────
    print('\n[1/5] HCV reference sequences → mini_panel.fasta')
    # 1a: base sequence ~50% GC, typical HCV composition
    seq_1a = random_seq(HCV_LEN, gc_frac=0.50, seed=SEED_1A)
    # 2b: 25% mutations from 1a (seed 43)
    seq_2b = mutate_seq(seq_1a, mutation_rate=0.25, seed=SEED_2B)
    # 3a: 25% mutations from 1a with different pattern (seed 44)
    seq_3a = mutate_seq(seq_1a, mutation_rate=0.25, seed=SEED_3A)

    mini_panel_path = os.path.join(outdir, 'mini_panel.fasta')
    write_fasta(mini_panel_path, [
        (f'1a_H77_synthetic length={HCV_LEN}', seq_1a),
        (f'2b_JFH1_synthetic length={HCV_LEN}', seq_2b),
        (f'3a_S52_synthetic length={HCV_LEN}', seq_3a),
    ])

    # ── 2. Host sequence ───────────────────────────────────────────────────────
    print('\n[2/5] Host reference → mini_host.fasta')
    # ~40% GC — lower than HCV to allow minimap2 discrimination
    seq_host = random_seq(HOST_LEN, gc_frac=0.40, seed=SEED_HOST)
    mini_host_path = os.path.join(outdir, 'mini_host.fasta')
    write_fasta(mini_host_path, [
        (f'chr1_synthetic length={HOST_LEN}', seq_host),
    ])

    # ── 3. single_gt.fastq.gz ─────────────────────────────────────────────────
    # ~200 reads: ~196 from 1a (98%) + 4 from host (2%)
    print('\n[3/5] single_gt reads → single_gt.fastq.gz')
    single_counts = {'1a': 196, 'host': 4}
    single_refs   = {'1a': seq_1a, 'host': seq_host}
    single_lines  = simulate_reads(
        refs=single_refs,
        counts=single_counts,
        seed=SEED_SINGLE,
        min_len=150,
        max_len=8000,
        error_rate=0.05,
    )
    write_fastq_gz(os.path.join(outdir, 'single_gt.fastq.gz'), single_lines)

    # ── 4. mixed_gt.fastq.gz ──────────────────────────────────────────────────
    # ~300 reads: 207 from 1a (69%) + 90 from 3a (30%) + 3 from host (1%)
    print('\n[4/5] mixed_gt reads → mixed_gt.fastq.gz')
    mixed_counts = {'1a': 207, '3a': 90, 'host': 3}
    mixed_refs   = {'1a': seq_1a, '3a': seq_3a, 'host': seq_host}
    mixed_lines  = simulate_reads(
        refs=mixed_refs,
        counts=mixed_counts,
        seed=SEED_MIXED,
        min_len=150,
        max_len=8000,
        error_rate=0.05,
    )
    write_fastq_gz(os.path.join(outdir, 'mixed_gt.fastq.gz'), mixed_lines)

    # ── 5. Samplesheet ────────────────────────────────────────────────────────
    print('\n[5/5] Samplesheet → samplesheet.csv')
    ss_path = os.path.join(outdir, 'samplesheet.csv')
    with open(ss_path, 'w') as fh:
        fh.write('sample_id,fastq,metadata_json\n')
        fh.write(f'single_gt,test_data/single_gt.fastq.gz,\n')
        fh.write('mixed_gt,test_data/mixed_gt.fastq.gz,{"note":"mixed_infection_test"}\n')
    print(f'  Written: {ss_path}')

    # ── Expected outputs ───────────────────────────────────────────────────────
    print('\n[+] Expected outputs → test_data/expected_outputs/')

    import json

    single_exp = {
        "sample_id": "single_gt",
        "is_mixed": False,
        "primary_genotype": "1a",
        "branches_to_run": ["1a"],
        "_tolerances": {
            "total_mapped_reads": {"min": 150, "max": 250},
            "genotypes[0].fraction": {"min": 0.90, "max": 1.0}
        }
    }
    single_exp_path = os.path.join(expected_dir, 'single_gt_genotype_summary.json')
    with open(single_exp_path, 'w') as fh:
        json.dump(single_exp, fh, indent=2)
    print(f'  Written: {single_exp_path}')

    mixed_exp = {
        "sample_id": "mixed_gt",
        "is_mixed": True,
        "primary_genotype": "1a",
        "_tolerances": {
            "total_mapped_reads": {"min": 200, "max": 400},
            "genotypes[0].fraction": {"min": 0.55, "max": 0.85},
            "genotypes[1].fraction": {"min": 0.15, "max": 0.45}
        }
    }
    mixed_exp_path = os.path.join(expected_dir, 'mixed_gt_genotype_summary.json')
    with open(mixed_exp_path, 'w') as fh:
        json.dump(mixed_exp, fh, indent=2)
    print(f'  Written: {mixed_exp_path}')

    # ── README ────────────────────────────────────────────────────────────────
    readme_path = os.path.join(outdir, 'README.md')
    with open(readme_path, 'w') as fh:
        fh.write("""\
# test_data/ — Synthetic test dataset for QuIVER

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
""")
    print(f'  Written: {readme_path}')

    print('\nDone. All test data files created successfully.')


if __name__ == '__main__':
    main()
