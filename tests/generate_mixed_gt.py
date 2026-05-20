#!/usr/bin/env python3
"""
generate_mixed_gt.py — Regenerate test_data/mixed_gt.fastq.gz using badread.

Simulates ONT reads from five real HCV reference sequences drawn from
assets/hcv_references.fasta:

    1a haplotypes (3 natural strains, ~60 % of reads):
        1a_AF009606.1   9 646 bp
        1a_M67463.1     9 416 bp
        1a_EF407457.1   9 286 bp

    3a haplotypes (2 natural strains, ~40 % of reads):
        3a_D17763.1     9 456 bp
        3a_D28917.1     9 454 bp

The natural divergence between strains within each genotype (~5–10 %) gives
DEVIDER enough signal to reconstruct per-genotype haplotypes.

Usage:
    python3 tests/generate_mixed_gt.py \\
        [--panel  assets/hcv_references.fasta] \\
        [--output test_data/mixed_gt.fastq.gz] \\
        [--badread /path/to/badread]

Parameters (fixed for reproducibility):
    --quantity  100x    ≈ 590 reads across all five references
    --length    8000,2000
    --seed      200
    --random_reads 2    2 % junk reads (non-HCV contamination stand-in)
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

WANTED = [
    "1a_AF009606.1",
    "1a_M67463.1",
    "1a_EF407457.1",
    "3a_D17763.1",
    "3a_D28917.1",
]

BADREAD_DEFAULTS = [
    "/Users/z3533036/Documents/H2Seq/quasispecies_simulations/.pixi/envs/default/bin/badread",
]


def find_badread(hint: str | None) -> str:
    if hint:
        if not os.path.isfile(hint):
            sys.exit(f"ERROR: badread not found at {hint}")
        return hint
    which = shutil.which("badread")
    if which:
        return which
    for path in BADREAD_DEFAULTS:
        if os.path.isfile(path):
            return path
    sys.exit(
        "ERROR: badread not found on PATH or in default locations.\n"
        "Pass --badread /path/to/badread or activate the simulations pixi env."
    )


def extract_sequences(panel: str, wanted: list[str]) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur_id: str | None = None
    chunks: list[str] = []

    with open(panel) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_id and cur_id in wanted:
                    seqs[cur_id] = "".join(chunks)
                cur_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if cur_id and cur_id in wanted:
            seqs[cur_id] = "".join(chunks)

    missing = [w for w in wanted if w not in seqs]
    if missing:
        sys.exit(f"ERROR: Sequences not found in panel: {missing}")
    return seqs


def write_fasta(path: str, seqs: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for sid, seq in seqs.items():
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel",   default="assets/hcv_references.fasta",
                        help="Path to the HCV reference panel FASTA")
    parser.add_argument("--output",  default="test_data/mixed_gt.fastq.gz",
                        help="Output path for the gzipped FASTQ")
    parser.add_argument("--badread", default=None,
                        help="Path to the badread executable")
    args = parser.parse_args()

    badread = find_badread(args.badread)
    print(f"Using badread: {badread}")

    panel = os.path.abspath(args.panel)
    if not os.path.isfile(panel):
        sys.exit(f"ERROR: Panel not found: {panel}")

    print(f"Extracting sequences from {panel} ...")
    seqs = extract_sequences(panel, WANTED)
    for sid in WANTED:
        print(f"  {sid}  {len(seqs[sid])} bp")

    with tempfile.TemporaryDirectory() as tmpdir:
        ref_fasta = os.path.join(tmpdir, "mixed_ref.fasta")
        # Write sequences in the canonical order defined by WANTED
        write_fasta(ref_fasta, {k: seqs[k] for k in WANTED})

        raw_fastq = os.path.join(tmpdir, "reads_raw.fastq")
        cmd = [
            badread, "simulate",
            "--reference",    ref_fasta,
            "--quantity",     "100x",
            "--length",       "8000,2000",
            "--seed",         "200",
            "--random_reads", "2",
            "--error_model",  "nanopore2023",
            "--qscore_model", "nanopore2023",
        ]
        print(f"\nRunning: {' '.join(cmd)}")
        with open(raw_fastq, "w") as out_fh:
            result = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            sys.exit(f"ERROR: badread exited with code {result.returncode}")

        # Count reads
        n_reads = 0
        with open(raw_fastq) as fh:
            for line in fh:
                if line.startswith("@"):
                    n_reads += 1
        print(f"Simulated {n_reads} reads")

        # Gzip to output
        output = os.path.abspath(args.output)
        os.makedirs(os.path.dirname(output), exist_ok=True)
        print(f"Writing {output} ...")
        with open(raw_fastq, "rb") as src, gzip.open(output, "wb", compresslevel=6) as dst:
            shutil.copyfileobj(src, dst)

    size_mb = os.path.getsize(output) / 1e6
    print(f"Done. {output}  ({n_reads} reads, {size_mb:.1f} MB)")

    print("\nComposition summary (approximate, by reference weight):")
    total_len = sum(len(s) for s in seqs.values())
    for sid in WANTED:
        pct = 100 * len(seqs[sid]) / total_len
        print(f"  {sid}  {pct:.1f}%")


if __name__ == "__main__":
    main()
