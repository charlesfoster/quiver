#!/usr/bin/env python3
"""
nanoq-compatible FASTQ statistics — outputs JSON matching the nanoq 0.10 schema.

Replaces the nanoq binary for platforms where a native binary is unavailable
(e.g. osx-arm64). Uses only Python stdlib; no binary dependencies.

CLI is a drop-in for the flags used by the nanoq Nextflow module:
    nanoq_stats.py -i reads.fastq.gz -r output.nanoq.json

JSON schema (identical to nanoq 0.10):
    reads, bases, n50, longest, shortest,
    mean_length, median_length, mean_quality, median_quality
"""
import argparse
import gzip
import json
import math
import statistics


def _phred_mean(qual: str) -> float:
    """Mean Phred quality via error-probability averaging (matches nanoq behaviour)."""
    if not qual:
        return 0.0
    total_prob = sum(10 ** (-(ord(c) - 33) / 10.0) for c in qual)
    mean_prob = total_prob / len(qual)
    return -10.0 * math.log10(mean_prob) if mean_prob > 0 else 0.0


def _n50(lengths: list) -> int:
    if not lengths:
        return 0
    desc = sorted(lengths, reverse=True)
    half = sum(desc) / 2.0
    cumsum = 0
    for ln in desc:
        cumsum += ln
        if cumsum >= half:
            return ln
    return desc[-1]


def _parse_fastq(path: str):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            hdr = fh.readline()
            if not hdr:
                break
            seq  = fh.readline().strip()
            fh.readline()              # +
            qual = fh.readline().strip()
            yield len(seq), _phred_mean(qual)


def compute_stats(path: str) -> dict:
    lengths, quals = [], []
    for ln, q in _parse_fastq(path):
        lengths.append(ln)
        quals.append(q)

    if not lengths:
        return {k: 0 for k in (
            "reads", "bases", "n50", "longest", "shortest",
            "mean_length", "median_length", "mean_quality", "median_quality",
        )}

    return {
        "reads":          len(lengths),
        "bases":          sum(lengths),
        "n50":            _n50(lengths),
        "longest":        max(lengths),
        "shortest":       min(lengths),
        "mean_length":    round(statistics.mean(lengths), 2),
        "median_length":  round(statistics.median(lengths), 2),
        "mean_quality":   round(statistics.mean(quals), 2),
        "median_quality": round(statistics.median(quals), 2),
    }


def main():
    ap = argparse.ArgumentParser(description="nanoq-compatible FASTQ stats (pure Python)")
    ap.add_argument("-i", "--input",  required=True, help="Input FASTQ (plain or gzip)")
    ap.add_argument("-r", "--report", help="Write JSON stats to this file")
    ap.add_argument("--json", action="store_true", help="Emit JSON to stdout")
    args = ap.parse_args()

    result = compute_stats(args.input)
    out = json.dumps(result, indent=2)

    if args.report:
        with open(args.report, "w") as fh:
            fh.write(out + "\n")

    if args.json or not args.report:
        print(out)


if __name__ == "__main__":
    main()
