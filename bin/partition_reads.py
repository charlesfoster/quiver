#!/usr/bin/env python3
"""
partition_reads.py — Split reads from the Round 1 BAM into per-genotype FASTQ
files using the genotype assignments produced by ``classify_genotype.py``.

This script is only invoked when ``GENOTYPE_CLASSIFY`` flagged the sample as
mixed (``is_mixed = true`` in the genotype summary JSON). Single-genotype
samples bypass this step and feed the host-depleted FASTQ directly into the
consensus build.

Inputs:
    --bam            Round 1 BAM (sorted, indexed).  We read SEQ and QUAL
                     directly from the BAM rather than re-opening the original
                     FASTQ so the script has a single source of truth and so
                     the partitioning is deterministic with respect to what
                     was actually classified.
    --assignments    The per-read TSV emitted by ``classify_genotype.py``.
                     Columns:
                         read_name reference subtype genotype AS XS is_ambiguous
                     Only reads listed in this TSV are partitioned; primary
                     alignments missing from the TSV (e.g. unmatched-header
                     reads, or supplementary entries the classifier skipped)
                     are dropped silently.
    --sample-id      Sample identifier; used as the output filename prefix.
    --output-dir     Directory for the resulting FASTQ + summary JSON files.
    --genotypes      (Optional) comma-separated list of genotypes to emit
                     explicitly.  If omitted, the genotypes are derived from
                     the TSV.  Reads whose genotype is not in this list are
                     placed in the ambiguous pool (so callers using
                     ``branches_to_run`` from the summary JSON can pin the
                     fan-out to exactly that set).

Outputs:
    <sid>.<genotype>.fastq.gz   One file per requested genotype, containing
                                non-ambiguous reads assigned to that genotype.
                                Reads are emitted with the BAM SEQ field
                                exactly as stored — do NOT reverse-complement.
                                The SEQ field is always written in original
                                read orientation, and minimap2 will handle
                                strand again on re-mapping.
    <sid>.ambiguous.fastq.gz    Reads with ``is_ambiguous = True`` plus any
                                reads whose assigned genotype is not in the
                                requested genotype list.  Kept only for QC.
    partition_summary.json      Sidecar JSON with per-output read counts.

Algorithm:
    1. Parse the TSV into ``{read_name: (genotype, is_ambiguous)}``.
    2. Walk primary BAM alignments. Emit SEQ and QUAL verbatim from the BAM
       fields (do NOT reverse-complement). Per the Prompt 9 gotcha, the
       BAM SEQ for a reverse-strand alignment is stored relative to the
       reference forward strand; emitting it as-is yields a valid FASTQ
       and minimap2 will re-resolve strand on the next mapping, so this is
       both correct and the simplest implementation.
    3. Each read appears in exactly one output: the ambiguous pool OR one
       per-genotype FASTQ. Reads not found in the TSV are skipped silently
       (matches the Prompt 9 spec — they correspond to entries the
       classifier dropped, e.g. supplementary alignments).
    4. Compression uses Python's built-in ``gzip`` module so the script has
       no external tool dependency.

Edge cases:
    * Empty BAM (no primary alignments) — produces empty FASTQs and a summary
      with zero counts for every requested genotype.
    * Assignments TSV with zero rows — same; everything goes to ambiguous (or
      is dropped if the BAM has nothing).
    * A read appears more than once as a primary alignment — only the first
      occurrence is written (BAM should never contain this but we are
      defensive).
    * Quality string longer/shorter than the sequence — pysam can raise; we
      catch and skip that read with a warning to stderr.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Optional

import pysam


# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #

def log(msg: str) -> None:
    """Write a structured log line to stderr."""
    print(f"[partition_reads] {msg}", file=sys.stderr)


def parse_bool(value: str) -> bool:
    """Tolerant boolean parser — accepts True/False, true/false, 1/0."""
    v = value.strip().lower()
    if v in ("true", "1", "yes"):
        return True
    if v in ("false", "0", "no"):
        return False
    raise ValueError(f"Unrecognised boolean literal: {value!r}")


def load_assignments(path: str) -> dict[str, tuple[Optional[str], bool]]:
    """
    Load the assignments TSV into ``{read_name: (genotype_or_None, is_ambiguous)}``.

    Rows where genotype is ``"NA"`` (header didn't match the panel pattern)
    are represented as ``(None, is_ambiguous)`` and treated as ambiguous on
    the partitioning side: we can't assign them to a genotype branch.
    """
    out: dict[str, tuple[Optional[str], bool]] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # Expected columns:
        # read_name  reference  subtype  genotype  AS  XS  is_ambiguous
        try:
            idx_name = header.index("read_name")
            idx_gt   = header.index("genotype")
            idx_amb  = header.index("is_ambiguous")
        except ValueError as exc:
            log(
                f"ERROR: assignments TSV {path!r} is missing a required column: {exc}"
            )
            sys.exit(2)

        for line_no, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= max(idx_name, idx_gt, idx_amb):
                log(f"WARNING: skipping malformed row {line_no} in {path!r}")
                continue
            read_name = parts[idx_name]
            gt_raw    = parts[idx_gt]
            amb_raw   = parts[idx_amb]
            try:
                is_amb = parse_bool(amb_raw)
            except ValueError:
                log(
                    f"WARNING: row {line_no}: unrecognised is_ambiguous value "
                    f"{amb_raw!r}; treating as ambiguous."
                )
                is_amb = True
            genotype: Optional[str] = None if gt_raw in ("NA", "") else gt_raw
            out[read_name] = (genotype, is_amb)
    return out


# --------------------------------------------------------------------------- #
# Main partition routine                                                      #
# --------------------------------------------------------------------------- #

def partition(
    bam_path: str,
    assignments_path: str,
    sample_id: str,
    output_dir: str,
    requested_genotypes: Optional[list[str]],
) -> dict:
    """Partition the BAM into per-genotype FASTQ files."""

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    assignments = load_assignments(assignments_path)
    if not assignments:
        log(
            f"WARNING: assignments TSV {assignments_path!r} produced zero "
            "read-name entries; all output FASTQs will be empty."
        )

    # If the caller did not pin the genotype list, derive it from the TSV
    # (deterministic order: sorted ascending so output filenames are stable).
    if requested_genotypes is None:
        derived: set[str] = {
            gt for gt, amb in assignments.values()
            if gt is not None and not amb
        }
        requested_genotypes = sorted(derived)

    # Open one gzip writer per genotype + one for the ambiguous pool.
    writers: dict[str, gzip.GzipFile] = {}
    paths:   dict[str, Path]          = {}
    counts:  Counter[str]             = Counter()

    for gt in requested_genotypes:
        p = out_dir / f"{sample_id}.{gt}.fastq.gz"
        writers[gt] = gzip.open(p, "wb")
        paths[gt]   = p
        counts[gt]  = 0

    ambiguous_path = out_dir / f"{sample_id}.ambiguous.fastq.gz"
    ambiguous_writer = gzip.open(ambiguous_path, "wb")
    counts["__ambiguous__"] = 0

    seen_reads: set[str] = set()
    skipped_no_assignment   = 0
    skipped_no_seq          = 0
    skipped_quality_mismatch = 0

    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                # Primary alignments only — match the classifier's perspective.
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                read_name = read.query_name
                if read_name is None:
                    continue
                if read_name in seen_reads:
                    # Defensive: a primary alignment for the same read appearing
                    # twice would be malformed; emit once and move on.
                    continue
                seen_reads.add(read_name)

                assignment = assignments.get(read_name)
                if assignment is None:
                    # The classifier didn't record this read — skip silently
                    # (per docs/implementation_prompts.md Prompt 9 gotcha).
                    skipped_no_assignment += 1
                    continue

                genotype, is_amb = assignment

                # Pull SEQ and QUAL straight from the BAM.  Per
                # docs/implementation_prompts.md: keep BAM orientation and
                # let minimap2 re-resolve strand on re-mapping.
                seq = read.query_sequence
                if seq is None:
                    skipped_no_seq += 1
                    continue

                qual_arr = read.query_qualities
                if qual_arr is None:
                    qual = "!" * len(seq)
                else:
                    if len(qual_arr) != len(seq):
                        skipped_quality_mismatch += 1
                        log(
                            f"WARNING: read {read_name!r}: quality length "
                            f"{len(qual_arr)} != seq length {len(seq)}; skipping."
                        )
                        continue
                    qual = "".join(chr(q + 33) for q in qual_arr)

                record = (
                    f"@{read_name}\n{seq}\n+\n{qual}\n".encode("ascii")
                )

                # Routing decision:
                #   - is_ambiguous → ambiguous pool
                #   - genotype is None (unparseable ref header) → ambiguous pool
                #   - genotype not in requested set → ambiguous pool
                #   - otherwise → per-genotype writer
                if is_amb or genotype is None or genotype not in writers:
                    ambiguous_writer.write(record)
                    counts["__ambiguous__"] += 1
                else:
                    writers[genotype].write(record)
                    counts[genotype] += 1
    finally:
        for w in writers.values():
            w.close()
        ambiguous_writer.close()

    if skipped_no_assignment:
        log(
            f"Note: {skipped_no_assignment} primary BAM read(s) had no entry "
            "in the assignments TSV; skipped silently."
        )
    if skipped_no_seq:
        log(
            f"WARNING: {skipped_no_seq} primary BAM read(s) had no SEQ field; "
            "skipped."
        )
    if skipped_quality_mismatch:
        log(
            f"WARNING: {skipped_quality_mismatch} primary BAM read(s) had a "
            "quality/sequence length mismatch; skipped."
        )

    summary = {
        "sample_id":       sample_id,
        "genotypes": [
            {
                "genotype":      gt,
                "reads_written": counts[gt],
                "output_file":   str(paths[gt].name),
            }
            for gt in requested_genotypes
        ],
        "ambiguous_reads": counts["__ambiguous__"],
        "ambiguous_file":  ambiguous_path.name,
        "skipped": {
            "no_assignment":     skipped_no_assignment,
            "no_seq":            skipped_no_seq,
            "quality_mismatch":  skipped_quality_mismatch,
        },
    }

    summary_path = out_dir / "partition_summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")

    log(
        f"sample_id={sample_id}  "
        f"genotypes={ {gt: counts[gt] for gt in requested_genotypes} }  "
        f"ambiguous={counts['__ambiguous__']}"
    )

    return summary


# --------------------------------------------------------------------------- #
# CLI                                                                         #
# --------------------------------------------------------------------------- #

def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Partition Round-1 BAM reads into per-genotype FASTQ files using "
            "the assignments TSV produced by classify_genotype.py."
        ),
    )
    p.add_argument(
        "--bam",
        required=True,
        help="Sorted, indexed Round 1 BAM (output of MINIMAP2_ROUND1).",
    )
    p.add_argument(
        "--assignments",
        required=True,
        help="Per-read assignments TSV from classify_genotype.py.",
    )
    p.add_argument(
        "--sample-id",
        required=True,
        help="Sample identifier; used as the output filename prefix.",
    )
    p.add_argument(
        "--output-dir",
        required=True,
        help="Directory in which to write per-genotype FASTQs and summary JSON.",
    )
    p.add_argument(
        "--genotypes",
        required=False,
        default=None,
        help=(
            "Optional comma-separated list of genotypes to emit "
            "(e.g. '1,3').  Reads whose genotype is outside this list "
            "are routed to the ambiguous pool.  Defaults to all genotypes "
            "present in the assignments TSV."
        ),
    )
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)

    requested: Optional[list[str]] = None
    if args.genotypes is not None:
        requested = [g.strip() for g in args.genotypes.split(",") if g.strip()]
        if not requested:
            log(
                "ERROR: --genotypes given but produced an empty list after "
                "splitting on commas."
            )
            sys.exit(2)

    partition(
        bam_path=args.bam,
        assignments_path=args.assignments,
        sample_id=args.sample_id,
        output_dir=args.output_dir,
        requested_genotypes=requested,
    )


if __name__ == "__main__":
    main()
