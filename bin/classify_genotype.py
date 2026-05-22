#!/usr/bin/env python3
"""
classify_genotype.py — HCV subtype classification with mixed-infection detection.

Reads Round 1 BAM (competitive minimap2 mapping against the HCV reference panel)
and emits:

    1. ${sample_id}.read_assignments.tsv
       One row per primary mapped read with columns:
           read_name  reference  subtype  genotype  AS  XS  is_ambiguous

    2. ${sample_id}.genotype_summary.json
       Schema matches `docs/configuration.md` exactly:
           {
             "sample_id":          <str>,
             "total_mapped_reads": <int>,
             "ambiguous_reads":    <int>,
             "ambiguous_fraction": <float>,
             "genotypes": [
               {"genotype": <str>,       # subtype label, e.g. "1a", "3a", "6xj"
                "major_genotype": <str>, # major genotype digit(s), e.g. "1", "6"
                "fraction": <float>,
                "reads": <int>,
                "top_subtype": <str>,    # same as "genotype" field (kept for compat)
                "top_reference": <str>}, ...
             ],
             "is_mixed":             <bool>,
             "primary_genotype":     <str | null>,   # subtype of primary branch
             "secondary_genotypes":  [<str>, ...],   # subtypes of secondary branches
             "branches_to_run":      [<str>, ...]    # subtypes to process downstream
           }

Algorithm (see docs/data_flow.md Step 5.8 and
docs/architecture_reasoning.md Section 5):

  1. Walk primary alignments only (skip secondary, supplementary, unmapped).
  2. Parse subtype from the reference name with regex `^([0-9]+[a-z]?[a-z]?)_`;
     extract major genotype as the leading digits of that subtype.
  3. Flag a read as ambiguous if both AS and XS tags are present and the gap
     between them is below `--ambiguous-delta-as`.  Ambiguous reads are
     counted separately and excluded from the subtype fraction computation
     (assigning them to either subtype would bias counts).
  4. Compute per-subtype reads / fractions, identify the primary subtype,
     and flag any non-primary subtype with fraction >= `--min-secondary-fraction`
     as a secondary genotype.  `branches_to_run` is the ordered list of
     subtype branches the downstream workflow should expand.  Each subtype
     becomes its own independent processing branch — including subtypes within
     the same major genotype (e.g. 1a + 1b each get their own consensus,
     variant calls, and haplotypes).

Edge cases:
  * Zero primary mapped reads in the BAM    -> exit non-zero with clear error.
  * All primary reads flagged as ambiguous  -> total_mapped=0, is_mixed=false,
    primary_genotype=null, branches_to_run=[].  Warning logged.
  * No XS tags present (single-reference-like panel)  -> skip ambiguous detection
    gracefully; log an informational note.
  * Reference whose header does not match the regex    -> reads counted toward
    `total_mapped_reads` but excluded from any subtype assignment; warning logged.
  * Ties in primary subtype                -> deterministic break: numerically
    lowest major genotype wins, then alphabetically earliest subtype letter.
    Warning logged.

Usage:
    classify_genotype.py \\
        --bam round1.bam \\
        --sample-id P001 \\
        --out-tsv P001.read_assignments.tsv \\
        --out-json P001.genotype_summary.json \\
        [--min-secondary-fraction 0.05] \\
        [--ambiguous-delta-as 20]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from typing import Optional

import pysam


# Regex matching the HCV panel header prefix.
# Group 1 = subtype token (e.g. "1a", "2b", "6xj"); the major genotype is
# the leading digit run of that token.
# The panel FASTA (assets/hcv_references.fasta) contains only named-subtype
# sequences — unsubtyped bare-major-genotype references have been removed.
SUBTYPE_RE = re.compile(r"^([0-9]+[a-z]?[a-z]?)_")


# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #

def log(msg: str) -> None:
    """Write a structured log line to stderr."""
    print(f"[classify_genotype] {msg}", file=sys.stderr)


def parse_subtype(ref_name: str) -> Optional[tuple[str, str]]:
    """
    Parse (major_genotype, subtype) from a panel reference name.

    Returns (major_genotype, subtype) on match, or None if the header does not
    follow the expected `<digits>[<letter>[<letter>]]_<accession>` convention.

    Examples:
        '1a_M62321.1'    -> ('1', '1a')
        '6xj_EF589068.1' -> ('6', '6xj')
    """
    m = SUBTYPE_RE.match(ref_name)
    if not m:
        return None
    subtype = m.group(1)
    genotype_match = re.match(r"^(\d+)", subtype)
    if not genotype_match:
        return None  # defensive; should never trigger given the outer regex
    return genotype_match.group(1), subtype


def get_tag_or_none(read: pysam.AlignedSegment, tag: str) -> Optional[int]:
    """Return the integer value of a BAM aux tag, or None if absent."""
    try:
        return int(read.get_tag(tag))
    except KeyError:
        return None


def _subtype_sort_key(item: tuple[str, int]) -> tuple[int, str, int, str]:
    """
    Deterministic sort key for (subtype, read_count) pairs.

    Primary: read count descending.
    Tie-break: major genotype integer ascending, then full subtype string
    ascending (so "1a" beats "1b" beats "2a" etc.).
    """
    st, n = item
    major = re.match(r"^(\d+)", st)
    major_int = int(major.group(1)) if major else 999
    return (-n, major_int, st)


# --------------------------------------------------------------------------- #
# Main classification routine                                                 #
# --------------------------------------------------------------------------- #

def classify(
    bam_path: str,
    sample_id: str,
    min_secondary_fraction: float,
    ambiguous_delta_as: int,
    out_tsv: str,
    out_json: str,
) -> None:
    """Classify reads in *bam_path* by subtype and write the TSV + JSON outputs."""

    # ------------------------------------------------------------------ #
    # Per-read records (also written verbatim to the TSV)                #
    # ------------------------------------------------------------------ #
    # Tally counters keyed by subtype (e.g. "1a", "3a", "6xj").
    subtype_counter: Counter[str] = Counter()
    ref_by_subtype: dict[str, Counter[str]] = defaultdict(Counter)

    total_primary = 0           # primary, mapped, non-empty references
    ambiguous_reads = 0
    unmatched_header_reads = 0  # primary reads whose reference header could not be parsed
    xs_seen = False             # for the "no XS tags" diagnostic

    with pysam.AlignmentFile(bam_path, "rb") as bam, \
         open(out_tsv, "w") as tsv_fh:

        # TSV header — kept identical to prior schema for downstream compat.
        tsv_fh.write(
            "read_name\treference\tsubtype\tgenotype\tAS\tXS\tis_ambiguous\n"
        )

        for read in bam.fetch(until_eof=True):
            # Primary alignments only — minimap2 with --secondary=no already
            # restricts to a single primary per read, but defensive filtering
            # protects against alternative inputs.
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            total_primary += 1
            ref_name = read.reference_name or ""

            as_val = get_tag_or_none(read, "AS")
            xs_val = get_tag_or_none(read, "XS")
            if xs_val is not None:
                xs_seen = True

            # Ambiguous detection: requires both AS and XS, and a small gap.
            is_ambiguous = (
                as_val is not None
                and xs_val is not None
                and (as_val - xs_val) < ambiguous_delta_as
            )

            parsed = parse_subtype(ref_name)
            if parsed is None:
                unmatched_header_reads += 1
                tsv_fh.write(
                    f"{read.query_name}\t{ref_name}\tNA\tNA\t"
                    f"{as_val if as_val is not None else 'NA'}\t"
                    f"{xs_val if xs_val is not None else 'NA'}\t"
                    f"{is_ambiguous}\n"
                )
                continue

            major_genotype, subtype = parsed

            tsv_fh.write(
                f"{read.query_name}\t{ref_name}\t{subtype}\t{major_genotype}\t"
                f"{as_val if as_val is not None else 'NA'}\t"
                f"{xs_val if xs_val is not None else 'NA'}\t"
                f"{is_ambiguous}\n"
            )

            if is_ambiguous:
                ambiguous_reads += 1
                continue

            # Assign read to its best-hit subtype.
            subtype_counter[subtype] += 1
            ref_by_subtype[subtype][ref_name] += 1

    # ------------------------------------------------------------------ #
    # Sanity check: BAM contained at least one primary alignment.        #
    # ------------------------------------------------------------------ #
    if total_primary == 0:
        log(
            f"ERROR: BAM '{bam_path}' contains no primary mapped alignments. "
            "Refusing to emit a genotype summary."
        )
        sys.exit(2)

    if not xs_seen:
        log(
            "Note: no XS tags found in any primary alignment. Ambiguous-read "
            "detection was skipped. (Expected when the panel has effectively "
            "one matching reference per read.)"
        )

    if unmatched_header_reads:
        log(
            f"WARNING: {unmatched_header_reads} primary read(s) mapped to a "
            "reference whose header did not match the expected "
            "'<genotype><letters>_<accession>' pattern. These reads are "
            "excluded from subtype assignment but counted in "
            "total_mapped_reads."
        )

    # ------------------------------------------------------------------ #
    # Compute fractions / mixed flag.                                    #
    # ------------------------------------------------------------------ #
    total_mapped_reads = sum(subtype_counter.values())

    ambig_denominator = total_mapped_reads + ambiguous_reads
    ambiguous_fraction = (
        round(ambiguous_reads / ambig_denominator, 6)
        if ambig_denominator > 0
        else 0.0
    )

    genotypes_payload: list[dict] = []
    primary_genotype: Optional[str] = None   # holds the primary subtype
    secondary_genotypes: list[str] = []
    is_mixed = False

    if total_mapped_reads == 0:
        log(
            "WARNING: Zero non-ambiguous primary reads. "
            f"(ambiguous={ambiguous_reads}, unmatched_header={unmatched_header_reads}, "
            f"total_primary={total_primary}.) "
            "Emitting an empty genotype summary."
        )
    else:
        ranked = sorted(subtype_counter.items(), key=_subtype_sort_key)

        # Tie detection (informational warning).
        if len(ranked) >= 2 and ranked[0][1] == ranked[1][1]:
            log(
                f"WARNING: Tie in primary-subtype read count "
                f"(subtype {ranked[0][0]} == subtype {ranked[1][0]} "
                f"with {ranked[0][1]} reads each). "
                f"Assigned primary subtype = {ranked[0][0]} "
                "(lowest major genotype, then alphabetically earliest subtype)."
            )

        for st, n in ranked:
            fraction = round(n / total_mapped_reads, 6)
            major_gt = re.match(r"^(\d+)", st).group(1)
            top_reference, _ = ref_by_subtype[st].most_common(1)[0]
            genotypes_payload.append(
                {
                    "genotype":      st,        # subtype label used as the branch key
                    "major_genotype": major_gt,  # e.g. "1" for "1a"/"1b"
                    "fraction":      fraction,
                    "reads":         n,
                    "top_subtype":   st,         # alias kept for downstream compat
                    "top_reference": top_reference,
                }
            )

        primary_genotype = ranked[0][0]   # e.g. "1a"

        for st, n in ranked[1:]:
            if n / total_mapped_reads >= min_secondary_fraction:
                secondary_genotypes.append(st)

        is_mixed = bool(secondary_genotypes)

    branches_to_run: list[str] = (
        [primary_genotype] + secondary_genotypes if primary_genotype else []
    )

    # ------------------------------------------------------------------ #
    # Assemble JSON                                                      #
    # ------------------------------------------------------------------ #
    summary = {
        "sample_id":           sample_id,
        "total_mapped_reads":  total_mapped_reads,
        "ambiguous_reads":     ambiguous_reads,
        "ambiguous_fraction":  ambiguous_fraction,
        "genotypes":           genotypes_payload,
        "is_mixed":            is_mixed,
        "primary_genotype":    primary_genotype,
        "secondary_genotypes": secondary_genotypes,
        "branches_to_run":     branches_to_run,
    }

    with open(out_json, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")

    # ------------------------------------------------------------------ #
    # Console summary (stderr) for human inspection                      #
    # ------------------------------------------------------------------ #
    log(
        f"sample_id={sample_id}  total_primary={total_primary}  "
        f"non_ambiguous={total_mapped_reads}  ambiguous={ambiguous_reads}  "
        f"is_mixed={is_mixed}  primary={primary_genotype}  "
        f"secondary={secondary_genotypes}  branches={branches_to_run}"
    )


# --------------------------------------------------------------------------- #
# CLI                                                                         #
# --------------------------------------------------------------------------- #

def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse and return command-line arguments."""
    p = argparse.ArgumentParser(
        description=(
            "Classify HCV reads by subtype from a Round-1 competitive-mapping "
            "BAM and emit a per-sample summary JSON + per-read TSV."
        ),
    )
    p.add_argument(
        "--bam",
        required=True,
        help="Sorted, indexed Round 1 BAM (output of MINIMAP2_ROUND1).",
    )
    p.add_argument(
        "--sample-id",
        required=True,
        help="Sample identifier; copied into the JSON output verbatim.",
    )
    p.add_argument(
        "--out-tsv",
        required=True,
        help="Output path for the per-read assignments TSV.",
    )
    p.add_argument(
        "--out-json",
        required=True,
        help="Output path for the per-sample genotype summary JSON.",
    )
    p.add_argument(
        "--min-secondary-fraction",
        type=float,
        default=0.05,
        help=(
            "Minimum fraction of non-ambiguous reads a non-primary subtype "
            "must reach to trigger the mixed-infection flag (default: 0.05)."
        ),
    )
    p.add_argument(
        "--ambiguous-delta-as",
        type=int,
        default=20,
        help=(
            "A read is flagged as ambiguous when (AS - XS) < this value. "
            "Ambiguous reads are excluded from per-subtype counts "
            "(default: 20)."
        ),
    )
    p.add_argument(
        "--panel-fasta",
        required=False,
        default=None,
        help=(
            "(Optional) Panel FASTA. Accepted for invocation compatibility; "
            "subtype/genotype is parsed from reference names already present "
            "in the BAM header."
        ),
    )
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)
    if args.min_secondary_fraction < 0 or args.min_secondary_fraction > 1:
        log(
            f"ERROR: --min-secondary-fraction must be in [0, 1], got "
            f"{args.min_secondary_fraction}."
        )
        sys.exit(2)
    if args.ambiguous_delta_as < 0:
        log(
            f"ERROR: --ambiguous-delta-as must be >= 0, got "
            f"{args.ambiguous_delta_as}."
        )
        sys.exit(2)

    classify(
        bam_path=args.bam,
        sample_id=args.sample_id,
        min_secondary_fraction=args.min_secondary_fraction,
        ambiguous_delta_as=args.ambiguous_delta_as,
        out_tsv=args.out_tsv,
        out_json=args.out_json,
    )


if __name__ == "__main__":
    main()
