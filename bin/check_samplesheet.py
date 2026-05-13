#!/usr/bin/env python3
"""
check_samplesheet.py — validate the hcv-quasi samplesheet CSV.

Usage:
    check_samplesheet.py <samplesheet.csv>

Stdout: JSON array of validated sample records (consumed by the Nextflow process).
Stderr: table summary on success, clear error message on failure.
Exit 0 on success, non-zero on any violation (fail-fast: first error stops processing).
"""

import csv
import json
import os
import re
import sys

REQUIRED_COLUMNS = {"sample_id", "fastq"}
OPTIONAL_COLUMNS = {"metadata_json"}
VALID_ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")


def die(message: str) -> None:
    """Print an error message to stderr and exit non-zero."""
    print(f"\nERROR: {message}\n", file=sys.stderr)
    sys.exit(1)


def validate_samplesheet(path: str) -> list[dict]:
    """
    Parse and validate the samplesheet CSV at *path*.

    Returns a list of dicts matching the JSON schema:
        [{"id": "P001", "fastq": "/abs/path.fastq.gz", "metadata": {...}}, ...]

    Raises SystemExit on the first validation failure.
    """
    if not os.path.isfile(path):
        die(f"Samplesheet not found: {path}")
    if not os.access(path, os.R_OK):
        die(f"Samplesheet is not readable: {path}")

    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)

        # ------------------------------------------------------------------ #
        # Column validation                                                    #
        # ------------------------------------------------------------------ #
        if reader.fieldnames is None:
            die("Samplesheet is empty or has no header row.")

        # Strip BOM / whitespace from field names (common CSV export artefact)
        fieldnames = [f.strip().lstrip("﻿") for f in reader.fieldnames]
        missing = REQUIRED_COLUMNS - set(fieldnames)
        if missing:
            die(
                f"Samplesheet is missing required column(s): {', '.join(sorted(missing))}.\n"
                f"  Found columns: {', '.join(fieldnames)}\n"
                f"  Required columns: {', '.join(sorted(REQUIRED_COLUMNS))}"
            )

        unknown = set(fieldnames) - REQUIRED_COLUMNS - OPTIONAL_COLUMNS
        if unknown:
            # Warn but do not fail — future-proofing
            print(
                f"WARNING: Ignoring unrecognised column(s): {', '.join(sorted(unknown))}",
                file=sys.stderr,
            )

        # ------------------------------------------------------------------ #
        # Row validation (fail-fast)                                          #
        # ------------------------------------------------------------------ #
        seen_ids: set[str] = set()
        validated: list[dict] = []

        for row_num, raw_row in enumerate(reader, start=2):  # row 1 = header
            # Re-map field names to stripped versions
            row = {k.strip().lstrip("﻿"): (v.strip() if v else "") for k, v in raw_row.items()}

            sample_id = row.get("sample_id", "")
            fastq_raw = row.get("fastq", "")
            metadata_raw = row.get("metadata_json", "")

            # 1. sample_id must be present
            if not sample_id:
                die(f"Row {row_num}: 'sample_id' is empty.")

            # 2. sample_id character set
            if not VALID_ID_RE.match(sample_id):
                invalid_chars = set(re.findall(r"[^\w.\-]", sample_id))
                # Also catch / and \ explicitly since \w doesn't exclude them in all locales
                problem = ", ".join(sorted(repr(c) for c in set(sample_id) - set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-")))
                die(
                    f"Row {row_num}: 'sample_id' contains invalid characters: {problem}\n"
                    f"  sample_id = '{sample_id}'\n"
                    f"  Allowed characters: A-Z a-z 0-9 . _ -"
                )

            # 3. Duplicate sample IDs
            if sample_id in seen_ids:
                die(f"Row {row_num}: Duplicate sample_id '{sample_id}'.")
            seen_ids.add(sample_id)

            # 4. fastq path must be provided
            if not fastq_raw:
                die(f"Row {row_num} (sample '{sample_id}'): 'fastq' path is empty.")

            fastq_path = os.path.abspath(fastq_raw)

            # 5. FASTQ file must exist
            if not os.path.exists(fastq_path):
                die(
                    f"Row {row_num} (sample '{sample_id}'): FASTQ file does not exist: {fastq_path}"
                )

            # 6. FASTQ file must be a regular file (not a directory, symlink to dir, etc.)
            if not os.path.isfile(fastq_path):
                die(
                    f"Row {row_num} (sample '{sample_id}'): FASTQ path is not a regular file: {fastq_path}"
                )

            # 7. FASTQ must be readable
            if not os.access(fastq_path, os.R_OK):
                die(
                    f"Row {row_num} (sample '{sample_id}'): FASTQ file is not readable: {fastq_path}"
                )

            # 8. FASTQ must not be empty (0 bytes)
            fastq_size = os.path.getsize(fastq_path)
            if fastq_size == 0:
                die(
                    f"Row {row_num} (sample '{sample_id}'): FASTQ file is empty (0 bytes): {fastq_path}\n"
                    f"  Note: downstream processes will receive an EMPTY_INPUT sentinel for this sample."
                )

            # 9. Parse optional metadata JSON (if provided)
            metadata: dict = {}
            if metadata_raw:
                try:
                    metadata = json.loads(metadata_raw)
                    if not isinstance(metadata, dict):
                        die(
                            f"Row {row_num} (sample '{sample_id}'): 'metadata_json' must be a JSON object, "
                            f"got {type(metadata).__name__}."
                        )
                except json.JSONDecodeError as exc:
                    die(
                        f"Row {row_num} (sample '{sample_id}'): 'metadata_json' is not valid JSON: {exc}\n"
                        f"  Value: {metadata_raw!r}"
                    )

            validated.append(
                {
                    "id": sample_id,
                    "fastq": fastq_path,
                    "metadata": metadata,
                }
            )

        if not validated:
            die("Samplesheet contains no data rows (only a header).")

    return validated


def print_summary(samples: list[dict]) -> None:
    """Print a human-readable table of validated samples to stderr."""
    col_id   = max(len("sample_id"), max(len(s["id"]) for s in samples))
    col_fq   = max(len("fastq"), max(len(s["fastq"]) for s in samples))
    col_meta = len("metadata")

    header = f"{'sample_id':<{col_id}}  {'fastq':<{col_fq}}  {'metadata'}"
    sep    = "-" * (col_id + 2 + col_fq + 2 + 16)

    print("", file=sys.stderr)
    print("Samplesheet validation PASSED", file=sys.stderr)
    print(sep, file=sys.stderr)
    print(header, file=sys.stderr)
    print(sep, file=sys.stderr)
    for s in samples:
        meta_str = json.dumps(s["metadata"]) if s["metadata"] else "(none)"
        print(f"{s['id']:<{col_id}}  {s['fastq']:<{col_fq}}  {meta_str}", file=sys.stderr)
    print(sep, file=sys.stderr)
    print(f"Total samples validated: {len(samples)}", file=sys.stderr)
    print("", file=sys.stderr)


def main() -> None:
    if len(sys.argv) != 2:
        print(
            f"Usage: {os.path.basename(sys.argv[0])} <samplesheet.csv>",
            file=sys.stderr,
        )
        sys.exit(1)

    samplesheet_path = sys.argv[1]
    validated = validate_samplesheet(samplesheet_path)
    print_summary(validated)

    # Write validated records as JSON to stdout (consumed by the Nextflow process)
    print(json.dumps(validated, indent=2))


if __name__ == "__main__":
    main()
