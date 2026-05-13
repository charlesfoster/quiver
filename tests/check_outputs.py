#!/usr/bin/env python3
"""
check_outputs.py — Validate hcv-quasi pipeline outputs against expected values.

Usage:
    python3 tests/check_outputs.py <results_dir> <expected_dir>

Checks for each sample in the expected_dir:
  1.  genotype_summary.json exists under results/<sample_id>/genotyping/
  2.  Exact string fields match: sample_id, primary_genotype
  3.  Exact boolean field matches: is_mixed
  4.  Numeric range checks using _tolerances from expected JSON
  5.  secondary_genotypes is a superset of expected (if present in expected)
  6.  branches_to_run matches (if present in expected)
  7.  consensus/<GT>/consensus.fasta exists per branch in branches_to_run
  8.  variants/<GT>/lofreq.filtered.vcf.gz exists per branch in branches_to_run

Exit codes:
  0 — all checks passed
  1 — one or more checks failed
  2 — usage error
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any, Optional


PASS = "PASS"
FAIL = "FAIL"


def check(condition: bool, label: str, detail: str = "") -> tuple[str, str, str]:
    status = PASS if condition else FAIL
    msg = f"{status}: {label}"
    if not condition and detail:
        msg += f" — {detail}"
    return status, label, msg


def load_json(path: Path) -> tuple:
    if not path.exists():
        return None, f"file not found: {path}"
    try:
        with open(path) as fh:
            return json.load(fh), None
    except json.JSONDecodeError as exc:
        return None, f"JSON parse error in {path}: {exc}"


def validate_sample(
    sample_id: str,
    results_dir: Path,
    expected: dict,
) -> list[tuple[str, str, str]]:
    """Run all checks for one sample. Returns list of (status, label, message)."""
    results = []

    # ── Locate genotype_summary.json ──────────────────────────────────────────
    summary_path = results_dir / sample_id / "genotyping" / f"{sample_id}.genotype_summary.json"
    status, label, msg = check(
        summary_path.exists(),
        f"{sample_id}: genotype_summary.json exists",
        f"expected at {summary_path}",
    )
    results.append((status, label, msg))
    if status == FAIL:
        return results  # no point continuing

    actual, err = load_json(summary_path)
    if err:
        results.append((FAIL, f"{sample_id}: genotype_summary.json readable", err))
        return results

    tolerances = expected.get("_tolerances", {})

    # ── Exact string fields ────────────────────────────────────────────────────
    for field in ("sample_id", "primary_genotype"):
        exp_val = expected.get(field)
        if exp_val is None:
            continue
        act_val = actual.get(field)
        results.append(check(
            act_val == exp_val,
            f"{sample_id}: {field}",
            f"expected={exp_val!r} actual={act_val!r}",
        ))

    # ── Exact boolean field ────────────────────────────────────────────────────
    if "is_mixed" in expected:
        exp_val = expected["is_mixed"]
        act_val = actual.get("is_mixed")
        results.append(check(
            act_val == exp_val,
            f"{sample_id}: is_mixed",
            f"expected={exp_val} actual={act_val}",
        ))

    # ── Numeric range checks from _tolerances ─────────────────────────────────
    for key, bounds in tolerances.items():
        lo = bounds.get("min")
        hi = bounds.get("max")
        # Resolve dotted paths like "genotypes[0].fraction"
        val = _resolve_path(actual, key)
        if val is None:
            results.append((FAIL, f"{sample_id}: {key} exists", f"path {key!r} not found in actual JSON"))
            continue
        in_range = True
        if lo is not None and val < lo:
            in_range = False
        if hi is not None and val > hi:
            in_range = False
        results.append(check(
            in_range,
            f"{sample_id}: {key} in [{lo}, {hi}]",
            f"actual={val}",
        ))

    # ── secondary_genotypes superset check ────────────────────────────────────
    if "secondary_genotypes" in expected:
        exp_secondary = set(expected["secondary_genotypes"])
        act_secondary = set(actual.get("secondary_genotypes", []))
        results.append(check(
            exp_secondary.issubset(act_secondary),
            f"{sample_id}: secondary_genotypes superset",
            f"expected to contain {sorted(exp_secondary)}, got {sorted(act_secondary)}",
        ))

    # ── branches_to_run ───────────────────────────────────────────────────────
    if "branches_to_run" in expected:
        exp_branches = expected["branches_to_run"]
        act_branches = actual.get("branches_to_run", [])
        results.append(check(
            set(exp_branches) == set(act_branches),
            f"{sample_id}: branches_to_run",
            f"expected={sorted(exp_branches)} actual={sorted(act_branches)}",
        ))
    else:
        # Use actual branches_to_run from the pipeline output for file checks
        exp_branches = actual.get("branches_to_run", [])

    # ── Per-branch file checks ─────────────────────────────────────────────────
    # Use whichever branches we expect (from expected JSON or actual output)
    check_branches = expected.get("branches_to_run", actual.get("branches_to_run", []))
    for gt in check_branches:
        consensus_path = results_dir / sample_id / "consensus" / gt / "consensus.fasta"
        results.append(check(
            consensus_path.exists(),
            f"{sample_id}/{gt}: consensus.fasta exists",
            f"expected at {consensus_path}",
        ))
        vcf_path = results_dir / sample_id / "variants" / gt / "lofreq.filtered.vcf.gz"
        results.append(check(
            vcf_path.exists(),
            f"{sample_id}/{gt}: lofreq.filtered.vcf.gz exists",
            f"expected at {vcf_path}",
        ))

    return results


def _resolve_path(data: Any, path: str) -> Any:
    """
    Resolve a dotted/indexed path in a JSON-like structure.
    Examples:
        "total_mapped_reads"       -> data["total_mapped_reads"]
        "genotypes[0].fraction"    -> data["genotypes"][0]["fraction"]
    """
    import re
    parts = re.split(r'\.', path)
    current = data
    for part in parts:
        m = re.match(r'^(\w+)\[(\d+)\]$', part)
        if m:
            key = m.group(1)
            idx = int(m.group(2))
            if not isinstance(current, dict) or key not in current:
                return None
            arr = current[key]
            if not isinstance(arr, list) or idx >= len(arr):
                return None
            current = arr[idx]
        else:
            if not isinstance(current, dict) or part not in current:
                return None
            current = current[part]
    return current


def main(argv: Optional[list] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    if len(argv) != 2:
        print(f"Usage: {sys.argv[0]} <results_dir> <expected_dir>", file=sys.stderr)
        return 2

    results_dir  = Path(argv[0])
    expected_dir = Path(argv[1])

    if not results_dir.is_dir():
        print(f"ERROR: results directory not found: {results_dir}", file=sys.stderr)
        return 2
    if not expected_dir.is_dir():
        print(f"ERROR: expected outputs directory not found: {expected_dir}", file=sys.stderr)
        return 2

    # Discover expected output files
    expected_files = sorted(expected_dir.glob("*_genotype_summary.json"))
    if not expected_files:
        print(f"ERROR: no *_genotype_summary.json files in {expected_dir}", file=sys.stderr)
        return 2

    all_results = []
    for exp_path in expected_files:
        # Derive sample_id from filename: <sample_id>_genotype_summary.json
        filename = exp_path.name
        if filename.endswith("_genotype_summary.json"):
            sample_id = filename[: -len("_genotype_summary.json")]
        else:
            sample_id = filename.split("_genotype_summary")[0]

        expected, err = load_json(exp_path)
        if err:
            all_results.append((FAIL, f"{sample_id}: load expected JSON", err))
            continue

        sample_results = validate_sample(sample_id, results_dir, expected)
        all_results.extend(sample_results)

    # ── Print results ──────────────────────────────────────────────────────────
    n_pass = 0
    n_fail = 0
    for status, _label, msg in all_results:
        print(msg)
        if status == PASS:
            n_pass += 1
        else:
            n_fail += 1

    print(f"\nSummary: {n_pass} passed, {n_fail} failed")

    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
