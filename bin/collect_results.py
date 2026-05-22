#!/usr/bin/env python3
"""
collect_results.py — gather QuIVER outputs for sharing.

Copies HTML reports (run summary, per-sample summaries, MultiQC) and
concatenates all per-sample per-genotype haplotype FASTAs into a single
multifasta.  The output directory preserves the relative link structure
between reports so they remain clickable without a web server.

Usage:
    python collect_results.py <results_dir> [--output <share_dir>] [--zip]
"""

import argparse
import shutil
import sys
import zipfile
from pathlib import Path

_SKIP_DIRS = {"reports", "reference", "pipeline_info", ".claude"}


def find_results_dir(path: Path) -> Path:
    """Accept either the run results dir or its parent (auto-detects one level down)."""
    if (path / "reports" / "run_summary.html").exists():
        return path
    candidates = [
        d for d in path.iterdir()
        if d.is_dir() and (d / "reports" / "run_summary.html").exists()
    ]
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        print(
            f"Multiple run directories found under {path}. Pass the specific run directory.",
            file=sys.stderr,
        )
        for c in candidates:
            print(f"  {c}", file=sys.stderr)
        sys.exit(1)
    print(
        f"No run_summary.html found under {path}/reports/. Is this a QuIVER results directory?",
        file=sys.stderr,
    )
    sys.exit(1)


def sample_dirs(results_dir: Path):
    """Yield per-sample subdirectories, skipping pipeline-level dirs."""
    for d in sorted(results_dir.iterdir()):
        if d.is_dir() and d.name not in _SKIP_DIRS:
            yield d


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------

def copy_reports(results_dir: Path, out_dir: Path) -> None:
    """
    Copy HTML reports into out_dir, preserving the relative link structure:

        out_dir/
          reports/run_summary.html
          reports/multiqc_report/multiqc_report.html   (+ multiqc_data/)
          {sample_id}/reports/{sample_id}_summary.html
    """
    src_reports = results_dir / "reports"
    dst_reports = out_dir / "reports"

    # Run-level summary
    run_summary = src_reports / "run_summary.html"
    if run_summary.exists():
        dst = dst_reports / "run_summary.html"
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(run_summary, dst)
        print(f"  [run summary]  {dst.relative_to(out_dir)}")
    else:
        print(f"Warning: run_summary.html not found at {run_summary}", file=sys.stderr)

    # MultiQC report directory
    multiqc_src = src_reports / "multiqc_report"
    if multiqc_src.is_dir():
        multiqc_dst = dst_reports / "multiqc_report"
        if multiqc_dst.exists():
            shutil.rmtree(multiqc_dst)
        shutil.copytree(multiqc_src, multiqc_dst)
        print(f"  [multiqc]      {(multiqc_dst / 'multiqc_report.html').relative_to(out_dir)}")

    # Per-sample summaries
    sample_count = 0
    for sample_dir in sample_dirs(results_dir):
        sample_report = sample_dir / "reports" / f"{sample_dir.name}_summary.html"
        if not sample_report.exists():
            continue
        dst = out_dir / sample_dir.name / "reports" / sample_report.name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(sample_report, dst)
        sample_count += 1
        print(f"  [sample]       {dst.relative_to(out_dir)}")

    if sample_count == 0:
        print("Warning: no per-sample summary reports found.", file=sys.stderr)


# ---------------------------------------------------------------------------
# Haplotypes
# ---------------------------------------------------------------------------

def collect_haplotypes(results_dir: Path, out_dir: Path) -> None:
    """
    Find all {sample_id}_{GT}_haplotypes.fasta files under
    {results_dir}/{sample_id}/haplotypes/{GT}/ and concatenate them into
    out_dir/haplotypes/all_haplotypes.fasta.

    Empty FASTAs (e.g. LOW_COVERAGE samples where DEVIDER was skipped) are
    silently skipped — their absence is reflected in the per-sample HTML report.
    """
    out_fasta = out_dir / "haplotypes" / "all_haplotypes.fasta"
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    total_records = 0
    skipped = 0

    with out_fasta.open("w") as out_fh:
        for sample_dir in sample_dirs(results_dir):
            hap_root = sample_dir / "haplotypes"
            if not hap_root.is_dir():
                continue
            for gt_dir in sorted(hap_root.iterdir()):
                if not gt_dir.is_dir():
                    continue
                fasta = gt_dir / f"{sample_dir.name}_{gt_dir.name}_haplotypes.fasta"
                if not fasta.exists():
                    continue
                content = fasta.read_text()
                n_records = content.count("\n>") + (1 if content.startswith(">") else 0)
                if n_records == 0:
                    skipped += 1
                    print(
                        f"  [haplotypes]   skipped (empty) {sample_dir.name}/{gt_dir.name}",
                        file=sys.stderr,
                    )
                    continue
                # Ensure file ends with a newline before concatenating
                if content and not content.endswith("\n"):
                    content += "\n"
                out_fh.write(content)
                total_records += n_records
                print(
                    f"  [haplotypes]   {sample_dir.name}/{gt_dir.name}"
                    f"  ({n_records} haplotype{'s' if n_records != 1 else ''})"
                )

    if total_records == 0:
        out_fasta.unlink()
        out_fasta.parent.rmdir()
        print("Warning: no haplotype sequences found — haplotypes/ not created.", file=sys.stderr)
    else:
        print(
            f"  [haplotypes]   → {out_fasta.relative_to(out_dir)}"
            f"  ({total_records} total sequences)"
        )
        if skipped:
            print(f"  [haplotypes]   {skipped} branch(es) skipped (empty/LOW_COVERAGE)")


# ---------------------------------------------------------------------------
# Zip
# ---------------------------------------------------------------------------

def zip_directory(src_dir: Path, zip_path: Path) -> None:
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for file in src_dir.rglob("*"):
            if file.is_file():
                zf.write(file, file.relative_to(src_dir.parent))
    print(f"\nZipped → {zip_path}  ({zip_path.stat().st_size / 1024:.0f} KB)")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Collect QuIVER reports and haplotype FASTAs for sharing.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("results_dir", help="Path to the pipeline results directory (or parent)")
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output directory (default: <results_dir_name>_results next to results_dir)",
    )
    parser.add_argument(
        "--zip", "-z",
        action="store_true",
        help="Also create a zip archive of the output directory",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output directory if it already exists",
    )
    args = parser.parse_args()

    results_path = Path(args.results_dir).resolve()
    results_dir = find_results_dir(results_path)

    out_dir = Path(args.output).resolve() if args.output else (
        results_path.parent / f"{results_dir.name}_results"
    )

    if out_dir.exists():
        if args.overwrite:
            shutil.rmtree(out_dir)
        else:
            print(
                f"Output directory already exists: {out_dir}\nUse --overwrite to replace it.",
                file=sys.stderr,
            )
            sys.exit(1)

    out_dir.mkdir(parents=True)
    print(f"Collecting results from: {results_dir}")
    print(f"Output directory:        {out_dir}\n")

    print("Reports:")
    copy_reports(results_dir, out_dir)

    print("\nHaplotypes:")
    collect_haplotypes(results_dir, out_dir)

    if args.zip:
        zip_path = out_dir.parent / f"{out_dir.name}.zip"
        zip_directory(out_dir, zip_path)

    print(f"\nDone. Open {out_dir / 'reports' / 'run_summary.html'} to start.")


if __name__ == "__main__":
    main()
