#!/usr/bin/env python3
"""
render_run_summary.py — Run-level HTML and JSON aggregation report.

Implements Step 5.22 (MultiQC aggregation) supporting layer: aggregates all
per-sample JSON summaries produced by render_sample_report.py into a single
run-level dashboard HTML and JSON file.

Usage
-----
    render_run_summary.py \\
        --sample-jsons <sid1_summary.json> [<sid2_summary.json> ...] \\
        --output-html <run_summary.html> \\
        --output-json <run_summary.json> \\
        [--pipeline-version <string>] \\
        [--report-dir <path>]   # relative base for per-sample HTML links

Output HTML
-----------
Self-contained, single-file HTML with inline CSS.  No external CDN dependencies.
Contains:
    - Run statistics block: total samples, pass rate, mean/median coverage
    - One-row-per-sample table with: sample_id, status, primary_genotype,
      is_mixed, mean_coverage, variants_called, haplotypes_reconstructed,
      link to per-sample report
    - Status colour-coding: green (PASS), orange (MIXED), red (NO_HCV /
      LOW_COVERAGE)

Output JSON schema
------------------
{
  "run_date":         str,
  "pipeline_version": str,
  "total_samples":    int,
  "pass_count":       int,
  "pass_rate":        float,
  "mean_coverage":    float | null,
  "median_coverage":  float | null,
  "samples": [
    {
      "sample_id":              str,
      "status":                 str,
      "primary_genotype":       str | null,
      "is_mixed":               bool,
      "mean_coverage":          float | null,
      "variants_called":        int | null,
      "haplotypes_reconstructed": int | null,
      "report_link":            str | null,
      "flags":                  list[str]
    }
  ]
}
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from datetime import date
from pathlib import Path


# ---------------------------------------------------------------------------
# Status colour map
# ---------------------------------------------------------------------------
STATUS_COLOURS: dict[str, dict[str, str]] = {
    "PASS":         {"bg": "#d4edda", "text": "#155724", "border": "#c3e6cb", "dot": "#28a745"},
    "MIXED":        {"bg": "#fde8dc", "text": "#7c2d12", "border": "#f1986a", "dot": "#fd7e14"},
    "NO_HCV":       {"bg": "#f8d7da", "text": "#721c24", "border": "#f5c6cb", "dot": "#dc3545"},
    "LOW_COVERAGE": {"bg": "#fff3cd", "text": "#856404", "border": "#ffc107", "dot": "#ffc107"},
    "UNKNOWN":      {"bg": "#e2e3e5", "text": "#383d41", "border": "#d6d8db", "dot": "#6c757d"},
}

STATUS_SORT_ORDER = {"NO_HCV": 0, "LOW_COVERAGE": 1, "MIXED": 2, "PASS": 3, "UNKNOWN": 4}


# ---------------------------------------------------------------------------
# Load per-sample JSON
# ---------------------------------------------------------------------------

def load_sample_json(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text())
    except Exception as exc:
        print(f"WARNING: Cannot parse {path}: {exc}", file=sys.stderr)
        return None


def extract_sample_row(data: dict) -> dict:
    """Build a flat row dict for the run summary table."""
    sid = data.get("sample_id", "unknown")
    status = data.get("overall_status", "UNKNOWN")
    primary_gt = data.get("primary_genotype")
    is_mixed = data.get("is_mixed", False)
    flags = [f["name"] for f in data.get("flags", [])]

    # Aggregate coverage across branches
    cov_values: list[float] = []
    total_variants = 0
    total_haplotypes = 0
    for branch in data.get("branches", []):
        c = branch.get("coverage_mean")
        if c is not None:
            cov_values.append(float(c))
        vc = branch.get("variant_count")
        if vc is not None:
            total_variants += int(vc)
        hc = branch.get("haplotype_count")
        if hc is not None:
            total_haplotypes += int(hc)

    mean_cov = round(sum(cov_values) / len(cov_values), 1) if cov_values else None

    # Build relative link to per-sample report.
    # run_summary.html lives in <outdir>/reports/; sample reports live in
    # <outdir>/<sample_id>/reports/ — so the relative path is always one
    # directory up, then into the sample's reports/ subdirectory.
    report_link = f"../{sid}/reports/{sid}_summary.html"

    return {
        "sample_id":               sid,
        "status":                  status,
        "primary_genotype":        primary_gt,
        "is_mixed":                is_mixed,
        "mean_coverage":           mean_cov,
        "variants_called":         total_variants if total_variants > 0 else None,
        "haplotypes_reconstructed": total_haplotypes if total_haplotypes > 0 else None,
        "report_link":             report_link,
        "flags":                   flags,
    }


# ---------------------------------------------------------------------------
# Run-level statistics
# ---------------------------------------------------------------------------

def compute_run_stats(rows: list[dict]) -> dict:
    total = len(rows)
    pass_count = sum(1 for r in rows if r["status"] == "PASS")
    pass_rate  = round(pass_count / total, 3) if total > 0 else 0.0

    cov_values = [r["mean_coverage"] for r in rows if r["mean_coverage"] is not None]
    mean_cov: float | None = None
    median_cov: float | None = None
    if cov_values:
        mean_cov   = round(sum(cov_values) / len(cov_values), 1)
        median_cov = round(statistics.median(cov_values), 1)

    mixed_count    = sum(1 for r in rows if r["status"] == "MIXED")
    no_hcv_count   = sum(1 for r in rows if r["status"] == "NO_HCV")
    low_cov_count  = sum(1 for r in rows if r["status"] == "LOW_COVERAGE")

    return {
        "total_samples":  total,
        "pass_count":     pass_count,
        "pass_rate":      pass_rate,
        "mixed_count":    mixed_count,
        "no_hcv_count":   no_hcv_count,
        "low_cov_count":  low_cov_count,
        "mean_coverage":  mean_cov,
        "median_coverage": median_cov,
    }


# ---------------------------------------------------------------------------
# HTML generation
# ---------------------------------------------------------------------------

HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>QuIVER | Run Summary</title>
<style>
*, *::before, *::after {{ box-sizing: border-box; }}
body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                 Helvetica, Arial, sans-serif;
    font-size: 14px;
    line-height: 1.5;
    color: #1a1a2e;
    background: #f5f5f5;
    margin: 0;
    padding: 0;
}}
.page-wrapper {{
    max-width: 1400px;
    margin: 0 auto;
    padding: 24px 16px 48px;
}}
.report-header {{
    background: #1a1a2e;
    color: #e8e8f0;
    padding: 24px 32px 20px;
    border-radius: 8px 8px 0 0;
}}
.report-header h1 {{
    margin: 0 0 6px;
    font-size: 22px;
    font-weight: 700;
}}
.report-header .meta-line {{
    font-size: 12px;
    color: #a0a0c0;
}}
.report-header .meta-line span {{ margin-right: 18px; }}
.section {{
    background: #ffffff;
    border: 1px solid #dde1ea;
    margin-bottom: 2px;
    padding: 20px 24px;
}}
.section:last-child {{ border-radius: 0 0 8px 8px; margin-bottom: 0; }}
.section-title {{
    font-size: 16px;
    font-weight: 700;
    color: #1a1a2e;
    margin: 0 0 14px;
    padding-bottom: 6px;
    border-bottom: 2px solid #e8e8f0;
}}
.stat-grid {{
    display: flex;
    flex-wrap: wrap;
    gap: 12px;
    margin-bottom: 4px;
}}
.stat-card {{
    background: #f0f2f8;
    border-radius: 6px;
    padding: 10px 16px;
    min-width: 130px;
    flex: 1 1 130px;
}}
.stat-card .stat-value {{
    font-size: 24px;
    font-weight: 700;
    color: #1a1a2e;
    line-height: 1.2;
}}
.stat-card .stat-label {{
    font-size: 11px;
    color: #6060a0;
    text-transform: uppercase;
    letter-spacing: 0.5px;
}}
table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 13px;
}}
th {{
    background: #f0f2f8;
    font-weight: 600;
    text-align: left;
    padding: 8px 10px;
    border-bottom: 2px solid #c8cce0;
    color: #2a2a4e;
    white-space: nowrap;
}}
td {{
    padding: 7px 10px;
    border-bottom: 1px solid #e8eaf2;
    vertical-align: middle;
}}
tr:last-child td {{ border-bottom: none; }}
tr:hover td {{ background: #f8f9ff; }}
.status-badge {{
    display: inline-flex;
    align-items: center;
    gap: 6px;
    padding: 3px 10px;
    border-radius: 12px;
    font-size: 12px;
    font-weight: 700;
    white-space: nowrap;
}}
.status-dot {{
    width: 8px;
    height: 8px;
    border-radius: 50%;
    display: inline-block;
    flex-shrink: 0;
}}
.mixed-badge {{
    display: inline-block;
    background: #b5451b;
    color: #fff;
    font-size: 10px;
    padding: 1px 6px;
    border-radius: 10px;
    margin-left: 4px;
    font-weight: 600;
}}
a.report-link {{
    color: #3a5bd9;
    text-decoration: none;
    font-size: 12px;
}}
a.report-link:hover {{ text-decoration: underline; }}
.na {{ color: #aaa; font-style: italic; }}
.params-meta {{ font-size: 13px; margin-bottom: 12px; line-height: 1.8; }}
.params-meta code {{
    display: block;
    background: #f0f2f8;
    padding: 8px 12px;
    border-radius: 4px;
    font-size: 11px;
    word-break: break-all;
    margin: 4px 0 10px;
    border: 1px solid #dde1ea;
}}
.params-table {{ font-size: 12px; margin-top: 8px; width: auto; }}
.params-table th {{ font-size: 11px; }}
.params-table td:first-child {{ font-family: monospace; color: #3a5bd9; white-space: nowrap; padding-right: 20px; }}
details.params-details {{ margin-top: 12px; }}
details.params-details > summary {{
    cursor: pointer;
    font-size: 12px;
    color: #3a5bd9;
    padding: 4px 0;
    user-select: none;
}}
details.params-details > summary:hover {{ text-decoration: underline; }}
</style>
</head>
<body>
<div class="page-wrapper">
  <div class="report-header">
    <h1>QuIVER &mdash; Run Summary</h1>
    <div class="meta-line">
      <span>Run date: {run_date}</span>
      <span>Pipeline: {pipeline_version}</span>
      <span>Samples: {total_samples}</span>
    </div>
  </div>

  <div class="section">
    <div class="section-title">Run Statistics</div>
    <div class="stat-grid">
      <div class="stat-card">
        <div class="stat-value">{total_samples}</div>
        <div class="stat-label">Total Samples</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">{pass_count}</div>
        <div class="stat-label">PASS</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">{pass_rate_pct}%</div>
        <div class="stat-label">Pass Rate</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">{mixed_count}</div>
        <div class="stat-label">Mixed Infection</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">{no_hcv_count}</div>
        <div class="stat-label">No HCV Detected</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">{low_cov_count}</div>
        <div class="stat-label">Low Coverage</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">{mean_cov_display}</div>
        <div class="stat-label">Mean Coverage (all PASS)</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">{median_cov_display}</div>
        <div class="stat-label">Median Coverage (all PASS)</div>
      </div>
    </div>
  </div>

  <div class="section">
    <div class="section-title">Sample Summary Table</div>
    <table>
      <thead>
        <tr>
          <th>#</th>
          <th>Sample ID</th>
          <th>Status</th>
          <th>Primary Genotype</th>
          <th>Mixed?</th>
          <th>Mean Coverage</th>
          <th>Variants Called</th>
          <th>Haplotypes</th>
          <th>Flags</th>
          <th>Report</th>
        </tr>
      </thead>
      <tbody>
{table_rows}
      </tbody>
    </table>
  </div>
{params_section}
</div>
</body>
</html>
"""


_KEY_PARAMS = [
    ("min_length",             "Min read length (bp)"),
    ("max_length",             "Max read length (bp)"),
    ("min_qual",               "Min read quality (Phred)"),
    ("min_secondary_fraction", "Mixed infection threshold"),
    ("min_mean_coverage",      "Min mean coverage (QC pass)"),
    ("min_consensus_cov",      "Consensus masking depth"),
    ("lofreq_max_depth",       "LoFreq depth cap (rasusa)"),
    ("min_af",                 "Min allele frequency (filter)"),
    ("min_dp",                 "Min read depth (filter)"),
    ("devider_max_depth",      "DEVIDER depth cap (rasusa)"),
    ("devider_min_abund",      "DEVIDER min haplotype abundance"),
    ("run_clair3",             "Clair3 corroboration"),
    ("skip_host_depletion",    "Host depletion skipped"),
]


def _build_params_section(run_info: dict | None) -> str:
    if not run_info:
        return ""
    import html as _html
    params   = run_info.get("params", {})
    cmd      = _html.escape(str(run_info.get("command_line", "")))
    run_name = _html.escape(str(run_info.get("run_name", "")))
    nf_ver   = _html.escape(str(run_info.get("nextflow_version", "")))

    key_rows = "\n".join(
        f'        <tr><td>{k}</td>'
        f'<td>{_html.escape(str(params.get(k, "")))}</td>'
        f'<td style="color:#666;font-size:11px;">{label}</td></tr>'
        for k, label in _KEY_PARAMS
        if k in params
    )
    all_rows = "\n".join(
        f'        <tr><td>{_html.escape(str(k))}</td>'
        f'<td style="word-break:break-all;">{_html.escape(str(v))}</td></tr>'
        for k, v in sorted(params.items())
    )
    meta_lines = []
    if cmd:
        meta_lines.append(f'      <strong>Command:</strong>\n      <code>{cmd}</code>')
    if run_name:
        meta_lines.append(f'      <strong>Run name:</strong> {run_name}')
    meta_lines.append(f'      <strong>Nextflow:</strong> {nf_ver}')
    meta_html = '\n'.join(meta_lines)

    return (
        '  <div class="section">\n'
        '    <div class="section-title">Run Parameters</div>\n'
        '    <div class="params-meta">\n'
        f'{meta_html}\n'
        '    </div>\n'
        '    <table class="params-table">\n'
        '      <thead><tr><th>Parameter</th><th>Value</th><th>Description</th></tr></thead>\n'
        f'      <tbody>\n{key_rows}\n      </tbody>\n'
        '    </table>\n'
        f'    <details class="params-details">\n'
        f'      <summary>All parameters ({len(params)} total)</summary>\n'
        '      <table class="params-table" style="margin-top:8px;">\n'
        '        <thead><tr><th>Parameter</th><th>Value</th></tr></thead>\n'
        f'        <tbody>\n{all_rows}\n        </tbody>\n'
        '      </table>\n'
        '    </details>\n'
        '  </div>'
    )


def _status_badge_html(status: str) -> str:
    c = STATUS_COLOURS.get(status, STATUS_COLOURS["UNKNOWN"])
    return (
        f'<span class="status-badge" style="background:{c["bg"]};color:{c["text"]};'
        f'border:1px solid {c["border"]};">'
        f'<span class="status-dot" style="background:{c["dot"]};"></span>'
        f'{status}</span>'
    )


def _build_table_rows(rows: list[dict], sample_order: list[str] | None = None) -> str:
    lines: list[str] = []
    if sample_order:
        order_index = {sid: i for i, sid in enumerate(sample_order)}
        sorted_rows = sorted(
            rows,
            key=lambda r: (
                order_index.get(r["sample_id"], len(sample_order)),
                r["sample_id"],
            ),
        )
    else:
        sorted_rows = sorted(
            rows,
            key=lambda r: (
                STATUS_SORT_ORDER.get(r["status"], 99),
                r["sample_id"],
            ),
        )
    for i, r in enumerate(sorted_rows, 1):
        status_html = _status_badge_html(r["status"])
        mixed_html  = ('<span class="mixed-badge">MIXED</span>'
                       if r["is_mixed"] else '<span class="na">No</span>')
        cov_html    = (f'{r["mean_coverage"]:.1f}x'
                       if r["mean_coverage"] is not None
                       else '<span class="na">n/a</span>')
        var_html    = (str(r["variants_called"])
                       if r["variants_called"] is not None
                       else '<span class="na">n/a</span>')
        hap_html    = (str(r["haplotypes_reconstructed"])
                       if r["haplotypes_reconstructed"] is not None
                       else '<span class="na">n/a</span>')
        gt_html     = (r["primary_genotype"]
                       if r["primary_genotype"]
                       else '<span class="na">n/a</span>')
        flags_html  = (", ".join(r["flags"])
                       if r["flags"]
                       else '<span class="na">none</span>')

        if r["report_link"]:
            link_html = (
                f'<a class="report-link" href="{r["report_link"]}" '
                f'target="_blank">View report</a>'
            )
        else:
            link_html = '<span class="na">n/a</span>'

        lines.append(
            f"        <tr>"
            f"<td>{i}</td>"
            f"<td><strong>{r['sample_id']}</strong></td>"
            f"<td>{status_html}</td>"
            f"<td>{gt_html}</td>"
            f"<td>{mixed_html}</td>"
            f"<td>{cov_html}</td>"
            f"<td>{var_html}</td>"
            f"<td>{hap_html}</td>"
            f"<td style='font-size:11px;color:#666;'>{flags_html}</td>"
            f"<td>{link_html}</td>"
            f"</tr>"
        )
    return "\n".join(lines)


def render_html(
    rows: list[dict],
    stats: dict,
    run_date: str,
    pipeline_version: str,
    run_info: dict | None = None,
) -> str:
    # Coverage stats for PASS samples only (for representative display)
    pass_covs = [
        r["mean_coverage"]
        for r in rows
        if r["status"] == "PASS" and r["mean_coverage"] is not None
    ]
    mean_cov_display   = f"{sum(pass_covs)/len(pass_covs):.1f}x" if pass_covs else "n/a"
    median_cov_display = f"{statistics.median(pass_covs):.1f}x"  if pass_covs else "n/a"

    sample_order = (run_info or {}).get("sample_order") or None
    table_rows = _build_table_rows(rows, sample_order)
    pass_rate_pct = round(stats["pass_rate"] * 100, 1)

    return HTML_TEMPLATE.format(
        run_date=run_date,
        pipeline_version=pipeline_version,
        total_samples=stats["total_samples"],
        pass_count=stats["pass_count"],
        pass_rate_pct=pass_rate_pct,
        mixed_count=stats["mixed_count"],
        no_hcv_count=stats["no_hcv_count"],
        low_cov_count=stats["low_cov_count"],
        mean_cov_display=mean_cov_display,
        median_cov_display=median_cov_display,
        table_rows=table_rows,
        params_section=_build_params_section(run_info),
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render run-level HTML + JSON summary for QuIVER.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--sample-jsons",
        nargs="+",
        required=True,
        help="Per-sample *_summary.json files produced by render_sample_report.py",
    )
    parser.add_argument("--output-html",        required=True)
    parser.add_argument("--output-json",        required=True)
    parser.add_argument("--pipeline-version",   default=None)
    parser.add_argument("--run-info",           default=None,
                        help="JSON file with command_line, run_name, nextflow_version, params")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    run_date        = date.today().isoformat()
    pipeline_version = args.pipeline_version or "QuIVER (development)"

    run_info: dict | None = None
    if args.run_info:
        try:
            run_info = json.loads(Path(args.run_info).read_text())
        except Exception as exc:
            print(f"WARNING: Cannot parse --run-info {args.run_info}: {exc}", file=sys.stderr)

    rows: list[dict] = []
    for json_path_str in args.sample_jsons:
        p = Path(json_path_str)
        data = load_sample_json(p)
        if data is None:
            print(f"WARNING: Skipping unreadable JSON: {p}", file=sys.stderr)
            continue
        rows.append(extract_sample_row(data))

    if not rows:
        print("ERROR: No valid per-sample JSON files were loaded.", file=sys.stderr)
        sys.exit(1)

    stats    = compute_run_stats(rows)
    html_out = render_html(rows, stats, run_date, pipeline_version, run_info)

    json_out = {
        "run_date":         run_date,
        "pipeline_version": pipeline_version,
        **stats,
        "samples":          rows,
        "run_info":         run_info,
    }

    Path(args.output_html).write_text(html_out, encoding="utf-8")
    Path(args.output_json).write_text(
        json.dumps(json_out, indent=2, default=str),
        encoding="utf-8",
    )

    print(
        f"[render_run_summary] {len(rows)} samples written to "
        f"{args.output_html} and {args.output_json}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
