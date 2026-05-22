#!/usr/bin/env python3
"""
render_sample_report.py — Per-sample HTML and JSON report generator.

Implements Step 5.21 of the HCV quasispecies pipeline data flow specification.
Wraps all per-sample per-genotype outputs into a single self-contained HTML
file and a machine-readable JSON summary for the run-level aggregator.

All optional inputs are handled gracefully — the report renders even when
haplotype reconstruction failed, NanoPlot images are absent, or coverage
data is missing.

Usage
-----
    render_sample_report.py \\
        --sample-id <sid> \\
        --genotype-summary <genotype_summary.json> \\
        --nanoq-raw <raw.nanoq.json> \\
        [--nanoq-filtered <filtered.nanoq.json>] \\
        [--host-stats <host_stats.json>] \\
        [--mosdepth-summaries <file1> [<file2> ...]] \\
        [--variant-tsvs <file1> [<file2> ...]] \\
        [--flagstats <file1> [<file2> ...]] \\
        [--haplotype-reports <file1> [<file2> ...]] \\
        [--nanoplot-dirs <dir1> [<dir2> ...]] \\
        [--flags <flag1> [<flag2> ...]] \\
        --output-html <sid_summary.html> \\
        --output-json <sid_summary.json> \\
        [--template <path/to/sample_report.html.j2>] \\
        [--pipeline-version <string>]

Output JSON schema
------------------
{
  "sample_id": str,
  "run_date": str,
  "pipeline_version": str,
  "overall_status": str,        # PASS | MIXED | NO_HCV | LOW_COVERAGE
  "is_mixed": bool,
  "primary_genotype": str | null,
  "genotype_summary": { ... },  # raw genotype_summary.json contents
  "funnel": { ... },
  "branches": [ ... ],
  "flags": [ { "name": str, "explanation": str } ]
}
"""

from __future__ import annotations

import argparse
import base64
import json
import re
import sys
from datetime import date
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Optional Jinja2 import — fail with a clear message if not available.
# ---------------------------------------------------------------------------
try:
    from jinja2 import Environment, FileSystemLoader, select_autoescape
except ImportError:
    sys.exit(
        "ERROR: Jinja2 is required.  Install it with: pip install jinja2"
    )


# ---------------------------------------------------------------------------
# Flag explanations
# ---------------------------------------------------------------------------
FLAG_EXPLANATIONS: dict[str, str] = {
    "NO_HCV_DETECTED":          "Fewer than 100 reads mapped to the HCV panel in Round 1; "
                                 "sample may not contain HCV or read quality was insufficient.",
    "LOW_COVERAGE":             "Mean coverage in Round 2 mapping was below the threshold "
                                 "(default 100x); DEVIDER haplotype reconstruction was skipped.",
    "LOW_MAPPING_RATE":         "Fewer than 90% of primary reads mapped in Round 2 mapping "
                                 "against the sample-specific consensus.",
    "LOW_COVERAGE_CONSENSUS":   "Consensus quality check failed: contig length out of range or "
                                 ">30% masked positions.",
    "NO_VIRAL_READS_LIKELY":    ">=99.5% of reads mapped to the host reference; "
                                 "very few viral reads remain.",
    "ALL_READS_FILTERED":       "All reads were removed by the chopper quality/length filter.",
    "EMPTY_INPUT":              "Input FASTQ was empty (0 bytes).",
    "DEVIDER_FAILED":           "DEVIDER exited non-zero; haplotype reconstruction was skipped "
                                 "for this genotype branch.",
    "MIXED":                    "Mixed-genotype infection detected (>= 5% reads to a secondary "
                                 "major genotype).",
}


# ---------------------------------------------------------------------------
# Parsers for each input file type
# ---------------------------------------------------------------------------

def _safe_load_json(path: Path | None) -> dict | None:
    """Load a JSON file; return None on any error."""
    if path is None or not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def parse_nanoq_json(path: Path | None) -> dict | None:
    """Return nanoq stats dict with keys: reads, bases, mean_length, n50."""
    d = _safe_load_json(path)
    if d is None:
        return None
    return {
        "reads":       d.get("reads") or d.get("reads_passed") or d.get("num_reads"),
        "bases":       d.get("bases") or d.get("bases_passed") or d.get("total_bases"),
        "mean_length": d.get("mean_length") or d.get("mean_read_length"),
        "n50":         d.get("n50") or d.get("read_n50"),
    }


def parse_host_stats(path: Path | None) -> dict | None:
    """Return host stats dict."""
    d = _safe_load_json(path)
    if d is None:
        return None
    return {
        "total_reads":        d.get("total_reads"),
        "host_reads_removed": d.get("host_reads_removed"),
        "kept_reads":         d.get("kept_reads"),
        "host_fraction":      d.get("host_fraction"),
    }


def parse_genotype_summary(path: Path | None) -> dict | None:
    """Return genotype summary dict (verbatim schema from docs/configuration.md)."""
    return _safe_load_json(path)


def parse_mosdepth_summary(path: Path | None) -> dict | None:
    """
    Parse a *.mosdepth.summary.txt file.

    mosdepth summary format (tab-delimited):
        chrom   length  bases   mean    min     max

    The 'total_region' row supplies genome-wide stats.
    Returns dict with keys: mean, min, max, reference_length, regions (list).
    """
    if path is None or not path.is_file():
        return None
    lines = path.read_text().splitlines()
    result: dict[str, Any] = {
        "mean": None, "min": None, "max": None,
        "reference_length": None, "regions": []
    }
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        chrom = parts[0]
        # Columns: chrom, length, bases, mean, min, max
        # Skip header row
        if chrom == "chrom":
            continue
        length_str = parts[1]
        mean_val   = parts[3]
        min_val    = parts[4] if len(parts) > 4 else None
        max_val    = parts[5] if len(parts) > 5 else None
        try:
            mean_f = float(mean_val)
        except ValueError:
            continue
        if chrom in ("total_region", "total"):
            try:
                result["reference_length"] = int(length_str)
            except ValueError:
                pass
            result["mean"] = mean_f
            result["min"]  = float(min_val) if min_val else None
            result["max"]  = float(max_val) if max_val else None
    return result


def parse_mosdepth_bed(path: Path | None) -> list[dict] | None:
    """
    Parse a *.regions.bed.gz produced by mosdepth --by 100.

    Format (0-based, tab-delimited): chrom  start  end  coverage
    Returns a list of {chrom, start, end, mean} dicts for the coverage chart.
    """
    if path is None or not path.is_file():
        return None
    import gzip
    regions: list[dict] = []
    try:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                try:
                    regions.append({
                        "chrom": parts[0],
                        "start": int(parts[1]),
                        "end":   int(parts[2]),
                        "mean":  float(parts[3]),
                    })
                except (ValueError, IndexError):
                    continue
    except Exception:
        return None
    return regions if regions else None


def parse_variant_tsv(path: Path | None) -> list[dict] | None:
    """
    Parse a variants.tsv produced by VARIANT_FILTER.

    Columns: CHROM, POS, REF, ALT, AF, DP, SB
    Returns a list of dicts.  Returns empty list (not None) if file is present
    but has no data rows, so the template can distinguish 'no data' from
    'no variants above threshold'.
    """
    if path is None:
        return None
    if not path.is_file():
        return None
    variants: list[dict] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        try:
            variants.append({
                "chrom": parts[0],
                "pos":   int(parts[1]),
                "ref":   parts[2],
                "alt":   parts[3],
                "af":    float(parts[4]),
                "dp":    int(parts[5]) if len(parts) > 5 else None,
                "sb":    float(parts[6]) if len(parts) > 6 and parts[6] not in (".", "") else None,
            })
        except (ValueError, IndexError):
            continue
    return variants


def parse_flagstat(path: Path | None) -> dict | None:
    """
    Parse a samtools flagstat file.

    Returns dict with keys: total, primary, primary_mapped, mapping_rate (str).
    """
    if path is None or not path.is_file():
        return None
    text = path.read_text()
    result: dict[str, Any] = {
        "total": None,
        "primary": None,
        "primary_mapped": None,
        "mapping_rate": None,
    }
    for line in text.splitlines():
        # "N + 0 in total (QC-passed reads + QC-failed reads)"
        m = re.match(r"^(\d+) \+ \d+ in total", line)
        if m:
            result["total"] = int(m.group(1))
        # "N + 0 primary"  (line ends after 'primary')
        m = re.match(r"^(\d+) \+ \d+ primary$", line)
        if m:
            result["primary"] = int(m.group(1))
        # "N + 0 primary mapped (PCT% : N/A)"
        m = re.match(r"^(\d+) \+ \d+ primary mapped \(([^)]+)\)", line)
        if m:
            result["primary_mapped"] = int(m.group(1))
            result["mapping_rate"]   = m.group(2).strip()
    return result


def parse_haplotype_report(path: Path | None) -> dict | None:
    """
    Parse a *_haplotype_report.json produced by FORMAT_HAPLOTYPES.

    Normalises chain dicts to add template-friendly aliases and computes
    longest_chain_length.  JSON uses: id, length_bp, abundance_lower_bound,
    windows (list); template and JSON summary use: chain_id, total_length,
    abundance, n_windows.
    """
    d = _safe_load_json(path)
    if d is None:
        return None
    longest = None
    for chain in d.get("chains", []):
        # Add template aliases (setdefault preserves any already-correct keys).
        chain.setdefault("chain_id",     chain.get("id"))
        chain.setdefault("total_length", chain.get("length_bp"))
        chain.setdefault("abundance",    chain.get("abundance_lower_bound"))
        chain.setdefault("n_windows",    len(chain.get("windows") or []))
        tlen = chain.get("length_bp")
        if tlen is not None:
            longest = max(longest, tlen) if longest is not None else tlen
    d["longest_chain_length"] = longest
    return d


def collect_nanoplot_pngs(dirs: list[Path]) -> list[dict]:
    """
    Find PNG files in NanoPlot output directories and return them as
    base64-encoded dicts ready for embedding in HTML.

    Priority order for selection (per directory):
      1. *LengthvsQualityScatterPlot*.png
      2. *Weighted*.png
      3. Any other *.png
    Only the first 2 images per directory are embedded to keep file size sane.
    """
    images: list[dict] = []
    priority_patterns = [
        re.compile(r"LengthvsQuality", re.IGNORECASE),
        re.compile(r"WeightedHistogram", re.IGNORECASE),
        re.compile(r".*", re.IGNORECASE),  # fallback
    ]
    for nanodir in dirs:
        if not nanodir.is_dir():
            continue
        pngs = sorted(nanodir.glob("*.png"))
        selected: list[Path] = []
        for pattern in priority_patterns:
            for p in pngs:
                if pattern.search(p.name) and p not in selected:
                    selected.append(p)
                    if len(selected) >= 2:
                        break
            if len(selected) >= 2:
                break
        for png in selected[:2]:
            try:
                data = base64.b64encode(png.read_bytes()).decode("ascii")
                images.append({"name": png.name, "data": data})
            except OSError:
                continue
    return images


# ---------------------------------------------------------------------------
# SVG coverage chart
# ---------------------------------------------------------------------------

def _make_coverage_svg(regions: list[dict], max_width: int = 900, height: int = 80) -> str:
    """
    Generate an inline SVG bar chart of per-window coverage.

    Each 100-bp window is one bar.  Bars are colour-coded:
        dark-green  >= 500x
        green       100–500x
        orange      20–100x
        red         0–20x
    """
    if not regions:
        return ""
    max_cov = max((r["mean"] for r in regions), default=1.0) or 1.0
    # Cap visual max at 99th-percentile to avoid one outlier squashing everything
    sorted_covs = sorted(r["mean"] for r in regions)
    p99_idx = max(0, int(len(sorted_covs) * 0.99) - 1)
    vis_max = sorted_covs[p99_idx] if sorted_covs else max_cov
    if vis_max <= 0:
        vis_max = max_cov or 1.0

    n = len(regions)
    bar_w = max(1.0, max_width / n)
    svg_w = bar_w * n
    bars = []
    for i, r in enumerate(regions):
        cov = r["mean"]
        bar_h = max(1.0, (min(cov, vis_max) / vis_max) * height)
        y = height - bar_h
        x = i * bar_w
        if cov >= 500:
            colour = "#1a6e3f"
        elif cov >= 100:
            colour = "#2d9d5a"
        elif cov >= 20:
            colour = "#e67e22"
        else:
            colour = "#c0392b"
        tip = f"{r.get('chrom','')}: {r['start']}-{r['end']} | {cov:.1f}x"
        bars.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" '
            f'fill="{colour}"><title>{tip}</title></rect>'
        )

    # Axis labels
    labels = []
    genomic_span = regions[-1]["end"] - regions[0]["start"]
    for tick_frac in [0, 0.25, 0.5, 0.75, 1.0]:
        pos = int(regions[0]["start"] + genomic_span * tick_frac)
        x_pos = svg_w * tick_frac
        labels.append(
            f'<text x="{x_pos:.0f}" y="{height + 14}" '
            f'font-size="10" fill="#666" text-anchor="middle">{pos:,}</text>'
        )

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{min(svg_w, max_width):.0f}" height="{height + 20}" '
        f'viewBox="0 0 {svg_w:.1f} {height + 20}" preserveAspectRatio="none">'
        + "".join(bars)
        + "".join(labels)
        + "</svg>"
    )


# ---------------------------------------------------------------------------
# Genotype label extraction from filenames
# ---------------------------------------------------------------------------

def _extract_genotype_from_filename(filename: str, sample_id: str) -> str | None:
    """
    Extract the genotype label from a filename like:
        P001_1a_haplotype_report.json
        P001_2b.mosdepth.summary.txt
        P001_3a_variants.tsv
        P001_1a_round2.flagstat

    Strategy: strip the sample_id prefix plus underscore, then take the next
    token (up to the next underscore or dot).  Fallback: None.
    """
    stem = Path(filename).name
    prefix = sample_id + "_"
    if stem.startswith(prefix):
        remainder = stem[len(prefix):]
        # remainder might be "1a_haplotype_report.json" or "1a.mosdepth.summary.txt"
        token = re.split(r"[_.]", remainder)[0]
        return token if token else None
    # Try to find a genotype-like token anywhere in the name
    m = re.search(r"[_\-]([0-9]+[a-z]{0,2})[_\-.]", stem)
    return m.group(1) if m else None


# ---------------------------------------------------------------------------
# Status determination
# ---------------------------------------------------------------------------

def _determine_status(
    flags: list[dict],
    is_mixed: bool,
    primary_genotype: str | None,
) -> str:
    flag_names = {f["name"] for f in flags}
    if "NO_HCV_DETECTED" in flag_names or "EMPTY_INPUT" in flag_names or "ALL_READS_FILTERED" in flag_names:
        return "NO_HCV"
    if "LOW_COVERAGE" in flag_names:
        return "LOW_COVERAGE"
    if is_mixed:
        return "MIXED"
    return "PASS"


# ---------------------------------------------------------------------------
# Main report builder
# ---------------------------------------------------------------------------

def build_context(args: argparse.Namespace) -> dict:
    """Collect all data and return the Jinja2 template context dict."""

    run_date        = date.today().isoformat()
    pipeline_version = args.pipeline_version or "QuIVER (development)"
    sample_id       = args.sample_id

    # ---- Genotype summary ----
    gt_summary = parse_genotype_summary(
        Path(args.genotype_summary) if args.genotype_summary else None
    )
    if gt_summary is None:
        gt_summary = {
            "sample_id":           sample_id,
            "total_mapped_reads":  None,
            "ambiguous_reads":     None,
            "ambiguous_fraction":  0.0,
            "genotypes":           [],
            "is_mixed":            False,
            "primary_genotype":    None,
            "secondary_genotypes": [],
            "branches_to_run":     [],
        }

    is_mixed        = gt_summary.get("is_mixed", False)
    primary_genotype = gt_summary.get("primary_genotype")

    # Build a genotype → top_subtype lookup from the per-genotype entries.
    subtype_by_gt: dict[str, str] = {
        entry["genotype"]: entry.get("top_subtype", entry["genotype"])
        for entry in gt_summary.get("genotypes", [])
        if "genotype" in entry
    }
    primary_subtype = subtype_by_gt.get(primary_genotype) if primary_genotype else None
    detected_subtypes: list[str] = [
        subtype_by_gt.get(gt, gt)
        for gt in gt_summary.get("branches_to_run", [])
    ]

    # ---- Flags ----
    raw_flags: list[str] = []
    if args.flags:
        for f in args.flags:
            p = Path(f)
            if p.is_file():
                # Filename convention: {sample_id}.{FLAGNAME} or {sample_id}_{GT}.{FLAGNAME}
                # or just {FLAGNAME}
                parts = p.name.split(".")
                if len(parts) >= 2:
                    raw_flags.append(parts[-1])  # last extension = flag name
                else:
                    raw_flags.append(p.name)

    if is_mixed:
        raw_flags.append("MIXED")

    flags: list[dict] = []
    seen_flags: set[str] = set()
    for fname in raw_flags:
        if fname not in seen_flags:
            seen_flags.add(fname)
            flags.append({
                "name":        fname,
                "explanation": FLAG_EXPLANATIONS.get(fname, "Pipeline flag raised."),
            })

    # ---- Read funnel ----
    nanoq_raw      = parse_nanoq_json(Path(args.nanoq_raw) if args.nanoq_raw else None)
    nanoq_filtered = parse_nanoq_json(Path(args.nanoq_filtered) if args.nanoq_filtered else None)
    host_stats     = parse_host_stats(Path(args.host_stats) if args.host_stats else None)

    raw_reads      = nanoq_raw["reads"] if nanoq_raw else None
    filtered_reads = nanoq_filtered["reads"] if nanoq_filtered else None

    # host_depleted_reads comes from host_stats.kept_reads
    host_depleted_reads = host_stats["kept_reads"] if host_stats else None

    # mapped_reads: total_mapped_reads from genotype_summary
    mapped_reads = gt_summary.get("total_mapped_reads")

    funnel = {
        "raw_reads":           raw_reads,
        "filtered_reads":      filtered_reads,
        "host_depleted_reads": host_depleted_reads,
        "mapped_reads":        mapped_reads,
        "host_fraction":       host_stats["host_fraction"] if host_stats else None,
        "min_qual":            8,    # pipeline default; not in scope to pass dynamically
        "min_length":          200,
        "max_length":          10000,
    }

    # ---- Per-branch data ----
    # Index each per-branch file by genotype label (extracted from filename).
    branches_to_run: list[str] = gt_summary.get("branches_to_run", [])

    def _index_files(paths: list[str] | None) -> dict[str, Path]:
        """Map genotype string -> Path for a list of file paths."""
        result: dict[str, Path] = {}
        if not paths:
            return result
        for p_str in paths:
            p = Path(p_str)
            if not p.exists():
                continue
            gt = _extract_genotype_from_filename(p.name, sample_id)
            if gt:
                result[gt] = p
        return result

    mosdepth_by_gt  = _index_files(args.mosdepth_summaries)
    mosdepth_bed_by_gt = _index_files(args.mosdepth_beds)
    variants_by_gt  = _index_files(args.variant_tsvs)
    flagstat_by_gt  = _index_files(args.flagstats)
    haplotype_by_gt = _index_files(args.haplotype_reports)

    # Collect all distinct genotype keys from all input files + branches_to_run
    all_gt_keys: set[str] = set(branches_to_run)
    for d in [mosdepth_by_gt, mosdepth_bed_by_gt, variants_by_gt, flagstat_by_gt, haplotype_by_gt]:
        all_gt_keys |= set(d.keys())

    branches: list[dict] = []
    for gt in sorted(all_gt_keys):
        # Coverage — summary for stats, BED for per-window chart.
        mosdepth_data = parse_mosdepth_summary(mosdepth_by_gt.get(gt))
        if mosdepth_data is None:
            mosdepth_data = {"mean": None, "min": None, "max": None,
                             "reference_length": None, "regions": []}
        bed_regions = parse_mosdepth_bed(mosdepth_bed_by_gt.get(gt))
        if bed_regions:
            mosdepth_data["svg"] = _make_coverage_svg(bed_regions)
        else:
            mosdepth_data["svg"] = ""

        # Variants
        variant_data = parse_variant_tsv(variants_by_gt.get(gt))

        # Flagstat
        flagstat_data = parse_flagstat(flagstat_by_gt.get(gt))

        # Haplotypes
        stitch_data = parse_haplotype_report(haplotype_by_gt.get(gt))

        coverage_mean  = mosdepth_data["mean"] if mosdepth_data else None
        variant_count  = len(variant_data) if variant_data is not None else None
        haplotype_count = stitch_data.get("chains_stitched") if stitch_data else None

        branches.append({
            "genotype":       gt,
            "top_subtype":    subtype_by_gt.get(gt, gt),
            "coverage":       mosdepth_data,
            "coverage_mean":  coverage_mean,
            "variants":       variant_data,
            "variant_count":  variant_count,
            "flagstat":       flagstat_data,
            "haplotypes":     stitch_data,
            "haplotype_count": haplotype_count,
        })

    # ---- NanoPlot images ----
    nanoplot_dirs: list[Path] = []
    if args.nanoplot_dirs:
        for d in args.nanoplot_dirs:
            p = Path(d)
            if p.is_dir():
                nanoplot_dirs.append(p)

    nanoplot_images = collect_nanoplot_pngs(nanoplot_dirs)

    # ---- Overall status ----
    overall_status = _determine_status(flags, is_mixed, primary_genotype)

    run_info: dict | None = None
    if args.run_info:
        try:
            run_info = json.loads(Path(args.run_info).read_text())
        except Exception as exc:
            print(f"WARNING: Cannot parse --run-info {args.run_info}: {exc}", file=sys.stderr)

    context = {
        "sample_id":         sample_id,
        "run_date":          run_date,
        "pipeline_version":  pipeline_version,
        "overall_status":    overall_status,
        "is_mixed":          is_mixed,
        "primary_genotype":  primary_genotype,
        "primary_subtype":   primary_subtype,
        "detected_subtypes": detected_subtypes,
        "genotype_summary":  gt_summary,
        "funnel":            funnel,
        "branches":          branches,
        "nanoplot_images":   nanoplot_images,
        "flags":             flags,
        "run_info":          run_info,
    }
    return context


def build_json_summary(context: dict) -> dict:
    """Build the machine-readable JSON summary (subset of context, serialisable)."""
    # Strip SVG and base64 image data — those are HTML-only.
    branches_json: list[dict] = []
    for b in context.get("branches", []):
        cov = b.get("coverage")
        cov_out: dict | None = None
        if cov is not None:
            cov_out = {k: v for k, v in cov.items() if k not in ("regions", "svg")}
        hap = b.get("haplotypes")
        hap_out: dict | None = None
        if hap is not None:
            hap_out = {
                k: v for k, v in hap.items()
                if k not in ("chains",)
            }
            # Include minimal chain summaries
            hap_out["chains"] = [
                {
                    "chain_id":    c.get("id") or c.get("chain_id"),
                    "abundance":   c.get("abundance_lower_bound") or c.get("abundance"),
                    "total_length": c.get("length_bp") or c.get("total_length"),
                    "n_windows":   len(c.get("windows") or []) or c.get("n_windows"),
                }
                for c in (hap.get("chains") or [])
            ]

        # variants: include count and first 20 by AF
        variants = b.get("variants") or []
        variants_sorted = sorted(variants, key=lambda v: v.get("af", 0), reverse=True)

        branches_json.append({
            "genotype":        b["genotype"],
            "top_subtype":     b.get("top_subtype"),
            "coverage_mean":   b.get("coverage_mean"),
            "coverage":        cov_out,
            "variant_count":   b.get("variant_count"),
            "variants_top20":  variants_sorted[:20],
            "haplotype_count": b.get("haplotype_count"),
            "haplotypes":      hap_out,
        })

    return {
        "sample_id":         context["sample_id"],
        "run_date":          context["run_date"],
        "pipeline_version":  context["pipeline_version"],
        "overall_status":    context["overall_status"],
        "is_mixed":          context["is_mixed"],
        "primary_genotype":  context["primary_genotype"],
        "primary_subtype":   context.get("primary_subtype"),
        "detected_subtypes": context.get("detected_subtypes", []),
        "genotype_summary":  context["genotype_summary"],
        "funnel":            context["funnel"],
        "branches":          branches_json,
        "flags":             context["flags"],
        "run_info":          context.get("run_info"),
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render per-sample HTML + JSON report for QuIVER.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--sample-id",          required=True)
    parser.add_argument("--genotype-summary",    default=None,
                        help="Path to *_genotype_summary.json")
    parser.add_argument("--nanoq-raw",           default=None,
                        help="Path to raw nanoq JSON")
    parser.add_argument("--nanoq-filtered",      default=None,
                        help="Path to post-filter nanoq JSON")
    parser.add_argument("--host-stats",          default=None,
                        help="Path to host_stats.json")
    parser.add_argument("--mosdepth-summaries",  nargs="*", default=None,
                        help="Per-branch *.mosdepth.summary.txt files")
    parser.add_argument("--mosdepth-beds",       nargs="*", default=None,
                        help="Per-branch *.regions.bed.gz files (per-window coverage chart)")
    parser.add_argument("--variant-tsvs",        nargs="*", default=None,
                        help="Per-branch *_variants.tsv files")
    parser.add_argument("--flagstats",           nargs="*", default=None,
                        help="Per-branch *_round2.flagstat files")
    parser.add_argument("--haplotype-reports",    nargs="*", default=None,
                        help="Per-branch *_haplotype_report.json files")
    parser.add_argument("--nanoplot-dirs",       nargs="*", default=None,
                        help="NanoPlot output directories (PNG images embedded)")
    parser.add_argument("--flags",               nargs="*", default=None,
                        help="Sentinel flag files (e.g. *.LOW_COVERAGE)")
    parser.add_argument("--output-html",         required=True)
    parser.add_argument("--output-json",         required=True)
    parser.add_argument("--run-info",            default=None,
                        help="JSON file with command_line, run_name, nextflow_version, params")
    parser.add_argument("--template",            default=None,
                        help="Path to sample_report.html.j2 (auto-detected if absent)")
    parser.add_argument("--pipeline-version",    default=None,
                        help="Pipeline version string for report header")
    return parser.parse_args()


def find_template(explicit: str | None) -> Path:
    """Locate the Jinja2 template file."""
    if explicit:
        p = Path(explicit)
        if p.is_file():
            return p
        sys.exit(f"ERROR: --template path not found: {explicit}")

    # Search relative to this script's location and common project layouts.
    script_dir = Path(__file__).resolve().parent
    candidates = [
        script_dir.parent / "assets" / "templates" / "sample_report.html.j2",
        script_dir / "sample_report.html.j2",
        Path("assets") / "templates" / "sample_report.html.j2",
    ]
    for c in candidates:
        if c.is_file():
            return c

    sys.exit(
        "ERROR: Cannot find sample_report.html.j2.  "
        "Supply the path explicitly with --template."
    )


# ---------------------------------------------------------------------------
# MultiQC custom content TSV generation
# ---------------------------------------------------------------------------

_SAMPLE_MQC_HEADER = """\
# id: 'quiver_sample_stats'
# section_name: 'QuIVER: Sample Summary'
# description: 'Per-sample read processing and genotyping statistics from the QuIVER HCV pipeline.'
# plot_type: 'generalstats'
# pconfig:
#   - status:
#       title: 'Status'
#       description: 'Overall pipeline status for this sample'
#   - genotype:
#       title: 'Genotype'
#       description: 'Primary HCV genotype detected'
#   - subtype:
#       title: 'Subtype'
#       description: 'Top HCV subtype (e.g. 1a, 3a)'
#   - is_mixed:
#       title: 'Mixed'
#       description: 'Mixed genotype infection (≥5% reads to secondary genotype)'
#   - raw_reads:
#       title: 'Raw Reads'
#       description: 'Total raw reads from sequencer'
#       format: '{:,.0f}'
#       min: 0
#       scale: 'Greens'
#   - filtered_reads:
#       title: 'Filtered Reads'
#       description: 'Reads after chopper length/quality filtering'
#       format: '{:,.0f}'
#       min: 0
#       scale: 'Blues'
#   - hcv_reads:
#       title: 'HCV Reads'
#       description: 'Reads remaining after host depletion'
#       format: '{:,.0f}'
#       min: 0
#       scale: 'Blues'
#   - host_pct:
#       title: 'Host %'
#       description: 'Percentage of reads removed as host (GRCh38)'
#       max: 100
#       min: 0
#       suffix: '%'
#       scale: 'Reds'
"""

_BRANCH_MQC_HEADER = """\
# id: 'quiver_branch_stats'
# section_name: 'QuIVER: Genotype Branch Summary'
# description: 'Per-genotype variant calling and haplotype reconstruction statistics from the QuIVER HCV pipeline. Rows are keyed {sample}_{genotype} to align with mosdepth/flagstat columns.'
# plot_type: 'generalstats'
# pconfig:
#   - variant_count:
#       title: 'Variants'
#       description: 'Filtered LoFreq variants (AF ≥ 1%, DP ≥ 20)'
#       format: '{:,.0f}'
#       min: 0
#       scale: 'Oranges'
#   - haplotype_count:
#       title: 'Haplotypes'
#       description: 'Reconstructed haplotypes from DEVIDER'
#       format: '{:,.0f}'
#       min: 0
#       scale: 'Purples'
#   - dominant_hap_pct:
#       title: 'Dom. Hap. %'
#       description: 'Abundance of the most dominant haplotype'
#       max: 100
#       min: 0
#       suffix: '%'
#       scale: 'RdYlGn'
"""


def _fmt(v: Any) -> str:
    """Return v as string, or empty string for None."""
    return "" if v is None else str(v)


def write_multiqc_tsv(json_out: dict, sample_id: str) -> None:
    """Write *_quiver_sample_mqc.tsv and *_quiver_branch_mqc.tsv for MultiQC."""
    funnel   = json_out.get("funnel") or {}
    status   = json_out.get("overall_status", "")
    is_mixed = json_out.get("is_mixed", False)
    branches = json_out.get("branches") or []

    primary_gt      = json_out.get("primary_genotype") or ""
    primary_subtype = json_out.get("primary_subtype") or ""
    # Fall back to genotype_summary entries if primary_subtype wasn't serialised
    if not primary_subtype and primary_gt:
        gt_entries = (json_out.get("genotype_summary") or {}).get("genotypes") or []
        for entry in gt_entries:
            if str(entry.get("genotype")) == str(primary_gt):
                primary_subtype = entry.get("top_subtype") or ""
                break

    # For mixed infections show all genotypes joined
    if is_mixed and branches:
        gt_display = "+".join(b["genotype"] for b in branches)
    else:
        gt_display = primary_gt or "N/A"

    host_fraction = funnel.get("host_fraction") or 0
    host_pct      = f"{host_fraction * 100:.2f}"

    # --- sample-level TSV ---
    sample_tsv = Path(f"{sample_id}_quiver_sample_mqc.tsv")
    with sample_tsv.open("w") as fh:
        fh.write(_SAMPLE_MQC_HEADER)
        fh.write("Sample\tstatus\tgenotype\tsubtype\tis_mixed\traw_reads\tfiltered_reads\thcv_reads\thost_pct\n")
        fh.write(
            f"{sample_id}\t{status}\t{gt_display}\t{primary_subtype or 'N/A'}\t"
            f"{'Yes' if is_mixed else 'No'}\t"
            f"{_fmt(funnel.get('raw_reads'))}\t"
            f"{_fmt(funnel.get('filtered_reads'))}\t"
            f"{_fmt(funnel.get('host_depleted_reads'))}\t"
            f"{host_pct}\n"
        )

    # --- branch-level TSV (one row per genotype branch) ---
    branch_tsv = Path(f"{sample_id}_quiver_branch_mqc.tsv")
    with branch_tsv.open("w") as fh:
        fh.write(_BRANCH_MQC_HEADER)
        fh.write("Sample\tvariant_count\thaplotype_count\tdominant_hap_pct\n")
        for b in branches:
            gt      = b.get("genotype", "")
            row_id  = f"{sample_id}_{gt}"
            vc      = _fmt(b.get("variant_count"))
            hc      = _fmt(b.get("haplotype_count"))
            # Dominant haplotype abundance from haplotype chains
            hap_data = b.get("haplotypes") or {}
            chains   = hap_data.get("chains") or []
            if chains:
                dom_abund = max(
                    (c.get("abundance") or c.get("abundance_lower_bound") or 0)
                    for c in chains
                )
                dom_pct = f"{dom_abund * 100:.1f}"
            else:
                dom_pct = ""
            fh.write(f"{row_id}\t{vc}\t{hc}\t{dom_pct}\n")


def main() -> None:
    args = parse_args()

    template_path = find_template(args.template)

    env = Environment(
        loader=FileSystemLoader(str(template_path.parent)),
        autoescape=select_autoescape(["html"]),
    )
    template = env.get_template(template_path.name)

    context  = build_context(args)
    html_out = template.render(**context)
    json_out = build_json_summary(context)

    Path(args.output_html).write_text(html_out, encoding="utf-8")
    Path(args.output_json).write_text(
        json.dumps(json_out, indent=2, default=str),
        encoding="utf-8",
    )
    write_multiqc_tsv(json_out, args.sample_id)

    print(
        f"[render_sample_report] Done: {args.output_html}  {args.output_json}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
