#!/usr/bin/env python3
"""Regression checks for dual-reference DEVIDER report handling."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

from jinja2 import Environment, FileSystemLoader, select_autoescape


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


sample_report = load_module("render_sample_report", ROOT / "bin/render_sample_report.py")
run_summary = load_module("render_run_summary", ROOT / "bin/render_run_summary.py")
format_haplotypes = load_module("format_haplotypes", ROOT / "bin/format_haplotypes.py")


def haplotype_report(arm: str, count: int, primary: bool) -> dict:
    label = "Round 2 simple consensus" if primary else "Round 1 consensus"
    return {
        "sample_id": "S01",
        "genotype": "1a",
        "reference_arm": arm,
        "reference_label": label,
        "primary": primary,
        "devider_marker_positions": 45 if primary else 48,
        "windows_found": 1,
        "chains_stitched": count,
        "fallback_used": False,
        "unlinked_windows": [],
        "chains": [
            {
                "id": f"S01_1a_{arm}_haplotype_{i}",
                "length_bp": 9129 if primary else 9130,
                "n_bases": 4 if primary else 14,
                "abundance_lower_bound": 1 / count,
                "windows": [0],
            }
            for i in range(count)
        ],
    }


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp_name:
        tmp = Path(tmp_name)

        chain = format_haplotypes.Chain(
            windows=[0],
            hap_indices=[0],
            sequence="ACGTNN",
            spanning_reads_per_junction=[],
            abundance_lower_bound=0.75,
        )
        formatted_fasta = tmp / "round1.fasta"
        format_haplotypes.write_fasta(
            formatted_fasta, [chain], "S01", "1a", id_suffix="_round1_comparison"
        )
        assert formatted_fasta.read_text().startswith(
            ">S01_1a_round1_comparison_haplotype_0"
        )
        formatted_report = format_haplotypes.build_report(
            sample_id="S01",
            genotype="1a",
            haplotypes=[],
            junctions=[],
            chains=[chain],
            unparsed_headers=[],
            fallback_used=False,
            fallback_reason=None,
            bam_present=False,
            min_reads=5,
            reference_arm="round1",
            reference_label="Round 1 consensus",
            primary=False,
            marker_positions=48,
            id_suffix="_round1_comparison",
        )
        assert formatted_report["chains"][0]["n_bases"] == 2
        assert formatted_report["devider_marker_positions"] == 48

        devider_dir = tmp / "devider_output"
        devider_dir.mkdir()
        (devider_dir / "majority_vote_haplotypes.fasta").write_text(
            ">Contig:test,Range:ALL-ALL,Haplotype:0,Abundance:75.0,Depth:100\n"
            "ACGTNN\n"
        )
        (devider_dir / "hap_info.txt").write_text(
            "Contig:test,Range:ALL-ALL\tHaplotype:0\n"
            "100\t0:0.99\n200\t1:0.95\n"
        )
        cli_fasta = tmp / "cli_haplotypes.fasta"
        cli_json = tmp / "cli_report.json"
        cli_mapping = tmp / "cli_mapping.tsv"
        assert format_haplotypes.main([
            "--devider-dir", str(devider_dir),
            "--output-fasta", str(cli_fasta),
            "--output-json", str(cli_json),
            "--output-mapping", str(cli_mapping),
            "--sample-id", "S01",
            "--genotype", "1a",
            "--reference-arm", "round2_simple",
            "--reference-label", "Round 2 simple consensus",
            "--primary",
        ]) == 0
        cli_report = json.loads(cli_json.read_text())
        assert cli_report["primary"] is True
        assert cli_report["devider_marker_positions"] == 2
        assert cli_report["chains"][0]["n_bases"] == 2

        genotype_summary = tmp / "S01_genotype_summary.json"
        genotype_summary.write_text(json.dumps({
            "primary_genotype": "1a",
            "total_mapped_reads": 100,
            "ambiguous_fraction": 0.0,
            "is_mixed": False,
            "branches_to_run": ["1a"],
            "genotypes": [{
                "genotype": "1a",
                "top_subtype": "1a",
                "reads": 100,
                "fraction": 1.0,
                "top_reference": "test_reference",
            }],
        }))

        round1 = tmp / "S01_1a_round1_comparison_haplotype_report.json"
        round2 = tmp / "S01_1a_haplotype_report.json"
        round1.write_text(json.dumps(haplotype_report("round1", 39, False)))
        round2.write_text(json.dumps(haplotype_report("round2_simple", 21, True)))

        standard_flagstat = tmp / "S01_1a_round2.flagstat"
        round2_flagstat = tmp / "S01_1a_round2_simple_devider.flagstat"
        standard_flagstat.write_text("100 + 0 primary mapped (90.00% : N/A)\n")
        round2_flagstat.write_text("101 + 0 primary mapped (91.00% : N/A)\n")

        args = SimpleNamespace(
            sample_id="S01",
            genotype_summary=str(genotype_summary),
            nanoq_raw=None,
            nanoq_filtered=None,
            host_stats=None,
            mosdepth_summaries=None,
            mosdepth_beds=None,
            variant_tsvs=None,
            flagstats=[str(standard_flagstat), str(round2_flagstat)],
            haplotype_reports=[str(round1), str(round2)],
            nanoplot_dirs=None,
            flags=None,
            run_info=None,
            pipeline_version="test",
        )

        context = sample_report.build_context(args)
        branch = context["branches"][0]
        assert branch["primary_haplotype_reference"] == "round2_simple"
        assert branch["haplotype_count"] == 21
        assert branch["haplotypes_by_reference"]["round1"]["chains_stitched"] == 39
        assert branch["haplotypes_by_reference"]["round2_simple"]["mapping"]["primary_mapped"] == 101

        summary = sample_report.build_json_summary(context)
        assert summary["primary_haplotype_reference"] == "round2_simple"
        summary_branch = summary["branches"][0]
        assert summary_branch["haplotype_count"] == 21
        assert set(summary_branch["haplotypes_by_reference"]) == {"round1", "round2_simple"}
        assert run_summary.extract_sample_row(summary)["haplotypes_reconstructed"] == 21

        template_path = ROOT / "assets/templates/sample_report.html.j2"
        env = Environment(
            loader=FileSystemLoader(str(template_path.parent)),
            autoescape=select_autoescape(["html"]),
        )
        html = env.get_template(template_path.name).render(**context)
        assert "Round 2 simple consensus" in html
        assert "Round 1 comparison haplotypes (39)" in html
        assert "Round 2 simple consensus is the primary result" in html


if __name__ == "__main__":
    main()
