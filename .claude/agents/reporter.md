---
name: reporter
description: Implement the per-sample HTML report generator and the MultiQC run summary. Covers Prompts 22–23 from docs/implementation_prompts.md. Run after all upstream modules are implemented.
model: claude-sonnet-4-6
---

You are implementing **Prompts 22–23** of the `hcv-quasi` pipeline: per-sample HTML reporting and MultiQC aggregation.

## Context
Read before starting:
- `docs/data_flow.md` — Steps 5.21–5.22
- `docs/implementation_prompts.md` — Prompts 22–23 gotchas
- `CLAUDE.md` — Section 4 (container tags); Section 1 (output structure)

## Files to create

### `bin/render_sample_report.py`
Python 3.11 script. Takes a directory of per-sample outputs and renders an HTML summary.

**Invocation:** `python render_sample_report.py <sample_dir> <output_html> <output_json>`

**What to collect (all paths relative to the sample results dir):**

| Data | Source file |
|---|---|
| Read counts at each stage | `qc/raw/nanoq.json`, `qc/posthost/nanoq.json`, chopper log |
| Host depletion % | `genotyping/host_stats.json` |
| Genotype summary | `genotyping/genotype_summary.json` |
| Per-branch coverage | `qc/coverage/${GT}/mosdepth.summary.txt` |
| Variants table | `variants/${GT}/variants.tsv` |
| Haplotype summary | `haplotypes/${GT}/stitching_report.json` (or absent if DEVIDER failed) |
| Pipeline flags | `qc/flags.txt` |
| NanoPlot images | `qc/raw/${SID}/NanoPlot-*.png` (embed as base64) |

**HTML structure:**
1. Header: sample ID, run date, pipeline version
2. Status banner: PASS / MIXED / NO_HCV / LOW_COVERAGE (colour-coded)
3. Read processing funnel table: raw → filtered → host-depleted → mapped (counts + % retained)
4. Genotype section: primary genotype, is_mixed flag, fraction table for all detected genotypes
5. Per-genotype tabs (one tab per branch):
   - Coverage plot (mosdepth data as a simple inline SVG bar chart)
   - Variants table (sortable HTML table from `variants.tsv`): POS, REF, ALT, AF, DP, SB flag
   - Haplotype summary: count, AF distribution, longest contiguous haplotype, stitching notes
6. NanoPlot read length/quality distribution (embedded PNG)
7. Flags section: list any emitted sentinel flags with brief explanations

Use Jinja2 for templating. Keep all CSS inline (no external dependencies) so the HTML file is self-contained.

### `assets/templates/sample_report.html.j2`
Jinja2 template. Referenced by `render_sample_report.py`. Variables available in template context match the section structure above. Use a clean, readable style — no heavy frameworks.

### `modules/local/sample_report.nf`
Process wrapping `bin/render_sample_report.py`.
- Input: `[meta, sample_results_dir]` (all per-sample outputs collected as a directory)
- Output: `${meta.id}_summary.html`, `${meta.id}_summary.json`
- Container: `python:3.11-slim` (with Jinja2 installed: `pip install jinja2`)
- Label: `process_low`

### `bin/render_run_summary.py`
Python 3.11 script. Aggregates all per-sample JSON summaries into a run-level dashboard.

**Invocation:** `python render_run_summary.py <results_dir> <output_html> <output_json>`

**Content:**
1. Run summary table: one row per sample with columns: `sample_id`, `status`, `primary_genotype`, `is_mixed`, `mean_coverage`, `variants_called`, `haplotypes_reconstructed`, link to per-sample report
2. Status colour codes: green (PASS), orange (MIXED), red (NO_HCV / LOW_COVERAGE)
3. Cross-sample statistics: total samples, pass rate, mean/median coverage across samples

### `modules/local/multiqc.nf`
```bash
multiqc -f \
    --title "hcv-quasi run summary" \
    -o multiqc_out \
    ${qc_dirs.join(' ')}
```

Input: list of all QC directories from all samples (NanoPlot outputs, mosdepth outputs, samtools flagstat outputs).

MultiQC parsers needed (must all be present in the collected directories):
- NanoPlot: `NanoStats.txt` (tsv_stats output)
- mosdepth: `*.mosdepth.summary.txt`
- samtools flagstat: `*.flagstat`
- bcftools stats: `*.bcftools_stats.txt`

Container: multiqc 1.25.1
Label: `process_medium`

## Success criteria
- `${SID}_summary.html` is a self-contained HTML file that opens in a browser without external requests
- Status banner correctly reflects sample flags (e.g. MIXED when `is_mixed = true`)
- Variants table is populated when variants exist; shows "No variants above threshold" when empty
- For DEVIDER-failed samples: haplotype section shows "Haplotype reconstruction failed" — not an error
- `run_summary.html` lists all samples with correct status colours
- `multiqc_report.html` is produced and contains NanoPlot and mosdepth sections
