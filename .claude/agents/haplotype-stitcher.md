---
name: haplotype-stitcher
description: Design and implement the cross-window haplotype stitching algorithm — walks DEVIDER's per-window haplotype output and the haplotype-tagged BAM to find read-spanning evidence linking adjacent haplotypes into longer sequences. Covers Prompt 21 from docs/implementation_prompts.md. Core algorithmic logic requires careful reasoning.
model: claude-opus-4-7
---

You are implementing **Prompt 21** of the `hcv-quasi` pipeline: post-hoc haplotype stitching.

## Context
Read before starting:
- `docs/data_flow.md` — Step 5.20 (logic description)
- `docs/architecture_reasoning.md` — Section 9 (why windowed DEVIDER + post-hoc stitching)
- `docs/implementation_prompts.md` — Prompt 21 gotchas
- `CLAUDE.md` — Section 3 D11; Section 10 item 6 (stitching is heuristic — document it)

## Background: DEVIDER's output structure

DEVIDER v0.0.1 operates in genomic windows. For each window where it can phase SNPs, it produces:
- One or more haplotype sequences (FASTA records) labelled with the window coordinates and a haplotype index
- When `--output-reads` is used: a BAM file where each read is tagged with its assigned haplotype ID

The exact DEVIDER output directory structure may vary — **do not hard-code filenames**. Instead, discover the structure by:
1. `glob("devider_out/**/*.fasta")` — find all FASTA files (haplotype sequences)
2. `glob("devider_out/**/*.bam")` — find the haplotype-tagged BAM(s)

DEVIDER haplotype FASTA record headers encode the window and haplotype. Typical format (v0.0.1): `>WINDOW_START-WINDOW_END_HAPLOTYPE_INDEX` or `>region:START-END|hap:N`. Parse defensively — check what is actually present.

## Files to create

### `bin/stitch_haplotypes.py`

**Invocation:**
```
python stitch_haplotypes.py \
    --devider-dir devider_out/ \
    --consensus consensus.fasta \
    --min-reads 5 \
    --sample-id P001 \
    --genotype 1a \
    --out-fasta merged_haplotypes.fasta \
    --out-json stitching_report.json
```

#### Algorithm overview

**Phase 1 — Parse DEVIDER output**

Discover and parse all haplotype FASTAs in `devider_out/`. For each haplotype sequence, record:
- Window start and end coordinates (genomic positions on the consensus)
- Haplotype index within the window
- The nucleotide sequence
- SNV alleles (if `--allele-output` was used, alleles are in the sequence rather than 0/1 codes)

Sort windows by start position. Group haplotypes by window.

**Phase 2 — Parse the haplotype-tagged BAM**

Open the BAM from `--output-reads`. For each read:
- Extract the haplotype tag: look for tag `HP` (standard haplotype tag) or `YH` (DEVIDER-specific). Parse to get `(window_id, haplotype_index)`.
- Record: `{read_id: [(window_id, hap_idx), ...]}` — one read may span multiple windows.

**Phase 3 — Build junction evidence**

For each pair of adjacent windows `(W_i, W_{i+1})`:
1. Find reads that have assignments in BOTH `W_i` AND `W_{i+1}` (spanning reads).
2. For each spanning read: record the `(hap_idx_in_W_i, hap_idx_in_{i+1})` pair.
3. Count how many reads support each pairing.
4. A link `(hap_A_in_W_i) → (hap_B_in_W_{i+1})` is **supported** if `count >= params.stitch_min_reads`.

When there are no spanning reads for a junction, no link is formed — emit the haplotypes separately.

**Phase 4 — Greedy path assembly**

For each window `W_0` haplotype, attempt to extend greedily:
1. From `hap_A` in `W_0`, follow the strongest supported link to `hap_B` in `W_1`.
2. Continue through subsequent windows until the chain breaks (no supported link).
3. A chain of length ≥ 2 windows produces a merged haplotype. A chain of length 1 is emitted as-is (DEVIDER's native resolution).

When a junction has multiple competing links (e.g., hap A→X and hap A→Y both have ≥ `min_reads` support), emit both chains as separate merged haplotypes rather than picking one arbitrarily. Document this in the stitching report.

**Phase 5 — Sequence concatenation**

For each assembled chain, concatenate the haplotype sequences from consecutive windows. At junctions, use the consensus sequence to fill any gap between windows (if `window_end_i < window_start_{i+1}`). If windows overlap, trim to avoid duplication at the junction (use the midpoint of the overlap as the boundary).

**Abundance estimation:**
Report the fraction of spanning reads that support the merged haplotype path, relative to total spanning reads at each junction. Approximate abundance = min of support fractions across all junctions in the chain.

**Output FASTA:** One record per distinct chain. Header format:
`>{sample_id}_{genotype}_hap{N}_{start}-{end}_abund{abundance:.3f}`

**Output JSON (`stitching_report.json`):**
```json
{
  "sample_id": "P001",
  "genotype": "1a",
  "n_windows": 3,
  "n_input_haplotypes": 7,
  "n_merged_haplotypes": 3,
  "junctions": [
    {
      "window_pair": ["W0", "W1"],
      "spanning_reads": 42,
      "links": [
        {"from": "W0_hap0", "to": "W1_hap0", "reads": 28, "supported": true},
        {"from": "W0_hap1", "to": "W1_hap1", "reads": 14, "supported": true}
      ]
    }
  ],
  "merged_haplotypes": [
    {"id": "hap0", "windows": ["W0", "W1", "W2"], "length": 9420, "abundance": 0.67},
    {"id": "hap1", "windows": ["W0", "W1"], "length": 6210, "abundance": 0.33}
  ],
  "unstitched_haplotypes": 1
}
```

**`devider.failed` handling:**
If `devider_out/devider.failed` exists, emit an empty `merged_haplotypes.fasta` and a minimal JSON `{"status": "devider_failed", ...}`. Do NOT raise an exception.

### `modules/local/stitch_haplotypes.nf`
```groovy
process STITCH_HAPLOTYPES {
    label 'process_medium'

    input:
    tuple val(meta), path(devider_dir), path(consensus_fasta)

    output:
    tuple val(meta), path("${meta.id}.${meta.genotype}.merged_haplotypes.fasta"), emit: haplotypes
    tuple val(meta), path("${meta.id}.${meta.genotype}.stitching_report.json"), emit: report

    script:
    """
    python ${projectDir}/bin/stitch_haplotypes.py \\
        --devider-dir ${devider_dir} \\
        --consensus ${consensus_fasta} \\
        --min-reads ${params.stitch_min_reads} \\
        --sample-id ${meta.id} \\
        --genotype ${meta.genotype} \\
        --out-fasta ${meta.id}.${meta.genotype}.merged_haplotypes.fasta \\
        --out-json ${meta.id}.${meta.genotype}.stitching_report.json
    """
}
```

## Documentation requirement
This is bespoke algorithmic code. Add a module-level docstring to `bin/stitch_haplotypes.py` explaining:
1. Why window-stitching is needed (DEVIDER v0.0.1 has no built-in cross-region merging)
2. What "spanning read" means in this context
3. The greedy assembly approach and its limitations (the abundance estimate is a lower bound; the algorithm cannot handle ambiguous branching paths optimally — it enumerates all supported branches)
4. What `stitch_min_reads` controls and how to tune it

## Success criteria
- Script completes (exit 0) even when DEVIDER failed (marker file present)
- For DEVIDER output with 2+ windows and sufficient spanning reads: at least one merged haplotype spans >1 window
- For DEVIDER output with 1 window: all haplotypes emitted as-is, stitching report notes `n_windows = 1`
- JSON report `n_input_haplotypes` equals the total number of FASTA records found in `devider_out/`
- Concatenated sequence length is consistent with the window coordinates (±10 bp tolerance for boundary handling)
- No hard-coded DEVIDER output filenames — discovery is by glob
