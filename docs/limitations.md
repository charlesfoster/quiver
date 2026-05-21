# Known Limitations and Future Work — QuIVER

1. **Genotype determination by best-hit only.** Recombinants assigned to a single genotype; intra-host recombination not explicitly detected. Future: RDP4-style breakpoint analysis on assembled haplotypes.

2. **Single-genome reference per genotype branch.** Dominant panel reference used as scaffold — may introduce small bias artefacts for strains phylogenetically between two panel sequences. Future: multi-reference graph genome (vg, minigraph) for consensus step.

3. **No GPU acceleration.** Clair3 has GPU paths; not configured here. Future: add GPU profile.

4. **No drug-resistance annotation.** Variants reported as raw VCFs; no mapping to HCV DAA resistance positions (NS3, NS5A, NS5B). Future: integrate HCV-GLUE or custom annotation layer.

5. **Host depletion uses one human reference.** Non-human, non-viral reads pass through. Acceptable for sterile clinical samples. Future: optional multi-host depletion.

6. **minimap2 host depletion requires a local reference.** The default nohuman mode auto-downloads its Kraken2 database, but `--use_minimap2` requires a local GRCh38 FASTA or pre-built `.mmi`. Passing a persistent `.mmi` is the fastest repeat-run path with minimap2. Future: add fully automatic GRCh38 download with a `storeDir`-backed index cache.

7. **DEVIDER haplotype stitching is heuristic.** Hamming-distance support thresholds; not probabilistic. Future: probabilistic stitcher modelling sequencing error.

8. **No phylogenetic placement.** Sequences produced but not placed on HCV phylogeny. Future: Nextclade-style placement.

9. **No support for paired/replicate samples.** Each sample is independent. Future: cross-sample variant rescue and contamination detection.

10. **No structural variant calling.** Large deletions/insertions not specifically detected beyond what LoFreq emits as long indels.

11. **ARM64 Mac production runs slower than x86 HPC** due to Rosetta translation. Native arm64 builds for DEVIDER and x86-only biocontainers would be needed. Conda-via-pixi local profile is the recommended escape hatch.

16. **Docker Desktop on macOS cannot access files inside `~/Documents`, `~/Desktop`, or `~/Downloads` without Full Disk Access.** macOS enforces TCC (Transparency, Consent, and Control) privacy restrictions on these folders. Docker Desktop's VirtioFS layer propagates those restrictions into containers — files under `Documents` are visible in directory listings but fail with `Operation not permitted` when accessed (`stat`, `open`, `cp`). This also causes a secondary failure during container init if the Nextflow work directory is itself inside `Documents`:

    ```
    runc create failed: error during container init: mkdir /Users/<you>/Documents: file exists
    ```

    **Root fix — grant Docker Desktop Full Disk Access:**

    > System Settings → Privacy & Security → Full Disk Access → `+` → add `/Applications/Docker.app` → restart Docker Desktop.

    This is required whenever the pipeline directory, input data, or assets live inside a TCC-protected folder. After granting access, all Docker-based profiles work normally.

    **Workaround (if you cannot grant Full Disk Access):** keep the pipeline and all input data outside `~/Documents`, `~/Desktop`, and `~/Downloads` (e.g. directly in `~/H2Seq/`), and always set the Nextflow work directory outside those folders:

    ```sh
    export NXF_WORK="$HOME/nf-work"
    ```

12. **LoFreq strand-bias filter was Illumina-tuned.** ONT R10.4.1 produces both strands reliably but with subtly different error profiles. Mitigated by lowering `--min-bq` and adding indel qualities with `lofreq indelqual`; advanced users may wish to re-tune `--sig` and the strand-bias filter.

13. **No real-time monitoring.** Pipeline is batch only. Future: Epi2Me-style streaming variant.

14. **Mixed-genotype test data is synthetic.** The current `mixed_gt` test case uses small synthetic reads and does not reflect the read depth or complexity of a real mixed-genotype clinical sample. Future: replace with simulation-derived reads spanning two or more major genotypes at clinically relevant proportions.

15. **Within-genotype subtype co-infections (e.g. 1a/1b) are not separately analysed.** The pipeline branches at the major-genotype level (1, 2, 3…), so genotype-1a and genotype-1b reads are merged into a single "genotype 1" branch. The dominant subtype is recorded in `top_subtype` but the minority subtype receives no independent variant calling or haplotype reconstruction. Future: add subtype-level branching to handle 1a/1b co-infections, with per-subtype output directories (e.g. `mapping/1a/`, `mapping/1b/`) and separate variant and haplotype outputs for each.
