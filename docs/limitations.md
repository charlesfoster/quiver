# Known Limitations and Future Work — HCV Quasispecies Pipeline

1. **Genotype determination by best-hit only.** Recombinants assigned to a single genotype; intra-host recombination not explicitly detected. Future: RDP4-style breakpoint analysis on assembled haplotypes.

2. **Single-genome reference per genotype branch.** Dominant panel reference used as scaffold — may introduce small bias artefacts for strains phylogenetically between two panel sequences. Future: multi-reference graph genome (vg, minigraph) for consensus step.

3. **No GPU acceleration.** Clair3 has GPU paths; not configured here. Future: add GPU profile.

4. **No drug-resistance annotation.** Variants reported as raw VCFs; no mapping to HCV DAA resistance positions (NS3, NS5A, NS5B). Future: integrate HCV-GLUE or custom annotation layer.

5. **Host depletion uses one human reference.** Non-human, non-viral reads pass through. Acceptable for sterile clinical samples. Future: optional multi-host depletion.

6. **DEVIDER haplotype stitching is heuristic.** Hamming-distance support thresholds; not probabilistic. Future: probabilistic stitcher modelling sequencing error.

7. **No phylogenetic placement.** Sequences produced but not placed on HCV phylogeny. Future: Nextclade-style placement.

8. **No support for paired/replicate samples.** Each sample is independent. Future: cross-sample variant rescue and contamination detection.

9. **No structural variant calling.** Large deletions/insertions not specifically detected beyond what LoFreq emits as long indels.

10. **ARM64 Mac production runs slower than x86 HPC** due to Rosetta translation. Native arm64 builds for DEVIDER and x86-only biocontainers would be needed. Conda-via-pixi local profile is the recommended escape hatch.

11. **LoFreq strand-bias filter was Illumina-tuned.** ONT R10.4.1 produces both strands reliably but with subtly different error profiles. Mitigated by lowering `--min-bq` and `alnqual` recalibration; advanced users may wish to re-tune `--sig` and the strand-bias filter.

12. **No real-time monitoring.** Pipeline is batch only. Future: Epi2Me-style streaming variant.
