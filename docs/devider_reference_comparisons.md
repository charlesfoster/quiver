# DEVIDER reference and version comparisons

Run date: 2026-10-02. These experiments informed the dual-reference DEVIDER
implementation, with Round 2 simple as primary and Round 1 as the comparison.

## Question and design

QuIVER previously gave DEVIDER only the first sample-specific consensus. We
compared that reference with the Round 2 **simple** consensus, then reran
DEVIDER with version 0.1.0 on the same inputs used by version 0.0.1. The comparisons cover
the bundled synthetic `mixed_gt` subtype 1a test and two PromethION samples:
barcode16 (subtype 1b) and barcode61 (subtype 1a).

For each reference arm, the paired workflow started with the same subtype reads,
mapped them against that arm's FASTA, and made separate depth-capped BAMs and
LoFreq VCFs. Thus each DEVIDER BAM, VCF, and FASTA uses matching coordinates.
The original paired runs used bcftools 1.24, LoFreq 2.1.5, DEVIDER 0.0.1, and
Nextflow 24.10.5. The two real samples included the normal host depletion step;
the bundled synthetic test skipped host depletion.

The version 0.1.0 runs reused each arm's **exact preprocessed DEVIDER BAM,
DEVIDER-filtered VCF, and reference FASTA**. They did not rerun mapping,
subsampling, variant calling, filtering, or consensus construction. The pinned
image was `quay.io/biocontainers/devider:0.1.0--hb42e459_0` (verified by
`devider --version`; image digest
`sha256:6b6968fcb4ec7e05aff5ec69ce493c04567bea6dab9b198671c16332d22d61be`).
Each run used 4 threads, `--preset nanopore-r10`, `--min-cov 10`,
`--min-abund 0.25`, `--output-reads`, and `--allele-output`. The existing QuIVER
`format_haplotypes.py` script formatted the new `majority_vote_haplotypes.fasta`
files. Every DEVIDER run succeeded and produced one window.

**SNP positions** below means numeric marker rows in DEVIDER's `hap_info.txt`;
the column-header row is excluded. It is not the number of all filtered LoFreq
variants because DEVIDER uses a more selective phasing input and may omit
positions from its reported haplotypes.

## Fixed upstream comparison inputs

These counts are the same for both DEVIDER versions because the version 0.1.0
runs reused the original upstream files.

| Dataset | Reference | Reference length | Primary mapped reads | Filtered LoFreq records |
|---|---|---:|---:|---:|
| Synthetic `mixed_gt`, 1a | First | 9,416 | 464 | 1,910 |
| Synthetic `mixed_gt`, 1a | Round 2 simple | 9,416 | 465 | 1,909 |
| Barcode16, 1b | First | 9,413 | 442,158 | 628 |
| Barcode16, 1b | Round 2 simple | 9,413 | 442,158 | 623 |
| Barcode61, 1a | First | 9,286 | 328,901 | 881 |
| Barcode61, 1a | Round 2 simple | 9,285 | 329,068 | 884 |

The synthetic references differ at 26 bases. Barcode16's references differ at
six bases, five clustered around positions 1489–1503. Barcode61's references
have nine short alignment differences, including substitutions and indels, with
a net one-base length difference. Both real samples received QuIVER's low
mapping rate flag in both arms: 89.89% for barcode16 and 84.96%/85.01% for
barcode61 (first/simple).

## DEVIDER results

`N` median and maximum count ambiguous bases in each formatted full-length
haplotype. “Below 1%” counts haplotypes whose estimated abundance is below 1%;
DEVIDER itself used the same 0.25% minimum abundance setting in every run.

| Dataset | Reference | DEVIDER | SNP positions | Haplotypes | Length range (bp) | Median / max `N` | Top abundance | Below 1% |
|---|---|---|---:|---:|---:|---:|---:|---:|
| Synthetic 1a | First | 0.0.1 | 323 | 3 | 9,224–9,416 | 1 / 13 | 43.78% | 0 |
| Synthetic 1a | Simple | 0.0.1 | 323 | 3 | 9,224–9,416 | 1 / 13 | 43.72% | 0 |
| Synthetic 1a | First | 0.1.0 | 323 | 3 | 9,224–9,416 | 1 / 13 | 43.78% | 0 |
| Synthetic 1a | Simple | 0.1.0 | 323 | 3 | 9,224–9,416 | 1 / 13 | 43.72% | 0 |
| Barcode16 1b | First | 0.0.1 | 9 | 2 | 9,248–9,249 | 6 / 7 | 94.41% | 0 |
| Barcode16 1b | Simple | 0.0.1 | 4 | 6 | 9,245–9,249 | 8.5 / 15 | 47.09% | 0 |
| Barcode16 1b | First | 0.1.0 | 9 | 2 | 9,248–9,249 | 6 / 7 | 94.41% | 0 |
| Barcode16 1b | Simple | 0.1.0 | 4 | 6 | 9,245–9,249 | 8.5 / 15 | 47.09% | 0 |
| Barcode61 1a | First | 0.0.1 | 47 | 39 | 5,577–9,130 | 14 / 990 | 13.72% | 14 |
| Barcode61 1a | Simple | 0.0.1 | 44 | 19 | 9,128–9,129 | 4 / 61 | 24.20% | 0 |
| Barcode61 1a | First | 0.1.0 | 47 | 39 | 5,577–9,130 | 14 / 990 | 13.73% | 14 |
| Barcode61 1a | Simple | 0.1.0 | 44 | 21 | 9,128–9,129 | 4 / 1,954 | 24.45% | 0 |

### Synthetic test

Both references and both DEVIDER versions produced the **same three haplotype
sequences**. Their abundance estimates differed slightly by reference (first:
43.78%, 28.74%, 27.49%; simple: 43.72%, 28.77%, 27.51%). Changing DEVIDER
version did not change the sequences or the displayed abundances in either arm.
This is a small synthetic branch with only 464–465 primary mapped reads.

### Barcode16

Both DEVIDER versions gave two haplotypes with the first reference and six with
the simple reference. The dominant sequence from the first arm is also present
in the simple arm. The sequences within each arm were exactly the same between
DEVIDER versions; only a 0.01 percentage point rounding difference appeared in
the simple arm's top-five abundance sum. Changing the reference removed five
positions from DEVIDER's reported SNP set (9 to 4), including positions in the
short region where the references differ. The extra haplotypes therefore
cannot be interpreted as greater biological accuracy without independent
evidence.

### Barcode61

With the first reference, both DEVIDER versions gave 39 haplotypes, including
one 5,577 bp sequence with 990 `N` bases. Fourteen haplotypes were below 1%
estimated abundance (9.39% combined in 0.0.1; 9.38% in 0.1.0). Across the two
version runs, 32 of 39 full sequences matched exactly; the remaining sequence
differences show some run-to-run variation even when the input is unchanged.

With the simple reference, version 0.0.1 gave 19 haplotypes, all 9,128–9,129
bp long with at most 61 `N` bases. Version 0.1.0 gave 21 haplotypes, also all
9,128–9,129 bp long, but one 1.33% haplotype contained 1,954 `N` bases. A
second 0.1.0 DEVIDER run on the identical inputs again produced 21 haplotypes
and the same SNP haplotype sequences; the `N`-rich sequence recurred at 1.34%.
No haplotype was below 1% in either simple-reference run. The top five
haplotypes carried 68.35% of estimated abundance under 0.0.1 and 67.51% under
0.1.0, versus 51.52% and 51.53% with the first reference.

After aligning the two reference coordinate systems, all 44 simple-arm SNP
positions in `hap_info.txt` corresponded to positions in the first arm. The
first arm had three additional positions near the 3′ end. Exact FASTA string
equality between arms is not informative because their references differ in
length. An alignment of the version 0.1.0 haplotypes mapped all 21 simple-arm
sequences to first-arm sequences; the best alignments for the five most abundant
simple-arm sequences had 4–9 differences across roughly 9.1 kb aligned spans.

## Interpretation and limits

The Round 2 simple reference changes DEVIDER's output substantially in both
real samples, but in opposite directions: barcode16 goes from 2 to 6
haplotypes, while barcode61 goes from 39 to 19 or 21. Version 0.1.0 preserves
the first-reference counts and the barcode16 counts, but adds two barcode61
simple-reference haplotypes relative to 0.0.1. It also introduces an `N`-rich
barcode61 simple-reference haplotype that was absent from the 0.0.1 output.

Barcode61's simple-reference sequences are more consistently full length and
have fewer low-abundance calls, which explains why they can look “neater.”
These properties do not establish which reference recovers the true biological
haplotypes. The missing low-abundance calls could be noise or real diversity;
the additional barcode16 calls could be better resolution or over-splitting.
Neither sample has an independent truth set. The two samples do not support a
general switch of QuIVER's DEVIDER reference yet.

Each run produced one window, so cross-window linking was unnecessary. The
version 0.1.0 reruns executed the DEVIDER binary and the existing formatter;
they did not run QuIVER's separate `haplotag_bam` helper. Consequently, the
formatter reports `fallback_used: true` with reason “no haplotype-tagged BAM
available,” although it parsed and emitted the DEVIDER haplotypes normally.
The missing BAM does not affect the one-window haplotype count or sequences.
Temporary result directories may be cleared by the operating system.

