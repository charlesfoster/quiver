# Consensus genomes produced by QuIVER

QuIVER emits three consensus FASTAs per genotype branch, not necessarily only
three per sample. A mixed-genotype sample gets a separate set under each
`consensus/<GT>/` directory.

| Output | Meaning | Main use |
|---|---|---|
| `<sample>_<GT>_consensus.fasta` | First-pass, sample-specific majority consensus | Pipeline reference and coordinate system |
| `<sample>_<GT>_round2_consensus_simple.fasta` | Second-pass majority consensus | Final single representative genome |
| `<sample>_<GT>_round2_consensus_iupac.fasta` | Second-pass ambiguity-coded population consensus | Showing within-sample variation |
| `<sample>_<GT>_haplotypes.fasta` | Multiple reconstructed, phased genomes | Quasispecies and haplotype analysis |

## 1. Initial sample consensus

**File:** `<sample>_<GT>_consensus.fasta`

This consensus is constructed by:

1. Selecting the closest or dominant reference from the HCV panel.
2. Mapping that genotype's reads to the selected reference.
3. Calling variants with `bcftools` in haploid mode.
4. Applying the called alternate alleles to the panel reference.
5. Masking positions below `min_consensus_cov`, which defaults to 10x, with
   `N`.

The result contains ordinary A/C/G/T bases and Ns, without ambiguity codes. The
haploid caller is intended to select the dominant allele at variable positions;
only calls passing the depth gate are applied. See
[`bcftools_consensus_call/main.nf`](../modules/local/bcftools_consensus_call/main.nf)
and [`apply_consensus/main.nf`](../modules/local/apply_consensus/main.nf).

This is the operationally important consensus:

- Reads are remapped against it in Round 2.
- LoFreq calls variants relative to it.
- The Round 1 DEVIDER comparison reconstructs haplotypes relative to it.
- Reported VCF positions therefore use this sequence's coordinate system.

A `LOW_COVERAGE_CONSENSUS` flag is emitted if at least 30% of the sequence is
N-masked. A `DIVERGENT_CONSENSUS` flag indicates less than the configured
`min_consensus_identity`, default 90%, identity to the selected panel reference.
Both are warnings; the workflow continues.

## 2. Round 2 simple consensus

**File:** `<sample>_<GT>_round2_consensus_simple.fasta`

After remapping to the initial consensus, the pipeline uses LoFreq to call
variants more sensitively. By default, the LoFreq calls used here must meet:

- allele frequency at least 1%;
- depth at least 20;
- `FILTER=PASS`; and
- strand-bias score below 200.

For the simple consensus, only variants with allele frequency at least 50% are
applied. Positions with Round 2 coverage below 10x are again masked with `N`.
See [`variant_filter/main.nf`](../modules/local/variant_filter/main.nf) and
[`build_round2_consensus/main.nf`](../modules/local/build_round2_consensus/main.nf).

Conceptually, this is a second polishing pass:

```text
panel reference
   -> first-pass majority calling
initial sample consensus
   -> remapping and LoFreq majority corrections
round2 simple consensus
```

For most downstream uses requiring one representative sequence per sample, such
as phylogeny, genotyping, submission, or comparison between time points, this is
the most final single-genome product.

This FASTA is also the primary DEVIDER reference. The pipeline remaps the
genotype-assigned reads to it, makes a new coordinate-matched LoFreq VCF and
DEVIDER input BAM, and runs DEVIDER independently from the Round 1 comparison
arm. See [`workflows/quiver.nf`](../workflows/quiver.nf).

## 3. Round 2 IUPAC consensus

**File:** `<sample>_<GT>_round2_consensus_iupac.fasta`

This starts from the same initial consensus and the same filtered LoFreq calls,
but handles variants differently:

- At allele frequency at least 50%, the alternate allele is applied as an
  ordinary base.
- Below 50%, the variant is assigned a synthetic heterozygous genotype and
  `bcftools consensus -H I` emits an IUPAC ambiguity code representing the
  reference and alternate alleles.

For example, a reference `C` with a minority `T` variant becomes `Y`. The
implementation is in
[`build_round2_consensus/main.nf`](../modules/local/build_round2_consensus/main.nf).

Important interpretation details:

- By default, any passing SNV down to 1% allele frequency may introduce an
  ambiguity code. This is more permissive than using IUPAC only for an
  intermediate-frequency range such as 20-80%.
- Treatment is asymmetric at 50%: a 40% alternate allele becomes ambiguous,
  while a 60% alternate allele becomes pure alternate even though 40% of the
  reads still support the reference allele.
- LoFreq calls indels as well as SNVs, and those records remain in the input
  VCF. IUPAC has no natural representation for an indel. Consequently, this
  FASTA can change length and should not be understood as merely the simple
  consensus with ambiguous letters added.

This version is useful for visualising population-level diversity, but it is
less suitable for tools that require an unambiguous haploid genome.

## 4. Reconstructed haplotypes

**Primary file:** `haplotypes/<GT>/<sample>_<GT>_haplotypes.fasta`

**Comparison file:**
`haplotypes/<GT>/<sample>_<GT>_round1_comparison_haplotypes.fasta`

These are not consensus variants in the same sense. They are multiple phased
genomes reconstructed by DEVIDER, with headers such as:

```text
>sample_GT_haplotype_0 abund:31.98 depth:1392.79 length:9286
```

They aim to represent distinct viral lineages within the sample rather than
collapse the population into one sequence. They are ordered by estimated
abundance. The unqualified primary file comes from the Round 2 simple consensus
arm. The explicitly labelled comparison file comes from the initial consensus
arm. Each arm uses a FASTA, BAM, and VCF in its own matching coordinate system.

The formatter preferentially uses DEVIDER's full-length
`majority_vote_haplotypes.fasta`. The raw `snp_haplotypes.fasta` contains only
variable-site sequences and is not a complete genome. See
[`format_haplotypes.py`](../bin/format_haplotypes.py) and
[`format_haplotypes/main.nf`](../modules/local/format_haplotypes/main.nf).

Use these sequences when the biological question concerns quasispecies
composition, linkage, or minority lineages, rather than when a single canonical
sample sequence is required.

## Which output should be used?

- For one final genome per sample and genotype, use
  `round2_consensus_simple.fasta`.
- To preserve ambiguity from passing low-frequency SNVs, use
  `round2_consensus_iupac.fasta`, with awareness of its default 1% threshold and
  indel behaviour.
- To work with the main reported VCF coordinates or reproduce the standard
  Round 2 mapping, use the initial `consensus.fasta`.
- To study distinct within-host viral populations, use the primary formatted
  `haplotypes.fasta`; use the `round1_comparison` FASTA to assess reference
  sensitivity.

