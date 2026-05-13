/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BUILD_CONSENSUS_FASTA — Apply variants and low-coverage masking to produce
    the sample-specific consensus FASTA, FASTA index, and minimap2 index.

    Purpose:
        Step 5.10 sub-steps 4–5.  Generates the polished consensus that is
        the basis for Round 2 mapping and all downstream variant/haplotype
        work.  Two processes are used to avoid requiring a single container
        with mosdepth + bcftools + samtools + minimap2 all at once.

    ----------------------------------------------------------------
    MAKE_MASK_BED — Coverage mask generation (mosdepth container)
    ----------------------------------------------------------------
    Runs mosdepth in quantize mode to produce a BED of regions below
    params.min_consensus_cov.  These positions will be written as N
    in the final consensus FASTA (design decision D8).

    mosdepth --quantize 0:<threshold>:
        Produces a file with rows like:
            chr1  0    100  0:10      <- 0–10 coverage (below threshold)
            chr1  100  200  10:500    <- 10–500 coverage (above threshold)
        We keep rows where the label is "0:<threshold>" (the below-threshold bin).
        With threshold=10 (default) the below-threshold label is "0:10".

    The quantize string "0:<threshold>:" uses three colon-separated breakpoints
    that divide coverage into three bins:
        [0, 1)          — zero coverage
        [1, <threshold>) — below threshold (inadequate for consensus)
        [<threshold>, ∞) — adequate coverage
    We mask positions in both the zero-coverage and below-threshold bins,
    i.e. any bin label that does NOT end with ":<threshold>" as its left bound.

    Approach: extract rows whose 4th column is NOT the highest bin
    (i.e. label starts with "0:" or is the below-threshold range).
    The awk pattern matches labels where the numeric prefix is < threshold.

    ----------------------------------------------------------------
    APPLY_CONSENSUS — Apply variants + mask, rename header, index
    ----------------------------------------------------------------
    Uses the mulled bcftools+minimap2+samtools image. Steps:
        1. bcftools consensus: apply filtered VCF variants and mask.
        2. Python one-liner: rename the FASTA header to the stable
           identifier "${meta.id}_${meta.genotype}_consensus".
           (Python used instead of sed -i to avoid BSD vs GNU sed
           portability issues — see Prompt 11 gotchas.)
        3. samtools faidx: build .fai index.
        4. minimap2 -x map-ont -d: build .mmi index for Round 2.
        5. Pass/fail check: warn if length outside [4000,11000] bp;
           emit LOW_COVERAGE_CONSENSUS sentinel if N-fraction >= 30%.

    Container split:
        MAKE_MASK_BED — mosdepth 0.3.10
        APPLY_CONSENSUS — mulled bcftools 1.21 + minimap2 2.28 + samtools 1.21

    Labels:
        MAKE_MASK_BED  — process_low
        APPLY_CONSENSUS — process_medium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ============================================================================
// Process 1: Generate the low-coverage mask BED using mosdepth.
// ============================================================================

process MAKE_MASK_BED {

    label 'process_low'

    tag "${meta.id}:${meta.genotype}"

    container 'quay.io/biocontainers/mosdepth:0.3.10--h4e814b3_0'
    conda 'bioconda::mosdepth=0.3.10'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("mask.bed"), emit: mask_bed
    path "versions.yml",               emit: versions

    script:
    // Build the quantize string: "0:<threshold>:" — three breakpoints that
    // create three bins: [0,1), [1,threshold), [threshold,inf).
    // We want to mask any position in the first two bins.
    def threshold = params.min_consensus_cov
    """
    # ----------------------------------------------------------------
    # Run mosdepth in quantize mode.
    # --no-per-base suppresses the per-base .bed.gz (not needed here).
    # The quantize string "0:${threshold}:" creates bins:
    #   [0, ${threshold})  — low coverage (mask these)
    #   [${threshold}, ∞)  — adequate coverage
    # Note: mosdepth interprets the string as right-exclusive breakpoints
    # so "0:${threshold}:" means: 0 ≤ cov < ${threshold} and ${threshold} ≤ cov.
    # ----------------------------------------------------------------
    mosdepth \\
        --quantize 0:${threshold}: \\
        --no-per-base \\
        ${meta.id}_${meta.genotype}_cov \\
        ${bam}

    # ----------------------------------------------------------------
    # Extract positions in the below-threshold bin.
    # The quantized BED has column 4 = bin label, e.g. "0:10" or "10:500".
    # We keep rows whose label indicates coverage < threshold.
    # With "0:<threshold>:" the below-threshold label is "0:<threshold>".
    # ----------------------------------------------------------------
    LOW_BIN="0:${threshold}"

    zcat ${meta.id}_${meta.genotype}_cov.quantized.bed.gz \\
        | awk -v low_bin="\${LOW_BIN}" '\$4 == low_bin { print \$1"\\t"\$2"\\t"\$3 }' \\
        > mask.bed

    # mask.bed may be empty — that is fine; bcftools consensus accepts an empty -m.
    echo "Mask BED generated: \$(wc -l < mask.bed) regions below coverage ${threshold}." >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """

    stub:
    """
    touch mask.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: "0.3.10"
    END_VERSIONS
    """
}

// ============================================================================
// Process 2: Apply variants and mask; build FASTA + .fai + .mmi.
// ============================================================================

process APPLY_CONSENSUS {

    label 'process_medium'

    tag "${meta.id}:${meta.genotype}"

    // Mulled container: bcftools 1.21 + minimap2 2.28 + samtools 1.21.
    // This is the same mulled image as MINIMAP2_ROUND1 / MINIMAP2_CONSENSUS_MAP,
    // which also carries bcftools in the same image build.
    // If bcftools is not present in that image, the conda fallback installs it.
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a969b254c884786eea2ca247ced:3161f532a5ea6f1ade5f7b9af6e853a844a2d2a3-0'
    conda 'bioconda::bcftools=1.21 bioconda::minimap2=2.28 bioconda::samtools=1.21'

    publishDir (
        path: { "${params.outdir}/${meta.id}/consensus/${meta.genotype}/" },
        mode: 'copy',
        pattern: { "${meta.id}_${meta.genotype}_consensus.*" }
    )

    input:
    tuple val(meta), path(ref_fasta), path(vcf), path(vcf_csi), path(mask_bed)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_consensus.fasta"),
          path("${meta.id}_${meta.genotype}_consensus.fasta.fai"),
          path("${meta.id}_${meta.genotype}_consensus.mmi"),
          emit: consensus
    tuple val(meta), path("${meta.id}_${meta.genotype}.LOW_COVERAGE_CONSENSUS"),
          optional: true,
          emit: low_cov_sentinel
    path "versions.yml", emit: versions

    script:
    def consensus_name = "${meta.id}_${meta.genotype}_consensus"
    """
    # ----------------------------------------------------------------
    # Step 1: Apply variants from the filtered VCF to the dominant
    # reference, masking low-coverage positions as N.
    #
    # bcftools consensus:
    #   -f <ref>    reference FASTA to modify
    #   -m <bed>    BED of positions to mask with N
    #   <vcf>       filtered majority-allele VCF from BCFTOOLS_CONSENSUS_CALL
    # ----------------------------------------------------------------
    bcftools consensus \\
        -f ${ref_fasta} \\
        -m ${mask_bed} \\
        -o consensus_raw.fasta \\
        ${vcf}

    # ----------------------------------------------------------------
    # Step 2: Rename the FASTA header to the stable identifier.
    #
    # Use Python instead of `sed -i` to avoid BSD vs GNU sed portability
    # differences (macOS sed -i requires an empty backup suffix; GNU sed
    # does not — and containers are Linux/GNU but local dev is macOS).
    # ----------------------------------------------------------------
    python3 - <<'PYEOF'
    import re, sys
    stable_id = "${consensus_name}"
    with open("consensus_raw.fasta") as fh_in, \\
         open("${consensus_name}.fasta", "w") as fh_out:
        for line in fh_in:
            if line.startswith(">"):
                fh_out.write(">" + stable_id + "\\n")
            else:
                fh_out.write(line)
    PYEOF

    # ----------------------------------------------------------------
    # Step 3: Index the consensus FASTA with samtools faidx.
    # ----------------------------------------------------------------
    samtools faidx ${consensus_name}.fasta

    # ----------------------------------------------------------------
    # Step 4: Build a minimap2 binary index for Round 2 mapping.
    # Must use -x map-ont (same preset as all other minimap2 steps).
    # ----------------------------------------------------------------
    minimap2 \\
        -x map-ont \\
        -d ${consensus_name}.mmi \\
        ${consensus_name}.fasta

    # ----------------------------------------------------------------
    # Step 5: Pass/fail quality checks.
    #
    # (a) Length check — HCV genome is ~9,646 bp; complete consensus
    #     should be 4,000–11,000 bp (allowing partial and slightly
    #     over-long assemblies).  Out-of-range → WARNING to stderr.
    #
    # (b) N-fraction check — if >= 30% of positions are masked (N or n),
    #     the consensus is unreliable and we emit a sentinel file.
    #     The process itself exits 0 so the sample continues; the sentinel
    #     is routed to the per-sample report.
    # ----------------------------------------------------------------
    length=\$(grep -v '^>' ${consensus_name}.fasta | tr -d '\\n' | wc -c | tr -d ' ')
    n_count=\$(grep -v '^>' ${consensus_name}.fasta | tr -d '\\n' | tr -cd 'Nn' | wc -c | tr -d ' ')

    if [ "\${length}" -eq 0 ]; then
        echo "ERROR: consensus FASTA for ${meta.id} (${meta.genotype}) has zero length." >&2
        exit 1
    fi

    n_frac=\$(python3 -c "print(\${n_count} / \${length})")

    if [ "\${length}" -lt 4000 ] || [ "\${length}" -gt 11000 ]; then
        echo "WARNING: consensus length \${length} bp is outside the expected range [4000, 11000] for sample ${meta.id} genotype ${meta.genotype}." >&2
    fi

    # Emit LOW_COVERAGE_CONSENSUS sentinel when N-fraction >= 30%
    low_cov=\$(python3 -c "import sys; sys.exit(0 if \${n_count} / \${length} >= 0.30 else 1)") && touch ${meta.id}_${meta.genotype}.LOW_COVERAGE_CONSENSUS || true

    echo "Consensus QC — length: \${length} bp, N count: \${n_count}, N fraction: \${n_frac}" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def consensus_name_stub = "${meta.id}_${meta.genotype}_consensus"
    """
    printf '>${consensus_name_stub}\\nACGTACGTACGT\\n' > ${consensus_name_stub}.fasta
    touch ${consensus_name_stub}.fasta.fai
    touch ${consensus_name_stub}.mmi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "1.21"
        minimap2: "2.28"
        samtools: "1.21"
    END_VERSIONS
    """
}
