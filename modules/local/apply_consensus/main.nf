/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    APPLY_CONSENSUS — Apply variants and low-coverage masking to produce the
    sample-specific consensus FASTA.

    Purpose:
        Step 5.10 sub-step 4.  Applies the majority-allele VCF from
        BCFTOOLS_CONSENSUS_CALL to the dominant reference FASTA, masking
        low-coverage positions (from MAKE_MASK_BED) as N.  Renames the FASTA
        header to a stable identifier.  Emits a LOW_COVERAGE_CONSENSUS sentinel
        if >= 30% of positions are masked.

        The samtools faidx and minimap2 -d indexing steps have been moved to a
        separate INDEX_CONSENSUS process so that this process only requires the
        bcftools container.

    Steps:
        1. bcftools consensus: apply filtered VCF variants and mask.
        2. awk one-liner: rename the FASTA header to the stable identifier
           "${meta.id}_${meta.genotype}_consensus".
           (awk used instead of python3/sed for portability across containers.)
        3. Pass/fail check: warn if length outside [4000,11000] bp;
           emit LOW_COVERAGE_CONSENSUS sentinel if N-fraction >= 30%.

    Container: bcftools 1.21
    Label: process_medium

    Output published to:
        ${params.outdir}/${meta.id}/consensus/${meta.genotype}/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process APPLY_CONSENSUS {

    label 'process_medium'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.21--h8b25389_0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/consensus/${meta.genotype}/" },
        mode: 'copy',
        pattern: '*_consensus.fasta'
    )

    input:
    tuple val(meta), path(ref_fasta), path(vcf), path(vcf_csi), path(mask_bed)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_consensus.fasta"),
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
    # Use awk instead of python3 or sed -i to avoid portability issues.
    # NR==1 matches only the first line (the header); subsequent lines
    # (sequence) are passed through unchanged.
    # ${consensus_name} is a Groovy variable interpolated before the
    # script block runs — it is NOT a shell variable.
    # ----------------------------------------------------------------
    awk 'NR==1 { print ">${consensus_name}"; next } { print }' consensus_raw.fasta > ${consensus_name}.fasta

    # ----------------------------------------------------------------
    # Step 3: Pass/fail quality checks.
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

    n_frac=\$(awk "BEGIN { print \${n_count} / \${length} }")

    if [ "\${length}" -lt 4000 ] || [ "\${length}" -gt 11000 ]; then
        echo "WARNING: consensus length \${length} bp is outside the expected range [4000, 11000] for sample ${meta.id} genotype ${meta.genotype}." >&2
    fi

    # Emit LOW_COVERAGE_CONSENSUS sentinel when N-fraction >= 30%
    awk "BEGIN { exit (\${n_count} / \${length} >= 0.30 ? 0 : 1) }" && touch ${meta.id}_${meta.genotype}.LOW_COVERAGE_CONSENSUS || true

    echo "Consensus QC — length: \${length} bp, N count: \${n_count}, N fraction: \${n_frac}" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def consensus_name_stub = "${meta.id}_${meta.genotype}_consensus"
    """
    printf '>${consensus_name_stub}\\nACGTACGTACGT\\n' > ${consensus_name_stub}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "1.21"
    END_VERSIONS
    """
}
