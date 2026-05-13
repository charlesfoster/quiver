/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VARIANT_FILTER — Post-call AF and depth filtering of LoFreq VCF.

    Purpose:
        Implements Step 5.16 of the data flow specification.  Applies a minimum
        allele-frequency and minimum depth filter to the raw LoFreq VCF, producing
        a filtered VCF and a human-readable TSV for reporting.

    bcftools inclusion filter (not exclusion):
        The filter uses `bcftools view -i` (include sites matching the expression),
        NOT `-e` (exclude).  Using `-i` is idiomatic bcftools for AF-threshold
        filtering and avoids the double-negative logic of `-e`.

    LoFreq INFO fields used:
        AF      Allele frequency of the ALT allele.  Range: (0, 1].
                All LoFreq calls have AF > 0; the filter floor is params.min_report_af
                (default 0.01 = 1%).
        DP      Total depth at the position after quality filters.  The filter
                floor is params.min_variant_depth (default 20).
        SB      Strand-bias score (Phred-scaled).  Reported in the TSV but not
                used as a filter here; analysts can apply an SB threshold in
                downstream tools.
        DP4     Read counts in four categories: ref-fwd, ref-rev, alt-fwd, alt-rev.
                Not filtered here; available for manual inspection in the raw VCF.
        INDEL   Flag field (present for indel calls, absent for SNVs).

    Expression operator note:
        bcftools filter expression uses `&&` for AND (not `&`).  Both forms work
        in bcftools 1.21 for the `-i` flag, but `&&` is the documented form.

    TSV report columns:
        CHROM   reference sequence name
        POS     1-based genomic position
        REF     reference allele
        ALT     alternate allele
        AF      allele frequency (INFO/AF)
        DP      depth (INFO/DP)
        SB      strand bias (INFO/SB)

    Inputs:
        meta  — val map with `id` and `genotype` fields
        vcf   — bgzip-compressed, tabix-indexed LoFreq VCF from LOFREQ_CALL
        tbi   — tabix index (.tbi)

    Outputs:
        vcf      — [meta, "*_filtered.vcf.gz", "*_filtered.vcf.gz.tbi"]
        tsv      — [meta, "*_variants.tsv"]
        versions — versions.yml

    Container: quay.io/biocontainers/bcftools:1.21--h8b25389_0
    Label: process_low (1 CPU, 2 GB, 5 min — Step 5.16 resource spec).

    Output published to:
        ${params.outdir}/${meta.id}/variants/${meta.genotype}/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process VARIANT_FILTER {

    label 'process_low'

    tag "${meta.id}:${meta.genotype}"

    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    conda 'bioconda::bcftools=1.21'

    publishDir (
        path: { "${params.outdir}/${meta.id}/variants/${meta.genotype}/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_filtered.vcf.gz"),
          path("${meta.id}_${meta.genotype}_filtered.vcf.gz.tbi"),
          emit: vcf
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_variants.tsv"),
          emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Filter the raw LoFreq VCF.
    #
    # -i  inclusion expression (keep sites matching the filter).
    #     NOT -e (exclusion) — use -i for AF threshold filtering.
    #
    # LoFreq INFO fields:
    #   AF  allele frequency (float, range 0–1).
    #   DP  total depth at the position.
    #
    # Retains only variants with:
    #   AF >= params.min_report_af  (default 0.01 = 1%)
    #   DP >= params.min_variant_depth (default 20)
    #
    # -Oz  output gzip-compressed VCF.
    # ----------------------------------------------------------------
    bcftools view \\
        -i "AF>=${params.min_report_af} && DP>=${params.min_variant_depth}" \\
        -Oz \\
        -o ${meta.id}_${meta.genotype}_filtered.vcf.gz \\
        ${vcf}

    tabix -p vcf ${meta.id}_${meta.genotype}_filtered.vcf.gz

    # ----------------------------------------------------------------
    # TSV report — one row per passing variant.
    #
    # Columns: CHROM, POS, REF, ALT, AF, DP, SB
    # SB (strand-bias Phred score) is included for QC; analysts can
    # apply their own SB threshold in downstream tools.
    # ----------------------------------------------------------------
    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\t%INFO/DP\\t%INFO/SB\\n' \\
        ${meta.id}_${meta.genotype}_filtered.vcf.gz \\
        > ${meta.id}_${meta.genotype}_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.genotype}_filtered.vcf.gz
    touch ${meta.id}_${meta.genotype}_filtered.vcf.gz.tbi
    touch ${meta.id}_${meta.genotype}_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "1.21"
    END_VERSIONS
    """
}
