/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FILTER_VCF_FOR_DEVIDER — DEVIDER-specific AF threshold filter on the LoFreq VCF.

    Purpose:
        Applies a higher AF threshold filter to the VARIANT_FILTER output VCF
        before it is passed to DEVIDER.  The analysis VCF (VARIANT_FILTER output)
        uses params.min_report_af (default 1%) to report all low-AF variants.
        DEVIDER's de Bruijn graph-based phasing, however, is degraded by a large
        number of low-AF variants (many positions → reads look unique at every
        site → no coherent haplotype grouping → graph collapses to 1 haplotype).

        Setting params.devider_min_af (default 5%) removes the lowest-confidence
        variants from the phasing graph while still capturing variants at
        frequencies well above the DEVIDER joint detection floor (min-abund AND
        min-cov must both be satisfied).

    Why separate from VARIANT_FILTER:
        VARIANT_FILTER produces the analysis VCF used for reporting and consensus
        calling.  Changing its AF threshold would suppress reported variants.
        Keeping DEVIDER's input filter separate preserves the analysis VCF at 1%
        while tuning DEVIDER's phasing SNP set independently.

    Input:
        meta   — val map with `id` and `genotype` fields
        vcf    — VARIANT_FILTER output VCF (.vcf.gz, tabix-indexed; PASS-only, DP>=20)
        tbi    — tabix index

    Output:
        vcf    — [meta, "*_devider_input.vcf.gz", "*_devider_input.vcf.gz.tbi"]
        versions — versions.yml

    Container: quay.io/biocontainers/bcftools:1.21--h8b25389_0
    Label: process_low

    Output NOT published — this VCF is an intermediate consumed only by DEVIDER.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FILTER_VCF_FOR_DEVIDER {

    label 'process_low'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.21--h8b25389_0' }"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_devider_input.vcf.gz"),
          path("${meta.id}_${meta.genotype}_devider_input.vcf.gz.tbi"),
          emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    bcftools view \\
        -i "AF>=${params.devider_min_af}" \\
        -Oz \\
        -o ${meta.id}_${meta.genotype}_devider_input.vcf.gz \\
        ${vcf}

    tabix -p vcf ${meta.id}_${meta.genotype}_devider_input.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.genotype}_devider_input.vcf.gz
    touch ${meta.id}_${meta.genotype}_devider_input.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "1.21"
    END_VERSIONS
    """
}
