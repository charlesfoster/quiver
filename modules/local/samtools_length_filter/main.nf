/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SAMTOOLS_LENGTH_FILTER — Retain only reads >= min_length bp from an indexed BAM.

    Purpose:
        Pre-filters reads before DEVIDER haplotype reconstruction.  Short reads
        (<2,000 bp by default) span too few SNPs to contribute useful phasing
        information to DEVIDER's de Bruijn graph — they increase node count without
        improving graph connectivity, causing fragmentation at high coverage.

        Applied BEFORE rasusa_aln so the depth cap is drawn exclusively from reads
        that are long enough to phase.

    samtools flag notes:
        -e 'qlen>=N'   expression filter on query sequence length (unclipped)
        -b             output BAM
        -@             additional threads for compression

    Inputs:
        meta        — val map with `id` and `genotype` fields
        bam         — indexed BAM
        bai         — BAM index
        min_length  — integer minimum read length in bp

    Outputs:
        bam         — [meta, "*_lengthfiltered.bam", "*_lengthfiltered.bam.bai"]
        versions    — versions.yml

    Container: quay.io/biocontainers/samtools:1.21--h50ea8bc_0
    Label: process_low
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SAMTOOLS_LENGTH_FILTER {

    label 'process_low'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'quay.io/biocontainers/samtools:1.21--h50ea8bc_0' }"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam), path(bai), val(min_length)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_lengthfiltered.bam"),
          path("${meta.id}_${meta.genotype}_lengthfiltered.bam.bai"),
          emit: bam
    path "versions.yml", emit: versions

    script:
    """
    samtools view \\
        -e 'qlen>=${min_length}' \\
        -b \\
        -@ ${task.cpus} \\
        -o ${meta.id}_${meta.genotype}_lengthfiltered.bam \\
        ${bam}

    samtools index ${meta.id}_${meta.genotype}_lengthfiltered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.genotype}_lengthfiltered.bam
    touch ${meta.id}_${meta.genotype}_lengthfiltered.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "1.21"
    END_VERSIONS
    """
}
