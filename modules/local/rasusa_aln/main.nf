/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RASUSA_ALN — Coverage-accurate depth-capping of an indexed BAM.

    Uses `rasusa aln` which calculates coverage from actual per-position alignment
    depth rather than estimating from read count × genome size.  This is more
    accurate than `rasusa reads` and avoids a redundant minimap2 remapping step:
    the full-depth round 2 BAM is subsampled directly.

    rasusa aln flag notes:
        --coverage    Target depth; rasusa passes all reads through if already below.
        --seed        Random seed for reproducibility.
        -o            Output BAM path.
        <FILE>        Indexed BAM is a positional argument.

    The input BAM must be indexed (.bai alongside the BAM path).

    Inputs:
        meta      — val map with `id` and `genotype` fields
        bam       — full-depth indexed BAM
        bai       — BAM index
        coverage  — integer coverage cap
        seed      — integer random seed

    Outputs:
        bam       — [meta, "*_subsampled.bam", "*_subsampled.bam.bai"]
        versions  — versions.yml

    Container: quay.io/biocontainers/rasusa:2.1.0--hc1c3326_1
    Label: process_medium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RASUSA_ALN {

    label 'process_medium'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rasusa:2.1.0--hc1c3326_1' :
        'quay.io/biocontainers/rasusa:2.1.0--hc1c3326_1' }"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam), path(bai), val(coverage), val(seed)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_subsampled.bam"),
          path("${meta.id}_${meta.genotype}_subsampled.bam.bai"),
          emit: bam
    path "versions.yml", emit: versions

    script:
    """
    rasusa aln \\
        --coverage ${coverage} \\
        --seed ${seed} \\
        -o ${meta.id}_${meta.genotype}_subsampled.bam \\
        ${bam}

    samtools index ${meta.id}_${meta.genotype}_subsampled.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: \$(rasusa --version 2>&1 | sed 's/rasusa //')
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.genotype}_subsampled.bam
    touch ${meta.id}_${meta.genotype}_subsampled.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: "2.1.0"
        samtools: "1.21"
    END_VERSIONS
    """
}
