/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAKE_ROUND2_MASK — Coverage mask generation from the Round 2 BAM.

    Runs mosdepth in quantize mode to produce a BED of regions below
    params.min_consensus_cov.  These positions will be written as N
    in the Round 2 consensus FASTAs (design decision D8).

    Part of the BUILD_ROUND2_CONSENSUS workflow (Step 5.16b).

    Container: mosdepth 0.3.14
    Label: process_low
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MAKE_ROUND2_MASK {

    label 'process_low'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.14--h05c3d44_0' :
        'quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_0' }"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("r2_mask.bed"), emit: mask_bed
    path "versions.yml",                  emit: versions

    script:
    def threshold = params.min_consensus_cov
    """
    mosdepth \\
        --quantize 0:${threshold}: \\
        --no-per-base \\
        ${meta.id}_${meta.genotype}_r2cov \\
        ${bam}

    LOW_BIN="0:${threshold}"
    zcat ${meta.id}_${meta.genotype}_r2cov.quantized.bed.gz \\
        | awk -v low_bin="\${LOW_BIN}" '\$4 == low_bin { print \$1"\\t"\$2"\\t"\$3 }' \\
        > r2_mask.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """

    stub:
    """
    touch r2_mask.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: "0.3.14"
    END_VERSIONS
    """
}
