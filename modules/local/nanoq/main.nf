/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANOQ — Fast read statistics, nanoq-compatible JSON output.

    The nanoq binary has no osx-arm64 build for any release. This module
    therefore calls bin/nanoq_stats.py, a pure-Python drop-in that produces
    identical JSON (reads, bases, n50, longest, shortest, mean_length,
    median_length, mean_quality, median_quality). MultiQC's nanoq parser
    and the sample report both consume this schema unchanged.

    Output filename: ${prefix}.nanoq.json where prefix defaults to ${meta.id}.
    Aliased invocations (NANOQ_FILT, NANOQ_POSTHOST) set task.ext.prefix via
    conf/modules.config to avoid filename collisions when files land in MultiQC.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process NANOQ {

    label 'process_low'

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://python:3.11' :
        'python:3.11' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/qc/raw/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.nanoq.json"), emit: json
    path "versions.yml",                           emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    nanoq_stats.py \\
        -i ${reads} \\
        -r ${prefix}.nanoq.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq_stats: \$(python3 --version 2>&1 | sed 's/Python /python3-/')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"reads":0,"bases":0,"n50":0,"longest":0,"shortest":0,"mean_length":0.0,"median_length":0.0,"mean_quality":0.0,"median_quality":0.0}' \\
        > ${prefix}.nanoq.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq_stats: python3-3.11
    END_VERSIONS
    """
}
