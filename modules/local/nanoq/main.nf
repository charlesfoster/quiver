/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANOQ — Fast read statistics using nanoq (Rust implementation).

    nanoq is run with --json to emit machine-parseable stats.

    Output filename: ${prefix}.nanoq.json where prefix defaults to ${meta.id}.
    Aliased invocations (NANOQ_FILT, NANOQ_POSTHOST) set task.ext.prefix via
    conf/base.config to avoid filename collisions when files are passed to MultiQC.

    This module publishes to qc/raw/ by default.  Aliased invocations
    (NANOQ_FILT, NANOQ_POSTHOST) override publishDir via conf/base.config
    withName directives so each lands in the correct subdirectory.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process NANOQ {

    label 'process_low'

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoq:0.10.0--h031d066_2' :
        'quay.io/biocontainers/nanoq:0.10.0--h031d066_2' }"
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
    nanoq \\
        -i ${reads} \\
        --json \\
        -r ${prefix}.nanoq.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version 2>&1 | sed 's/nanoq //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"reads": 0, "bases": 0}' > ${prefix}.nanoq.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: "0.10.0"
    END_VERSIONS
    """
}
