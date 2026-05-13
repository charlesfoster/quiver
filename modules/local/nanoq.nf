/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANOQ — Fast read statistics using nanoq (Rust implementation).

    nanoq is run with --json to emit machine-parseable stats and --report to
    produce a summary for MultiQC.

    The output JSON is named ${meta.id}.nanoq.json.

    Output is published to:
        ${params.outdir}/${meta.id}/qc/raw/

    This process is designed to be called on both raw reads (RAW_QC subworkflow)
    and post-host-depleted reads (step 5.5). The publishDir is set by the calling
    context; this module publishes unconditionally to the raw QC location.
    For post-host QC, override via a separate process invocation or subworkflow.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process NANOQ {

    label 'process_low'

    tag "${meta.id}"

    container 'quay.io/biocontainers/nanoq:0.10.0--h031d066_2'
    conda 'bioconda::nanoq=0.10.0'

    publishDir (
        path: { "${params.outdir}/${meta.id}/qc/raw/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.nanoq.json"), emit: json
    path "versions.yml",                            emit: versions

    script:
    """
    nanoq \\
        -i ${reads} \\
        --json \\
        -o ${meta.id}.nanoq.json \\
        --report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version 2>&1 | sed 's/nanoq //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"reads": 0, "bases": 0}' > ${meta.id}.nanoq.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: "0.10.0"
    END_VERSIONS
    """
}
