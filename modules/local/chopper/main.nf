/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHOPPER — Length/quality read filter for Oxford Nanopore reads.

    Pipeline: zcat | chopper | gzip

    Parameters consumed from nextflow.config:
        params.min_qual    (default 8)
        params.min_length  (default 200)
        params.max_length  (default 10000)

    Chopper stderr is captured to chopper.log.

    Sentinel behaviour:
        If the output FASTQ contains zero reads (empty file or gzip of empty FASTQ),
        a sentinel file ${meta.id}.ALL_READS_FILTERED is emitted alongside the
        (empty) filtered FASTQ.  The process does NOT fail — downstream processes
        must gate on the absence of this sentinel.

    Output published to:
        ${params.outdir}/${meta.id}/reads/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CHOPPER {

    label 'process_medium'

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chopper:0.9.2--hcdda2d0_0' :
        'quay.io/biocontainers/chopper:0.9.2--hcdda2d0_0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/reads/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filtered.fastq.gz"), emit: reads
    tuple val(meta), path("chopper.log"),                  emit: log
    tuple val(meta), path("${meta.id}.ALL_READS_FILTERED"),
          optional: true,                                  emit: sentinel
    path "versions.yml",                                   emit: versions

    script:
    """
    # Run the filter pipeline; chopper stderr goes to chopper.log
    zcat ${reads} \\
        | chopper \\
            -q ${params.min_qual} \\
            --minlength ${params.min_length} \\
            --maxlength ${params.max_length} \\
            --threads ${task.cpus} \\
            2> chopper.log \\
        | gzip -c > ${meta.id}.filtered.fastq.gz

    # Detect whether any reads survived filtering.
    # zcat on an empty gzip returns nothing; awk 'NR==1' checks for at least
    # one line (the first FASTQ header).
    n_reads=\$(zcat ${meta.id}.filtered.fastq.gz | awk 'NR%4==1' | wc -l | tr -d ' ')
    if [ "\${n_reads}" -eq 0 ]; then
        echo "WARNING: All reads were filtered out for sample ${meta.id}." >&2
        touch ${meta.id}.ALL_READS_FILTERED
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | sed 's/chopper //')
        gzip: \$(gzip --version 2>&1 | head -1 | sed 's/gzip //')
    END_VERSIONS
    """

    stub:
    """
    echo -n "" | gzip -c > ${meta.id}.filtered.fastq.gz
    touch chopper.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "0.9.2"
        gzip: "1.12"
    END_VERSIONS
    """
}
