/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHOPPER — Length/quality read filter for Oxford Nanopore reads.

    Pipeline: zcat | chopper | pigz

    Parameters consumed from nextflow.config:
        params.min_qual    (default 8)
        params.min_length  (default 200)
        params.max_length  (default 10000)

    Chopper stderr is captured to chopper.log.  pigz is bundled in the chopper
    biocontainer; if it is not available in the container, replace with gzip.

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

    container 'quay.io/biocontainers/chopper:0.9.2--hdcf5f25_0'
    conda 'bioconda::chopper=0.9.2'

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
        | pigz -p ${task.cpus} > ${meta.id}.filtered.fastq.gz

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
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //')
    END_VERSIONS
    """

    stub:
    """
    echo -n "" | pigz > ${meta.id}.filtered.fastq.gz
    touch chopper.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "0.9.2"
        pigz: "2.8"
    END_VERSIONS
    """
}
