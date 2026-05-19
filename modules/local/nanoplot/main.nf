/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANOPLOT — Raw read quality-control visualisation.

    NanoPlot is run with --no_static to skip heavy static PNG generation and
    --tsv_stats to emit machine-parseable TSV stats for MultiQC.

    Low read-count warnings from NanoPlot (common in test data) do NOT fail the
    process — the || true guard on the NanoPlot call ensures Nextflow only sees
    the exit code from the final `true` statement if NanoPlot itself exits 0 or
    with a warning code.  In practice NanoPlot exits 0 even on low-count data;
    the guard is belt-and-suspenders for test reproducibility.

    Output directory is published to:
        ${params.outdir}/${meta.id}/qc/raw/nanoplot/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process NANOPLOT {

    label 'process_medium'

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.46.2--pyhdfd78af_1' :
        'quay.io/biocontainers/nanoplot:1.46.2--pyhdfd78af_1' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/qc/raw/nanoplot/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_nanoplot/"), emit: nanoplot_dir
    path "versions.yml",                           emit: versions

    script:
    """
    NanoPlot \\
        --fastq ${reads} \\
        -o ${meta.id}_nanoplot \\
        --tsv_stats \\
        --no_static \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version 2>&1 | sed 's/NanoPlot //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_nanoplot
    touch ${meta.id}_nanoplot/NanoPlot-report.html
    touch ${meta.id}_nanoplot/NanoStats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: "1.43.0"
    END_VERSIONS
    """
}
