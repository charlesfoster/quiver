/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INDEX_HOST — Build a minimap2 .mmi index from the host reference FASTA.

    Called once per run (Nextflow caches by content hash of the FASTA).

    If params.host_reference is already a .mmi file, this process is skipped
    and the file is passed through directly by the calling workflow.

    If params.host_reference is null, the get_host_reference.sh script is
    invoked to download GRCh38 no-alt from NCBI before indexing.

    The index is built with -x map-ont to match the host-depletion alignment
    parameters in HOST_DEPLETE_MINIMAP2.

    Output published to:
        ${params.outdir}/reference/host/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process INDEX_HOST {

    label 'process_high'

    tag "host_index"

    container 'quay.io/biocontainers/minimap2:2.28--he4a0461_0'
    conda 'bioconda::minimap2=2.28'

    publishDir (
        path: "${params.outdir}/reference/host/",
        mode: 'copy',
        saveAs: { filename -> filename.endsWith('.mmi') ? filename : null }
    )

    input:
    path fasta   // host reference FASTA (possibly .gz); may be a pre-built .mmi

    output:
    path "host.mmi",   emit: index
    path "versions.yml", emit: versions

    script:
    """
    minimap2 \\
        -x map-ont \\
        -t ${task.cpus} \\
        -d host.mmi \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch host.mmi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "2.28"
    END_VERSIONS
    """
}
