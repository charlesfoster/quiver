/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INDEX_PANEL — Build a minimap2 index and samtools FAI for the HCV reference panel.

    This process runs ONCE per pipeline run.  Nextflow caches the result by the
    content hash of the input FASTA, so re-runs that share the same panel file
    never re-index.

    IMPORTANT: The process input is only `path fasta` (no sample-specific context).
    Keeping the cache key independent of samples is design requirement D14.

    Steps performed:
      1. Copy the panel FASTA to the work directory so that samtools faidx can
         write the .fai index alongside it (samtools requires the .fai to live
         next to the FASTA it describes).
      2. Build the minimap2 index with -x map-ont (k=15, preset for ONT long reads).
      3. Index the FASTA with samtools faidx.

    Outputs emitted as a single tuple [panel_mmi, panel_fasta, panel_fai]:
        panel.mmi        — minimap2 binary index
        panel.fasta      — copy of the reference panel FASTA (in work dir)
        panel.fasta.fai  — samtools FASTA index

    Container: mulled image providing both minimap2 2.28 and samtools 1.21.
    Label: process_medium (4 CPU, 8 GB, matches Step 5.6 resource spec).

    Output published to:
        ${params.outdir}/reference/panel/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process INDEX_PANEL {

    label 'process_medium'

    tag "panel_index"

    // Mulled container providing minimap2 2.28 + samtools 1.21 in a single image.
    // This is the same image used by HOST_DEPLETE_MINIMAP2 and MINIMAP2_ROUND1.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: "${params.outdir}/reference/panel/",
        mode: 'copy'
    )

    input:
    path fasta   // HCV reference panel FASTA; content hash is the cache key

    output:
    tuple path("panel.mmi"), path("panel.fasta"), path("panel.fasta.fai"), emit: index
    path "versions.yml",                                                    emit: versions

    script:
    """
    # Step 1 — copy FASTA into the work directory so samtools faidx can write
    # the .fai alongside it.  A symlink is insufficient because samtools faidx
    # writes ${fasta}.fai relative to the file's location, which would be
    # outside the work directory if we used a symlink to the staging path.
    cp ${fasta} panel.fasta

    # Step 2 — build the minimap2 binary index with the map-ont preset.
    # -x map-ont sets k=15, w=10, which is the correct kmer size for ONT reads.
    # -d writes the index to the named file.
    minimap2 \\
        -x map-ont \\
        -t ${task.cpus} \\
        -d panel.mmi \\
        panel.fasta

    # Step 3 — build the samtools FASTA index (.fai).
    samtools faidx panel.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch panel.mmi
    cp ${fasta} panel.fasta
    touch panel.fasta.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "2.28"
        samtools: "1.21"
    END_VERSIONS
    """
}
