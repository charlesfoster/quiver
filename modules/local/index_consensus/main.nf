/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INDEX_CONSENSUS — Build samtools .fai and minimap2 .mmi indices for the
    sample-specific consensus FASTA.

    Purpose:
        Step 5.10 sub-step 5.  Indexes the consensus FASTA produced by
        APPLY_CONSENSUS so that MINIMAP2_ROUND2 (and downstream modules) can
        use it directly without re-computing the index.

    Steps:
        1. samtools faidx: build .fai index.
        2. minimap2 -x map-ont -d: build .mmi index for Round 2 mapping.

    The downstream workflow expects a tuple [meta, fasta, fai, mmi] from
    BUILD_CONSENSUS.out.consensus.  This process produces the fai and mmi
    components and emits the full tuple.

    Inputs:
        meta            — val map with `id` and `genotype` fields
        consensus_fasta — consensus FASTA from APPLY_CONSENSUS

    Outputs:
        consensus — [meta, fasta, fai, mmi]  single tuple for Round 2 input
        versions  — versions.yml

    Container: mulled minimap2 2.28 + samtools 1.21
    Tag: "${meta.id}:${meta.genotype}"
    Label: process_low (fast indexing step)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process INDEX_CONSENSUS {

    label 'process_low'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(consensus_fasta)

    output:
    tuple val(meta),
          path(consensus_fasta),
          path("${consensus_fasta}.fai"),
          path("${consensus_fasta.baseName}.mmi"),
          emit: consensus
    path "versions.yml", emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Step 1: Index the consensus FASTA with samtools faidx.
    # Produces ${consensus_fasta}.fai alongside the FASTA.
    # ----------------------------------------------------------------
    samtools faidx ${consensus_fasta}

    # ----------------------------------------------------------------
    # Step 2: Build a minimap2 binary index for Round 2 mapping.
    # Must use -x map-ont (same preset as all other minimap2 steps).
    # ----------------------------------------------------------------
    minimap2 \\
        -x map-ont \\
        -t ${task.cpus} \\
        -d ${consensus_fasta.baseName}.mmi \\
        ${consensus_fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch ${consensus_fasta}.fai
    touch ${consensus_fasta.baseName}.mmi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "1.21"
        minimap2: "2.28"
    END_VERSIONS
    """
}
