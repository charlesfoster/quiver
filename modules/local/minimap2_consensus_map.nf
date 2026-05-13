/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MINIMAP2_CONSENSUS_MAP — Map per-genotype reads to the dominant reference
    sequence in preparation for consensus building.

    Purpose:
        Step 5.10 sub-step 2.  Maps the per-genotype FASTQ against the single
        dominant panel reference extracted by EXTRACT_REF.  The resulting BAM
        feeds BCFTOOLS_CONSENSUS_CALL and BUILD_CONSENSUS_FASTA.

        Unlike Round 1 (competitive across 238 refs), this is a single-reference
        mapping; `-N 5 --secondary=no` are therefore dropped.  Flags retained:
            -Y      soft-clip supplementary (required by LoFreq in later steps)
            --MD    MD mismatch tag for IGV / bcftools
            --eqx   =/X CIGAR ops for bcftools mpileup accuracy

        The read-group tag (`-R`) uses the same SM field as Round 1 so that the
        BAM header is consistent across the pipeline.

    Inputs:
        meta        — val map with `id` and `genotype` fields
        ref_fasta   — single-record dominant reference FASTA from EXTRACT_REF
        reads       — per-genotype FASTQ (host-depleted, possibly partitioned)

    Outputs:
        bam         — [meta, "*_consensus_map.bam", "*_consensus_map.bam.bai"]
        versions    — versions.yml

    Container: mulled minimap2 2.28 + samtools 1.21 (same image as MINIMAP2_ROUND1).
    Label: process_high (16 CPU, 32 GB — Step 5.10 resource spec).

    Output published to:
        ${params.outdir}/${meta.id}/consensus/${meta.genotype}/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MINIMAP2_CONSENSUS_MAP {

    label 'process_high'

    tag "${meta.id}:${meta.genotype}"

    // Mulled container providing minimap2 2.28 + samtools 1.21.
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a969b254c884786eea2ca247ced:3161f532a5ea6f1ade5f7b9af6e853a844a2d2a3-0'
    conda 'bioconda::minimap2=2.28 bioconda::samtools=1.21'

    publishDir (
        path: { "${params.outdir}/${meta.id}/consensus/${meta.genotype}/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(ref_fasta), path(reads)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_consensus_map.bam"),
          path("${meta.id}_${meta.genotype}_consensus_map.bam.bai"),
          emit: bam
    path "versions.yml", emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Map per-genotype reads to the single dominant reference.
    # Pipe directly to samtools sort — no intermediate SAM on disk.
    #
    # Flag rationale (docs/architecture_reasoning.md §6):
    #   -ax map-ont   ONT preset (k=15, w=10)
    #   -Y            soft-clip supplementary (required by LoFreq later)
    #   --MD          mismatch tag for IGV / bcftools mpileup
    #   --eqx         =/X CIGAR for bcftools accuracy
    #   -R            read group — consistent with Round 1 header
    # ----------------------------------------------------------------
    minimap2 \\
        -ax map-ont \\
        -t ${task.cpus} \\
        -Y \\
        --MD \\
        --eqx \\
        -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ONT" \\
        ${ref_fasta} \\
        ${reads} \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -O bam \\
        -o ${meta.id}_${meta.genotype}_consensus_map.bam \\
        -

    # ----------------------------------------------------------------
    # Index the sorted BAM.
    # ----------------------------------------------------------------
    samtools index -@ ${task.cpus} ${meta.id}_${meta.genotype}_consensus_map.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.genotype}_consensus_map.bam
    touch ${meta.id}_${meta.genotype}_consensus_map.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "2.28"
        samtools: "1.21"
    END_VERSIONS
    """
}
