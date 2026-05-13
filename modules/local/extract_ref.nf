/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXTRACT_REF — Extract a single reference sequence from the HCV panel FASTA.

    Purpose:
        The consensus-build subworkflow (Step 5.10) needs to map per-genotype
        reads to only the single dominant reference for that genotype — not the
        full 238-sequence panel.  This process uses `samtools faidx` to extract
        exactly that one record and writes it to a fresh FASTA file.

    The `dominant_ref_id` value originates from classify_genotype.py's
    `top_reference` field in the genotype_summary JSON.  The value is
    a FASTA record name that matches the panel header exactly (e.g.
    "1a_M62321.1"), and samtools faidx looks it up by that name.

    Inputs:
        meta             — val map with `id` (sample ID) and `genotype` fields
        panel_fasta      — the full HCV reference panel FASTA
        panel_fai        — samtools index for the panel FASTA (required by faidx)
        dominant_ref_id  — val string: the FASTA record name to extract

    Outputs:
        ref_fasta   — [meta, "${meta.id}_${meta.genotype}_dominant_ref.fasta"]
        versions    — versions.yml

    Container: samtools 1.21 (CLAUDE.md §4)
    Label: process_low (fast single-record extraction).

    Output published to:
        ${params.outdir}/${meta.id}/consensus/${meta.genotype}/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EXTRACT_REF {

    label 'process_low'

    tag "${meta.id}:${meta.genotype}"

    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    conda 'bioconda::samtools=1.21'

    publishDir (
        path: { "${params.outdir}/${meta.id}/consensus/${meta.genotype}/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(panel_fasta), path(panel_fai), val(dominant_ref_id)

    output:
    tuple val(meta), path("${meta.id}_${meta.genotype}_dominant_ref.fasta"), emit: ref_fasta
    path "versions.yml",                                                      emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Extract the single dominant reference record from the panel FASTA.
    # samtools faidx uses the .fai index; both panel_fasta and panel_fai
    # must be present in the work directory (declared as inputs above).
    # ----------------------------------------------------------------
    samtools faidx \\
        ${panel_fasta} \\
        "${dominant_ref_id}" \\
        > ${meta.id}_${meta.genotype}_dominant_ref.fasta

    # Guard: die loudly if the record was not found (empty output).
    if [ ! -s ${meta.id}_${meta.genotype}_dominant_ref.fasta ]; then
        echo "ERROR: samtools faidx produced an empty file for reference '${dominant_ref_id}'. " \\
             "Check that this ID exists verbatim in ${panel_fasta}." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    printf '>%s\\nACGTACGT\\n' "${dominant_ref_id}" \\
        > ${meta.id}_${meta.genotype}_dominant_ref.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "1.21"
    END_VERSIONS
    """
}
