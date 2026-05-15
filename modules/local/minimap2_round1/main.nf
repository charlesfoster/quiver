/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MINIMAP2_ROUND1 — Round 1 competitive mapping of host-depleted reads against
    the full HCV reference panel.

    Purpose:
        Assign each read to its best-matching HCV reference across all 238 panel
        sequences.  The resulting BAM feeds the genotype-classification step
        (bin/classify_genotype.py) which detects mixed-genotype infections and
        assigns reads to genotype branches.

    Key flag rationale (design decision D3 + docs/architecture_reasoning.md §6):
        -ax map-ont     ONT long-read preset (k=15, w=10).
        --secondary=no  Emit only the single primary alignment per read.
        -N 5            Retain up to 5 candidate alignments internally before
                        selecting the best primary; together with --secondary=no
                        this achieves competitive mapping: the scorer considers
                        multiple references but the BAM carries only the winner.
        -Y              Soft-clip supplementary alignments — required by LoFreq
                        (hard clipping corrupts its per-base quality recalibration).
        --MD            Emit the MD tag (mismatch string) for IGV and LoFreq.
        --eqx           Use =/X CIGAR ops instead of M — required by LoFreq's
                        alnqual step.
        -R              Read group tag — LoFreq errors without it.

    Pipeline:
        minimap2 → samtools sort (no intermediate SAM) → samtools index
        samtools flagstat → sentinel gate

    Sentinel logic:
        Count primary-mapped reads from flagstat output.
        If count < params.min_round1_mapped (default 100), write an empty
        ${meta.id}.NO_HCV_DETECTED file alongside the BAM.
        The process still exits 0 — the calling workflow routes flagged samples
        to EMIT_FAILURE_REPORT rather than downstream analysis (Prompt 24).

    Outputs:
        bam        — [meta, round1.bam, round1.bam.bai]
        flagstat   — [meta, round1.flagstat]
        sentinel   — [meta, NO_HCV_DETECTED]  (optional: true)
        versions   — versions.yml

    Container: mulled image providing both minimap2 2.28 and samtools 1.21.
    Label: process_high (16 CPU, 32 GB, 60 min — matches Step 5.7 resource spec).

    Output published to:
        ${params.outdir}/${meta.id}/genotyping/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MINIMAP2_ROUND1 {

    label 'process_high'

    tag "${meta.id}"

    // Mulled container providing minimap2 2.28 + samtools 1.21 in a single image.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/genotyping/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(reads)
    tuple path(panel_mmi), path(panel_fasta), path(panel_fai)

    output:
    tuple val(meta), path("${meta.id}.round1.bam"), path("${meta.id}.round1.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.round1.flagstat"),                                emit: flagstat
    tuple val(meta), path("${meta.id}.NO_HCV_DETECTED"),
          optional: true,                                                               emit: sentinel
    path "versions.yml",                                                               emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Round 1 competitive mapping.
    # Pipe directly to samtools sort — no intermediate SAM file.
    # ----------------------------------------------------------------
    minimap2 \\
        -ax map-ont \\
        -t ${task.cpus} \\
        --secondary=no \\
        -N 5 \\
        -Y \\
        --MD \\
        --eqx \\
        -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ONT" \\
        ${panel_mmi} \\
        ${reads} \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -O bam \\
        -o ${meta.id}.round1.bam \\
        -

    # ----------------------------------------------------------------
    # Index the sorted BAM
    # ----------------------------------------------------------------
    samtools index -@ ${task.cpus} ${meta.id}.round1.bam

    # ----------------------------------------------------------------
    # Flagstat — written to file for MultiQC and downstream gate check
    # ----------------------------------------------------------------
    samtools flagstat -@ ${task.cpus} ${meta.id}.round1.bam > ${meta.id}.round1.flagstat

    # ----------------------------------------------------------------
    # Gate check — emit NO_HCV_DETECTED sentinel if fewer than
    # params.min_round1_mapped primary reads mapped.
    # The process still exits 0 so the sample is routed, not failed.
    # ----------------------------------------------------------------
    mapped=\$(grep "primary mapped" ${meta.id}.round1.flagstat | awk '{print \$1}')
    if [ "\${mapped}" -lt "${params.min_round1_mapped}" ]; then
        echo "WARNING: Only \${mapped} primary reads mapped to HCV panel for sample ${meta.id} (threshold: ${params.min_round1_mapped}). Emitting NO_HCV_DETECTED sentinel." >&2
        touch ${meta.id}.NO_HCV_DETECTED
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.round1.bam
    touch ${meta.id}.round1.bam.bai
    printf "0 + 0 in total (QC-passed reads + QC-failed reads)\n0 + 0 primary\n0 + 0 secondary\n0 + 0 supplementary\n0 + 0 duplicates\n0 + 0 primary mapped (0.00%% : N/A)\n" \\
        > ${meta.id}.round1.flagstat
    touch ${meta.id}.NO_HCV_DETECTED

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "2.28"
        samtools: "1.21"
    END_VERSIONS
    """
}
