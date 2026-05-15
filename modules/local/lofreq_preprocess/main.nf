/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOFREQ_PREPROCESS — Prepare a BAM for LoFreq variant calling via indel-quality
    scoring.

    Purpose:
        Implements Step 5.12 of the data flow specification.  LoFreq was designed for
        Illumina data and does not natively handle ONT-specific error profiles.  This
        preprocessing step adds indel quality tags required for indel calling:

        `lofreq indelqual --dindel`
        Calibrates per-read indel quality scores using the DINDEL algorithm.
        Without this step, homopolymer errors in ONT reads generate explosive
        indel false-positives because LoFreq interprets the raw base qualities
        as reliable indel evidence.  It produces BAM tags BI (indel base quality)
        and BD (deletion base quality) that LoFreq's statistical model uses during
        variant calling.

    Alignment-quality note:
        The previous pipeline ran `lofreq alnqual -b` after indelqual.  On Docker
        Desktop for Apple Silicon this amd64 LoFreq command was repeatedly killed
        under emulation, so the integration has been removed.  The caller consumes
        the indelqual BAM directly.

    Prerequisites for LoFreq preprocessing to work (all guaranteed upstream):
        - BAM is sorted and indexed.
        - minimap2 was run with `-Y` (soft-clip supplementary).
        - minimap2 was run with `--MD --eqx` (MD tag and =/X CIGAR).
        - BAM has a read group (`-R` flag in minimap2).
        See docs/architecture_reasoning.md §7 and CLAUDE.md D3, D9.

    Inputs:
        meta        — val map with `id` and `genotype` fields
        bam         — sorted, indexed Round 2 BAM (subsampled for LoFreq)
        bai         — BAM index
        ref_fasta   — per-genotype consensus FASTA (same reference used in mapping)

    Outputs:
        bam         — [meta, "*_preprocessed.bam", "*_preprocessed.bam.bai"]
                      The indelqual BAM, ready for lofreq call.
        versions    — versions.yml

    Container: quay.io/biocontainers/lofreq:2.1.5--py310h4966b78_15
    Label: process_medium (4 CPU, 8 GB, 30 min — Step 5.12 resource spec).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process LOFREQ_PREPROCESS {

    label 'process_medium'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py310h4966b78_15' :
        'quay.io/biocontainers/lofreq:2.1.5--py310h4966b78_15' }"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam), path(bai), path(ref_fasta)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_preprocessed.bam"),
          path("${meta.id}_${meta.genotype}_preprocessed.bam.bai"),
          emit: bam
    path "versions.yml", emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Step 1: Indel quality scoring with DINDEL algorithm.
    #
    # lofreq indelqual --dindel adds per-read BI (indel base quality)
    # and BD (deletion base quality) tags.  Without these tags,
    # lofreq call will refuse to run with --call-indels.
    #
    # Prerequisite: BAM must have been sorted, indexed, and mapped with
    # -Y --MD --eqx (guaranteed by MINIMAP2_ROUND2).
    #
    # Output: ${meta.id}_${meta.genotype}_iq.bam (unsorted; already sorted)
    # ----------------------------------------------------------------
    lofreq indelqual \\
        --dindel \\
        -f ${ref_fasta} \\
        -o ${meta.id}_${meta.genotype}_iq.bam \\
        ${bam}

    samtools index ${meta.id}_${meta.genotype}_iq.bam

    # ----------------------------------------------------------------
    # Rename to the canonical output name consumed downstream.
    # ----------------------------------------------------------------
    mv ${meta.id}_${meta.genotype}_iq.bam     ${meta.id}_${meta.genotype}_preprocessed.bam
    mv ${meta.id}_${meta.genotype}_iq.bam.bai ${meta.id}_${meta.genotype}_preprocessed.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(lofreq version 2>&1 | head -1 | sed 's/lofreq version //')
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.genotype}_preprocessed.bam
    touch ${meta.id}_${meta.genotype}_preprocessed.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: "2.1.5"
        samtools: "1.21"
    END_VERSIONS
    """
}
