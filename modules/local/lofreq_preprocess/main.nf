/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOFREQ_PREPROCESS — Prepare a BAM for LoFreq variant calling via indel-quality
    scoring and optional alignment-quality recalibration.

    Purpose:
        Implements Step 5.12 of the data flow specification.  LoFreq was designed for
        Illumina data and does not natively handle ONT-specific error profiles.  Two
        preprocessing steps bridge this gap:

        1. `lofreq indelqual --dindel`
           Calibrates per-read indel quality scores using the DINDEL algorithm.
           Without this step, homopolymer errors in ONT reads generate explosive
           indel false-positives because LoFreq interprets the raw base qualities
           as reliable indel evidence.
           Produces BAM tags BI (indel base quality) and BD (deletion base quality)
           that LoFreq's statistical model uses during variant calling.

        2. `lofreq alnqual -b`
           Recalibrates per-base alignment qualities by considering the local
           alignment context.  This dampens systematic ONT error sites (e.g.
           positions adjacent to homopolymers that generate correlated base errors
           at a fixed frequency).
           This step is optional in the sense that if it fails (non-zero exit),
           the pipeline falls back to the indelqual-only BAM.  In practice alnqual
           rarely fails when --eqx and --MD were used during mapping (D3).

    Fallback rationale:
        `lofreq alnqual` can fail on malformed CIGAR strings or unsupported MD
        tags.  Since MINIMAP2_ROUND2 uses `--MD --eqx`, this should not occur in
        normal operation.  However, the fallback (cp indelqual.bam → alnqual.bam)
        ensures the pipeline continues gracefully rather than losing the entire
        sample.  A WARNING is emitted to stderr so the operator is alerted.

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
                      The alnqual (or fallback indelqual) BAM, ready for lofreq call.
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
    # lofreq call-parallel will refuse to run with --call-indels.
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
    # Step 2: Alignment quality recalibration.
    #
    # lofreq alnqual -b recalibrates per-base alignment qualities using
    # the alignment context, dampening systematic ONT error sites.
    # The -b flag writes BAM output (stdout) rather than SAM.
    #
    # Fallback: if alnqual exits non-zero (can happen with unusual CIGAR
    # strings), copy the indelqual BAM unchanged and log a warning.
    # The downstream LoFreq call will still benefit from indelqual.
    # ----------------------------------------------------------------
    lofreq alnqual -b \\
        ${meta.id}_${meta.genotype}_iq.bam \\
        ${ref_fasta} \\
        > ${meta.id}_${meta.genotype}_iq.alnq.bam \\
    || {
        echo "WARNING: lofreq alnqual failed for ${meta.id}:${meta.genotype}; using indelqual-only BAM for downstream calling." >&2
        cp ${meta.id}_${meta.genotype}_iq.bam ${meta.id}_${meta.genotype}_iq.alnq.bam
    }

    samtools index ${meta.id}_${meta.genotype}_iq.alnq.bam

    # ----------------------------------------------------------------
    # Rename to the canonical output name consumed downstream.
    # ----------------------------------------------------------------
    mv ${meta.id}_${meta.genotype}_iq.alnq.bam     ${meta.id}_${meta.genotype}_preprocessed.bam
    mv ${meta.id}_${meta.genotype}_iq.alnq.bam.bai ${meta.id}_${meta.genotype}_preprocessed.bam.bai

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
