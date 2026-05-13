/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOFREQ_CALL — Low-frequency variant calling with LoFreq.

    Purpose:
        Implements Step 5.15 of the data flow specification.  Calls SNVs and indels
        in the preprocessed, depth-capped BAM using LoFreq's parallel caller.

    Scientific rationale (CLAUDE.md D9, docs/architecture_reasoning.md §7):
        LoFreq is the primary variant caller for this pipeline.  It was designed for
        low-allele-frequency detection with rigorous Bonferroni/FDR control.  The
        ONT preprocessing chain (indelqual + alnqual) from LOFREQ_PREPROCESS closes
        the gap between LoFreq's Illumina-oriented model and ONT error profiles.
        See D9 for the full justification.

    Key flags:
        call-parallel   Multi-threaded caller (threads via --pp-threads, not --threads).
        --pp-threads N  Binds to task.cpus.  Note: --threads is a different flag that
                        controls pre-processing threads; --pp-threads controls the
                        parallel calling pool.
        --call-indels   Enables indel calling.  Requires BI/BD tags in the BAM
                        (provided by lofreq indelqual --dindel in LOFREQ_PREPROCESS).
                        Must be explicit — LoFreq defaults to SNV-only.
        --min-mq        Minimum mapping quality (params.min_mq, default 20).
        --min-bq        Minimum base quality (params.min_bq, default 7).
        --sig           Bonferroni p-value significance threshold (params.lofreq_sig,
                        default 0.01).
        --min-cov 1     Minimum site coverage to attempt calling.  The downstream
                        VARIANT_FILTER applies the final depth gate (params.min_variant_depth,
                        default 20), so the raw calls are unrestricted here to allow
                        the TSV to report all attempted sites.
        -f              Reference FASTA (per-genotype consensus).
        -o              Output VCF (plain text; bgzip + tabix follow).

    Post-call compression:
        bgzip + tabix produce the standard .vcf.gz + .tbi pair consumed by bcftools
        downstream (VARIANT_FILTER, Prompt 17) and optionally DEVIDER (Prompt 20).

    Inputs:
        meta        — val map with `id` and `genotype` fields
        bam         — preprocessed BAM from LOFREQ_PREPROCESS (indelqual + alnqual)
        bai         — BAM index
        ref_fasta   — per-genotype consensus FASTA

    Outputs:
        vcf         — [meta, "*.vcf.gz", "*.vcf.gz.tbi"]
        versions    — versions.yml

    Container: quay.io/biocontainers/lofreq:2.1.5--py310h0dbaff4_3
    Label: process_high (16 CPU, 16 GB, 90 min — Step 5.15 resource spec).

    Output published to:
        ${params.outdir}/${meta.id}/variants/${meta.genotype}/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process LOFREQ_CALL {

    label 'process_high'

    tag "${meta.id}:${meta.genotype}"

    container 'quay.io/biocontainers/lofreq:2.1.5--py310h0dbaff4_3'
    conda 'bioconda::lofreq=2.1.5'

    publishDir (
        path: { "${params.outdir}/${meta.id}/variants/${meta.genotype}/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(bam), path(bai), path(ref_fasta)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_lofreq.vcf.gz"),
          path("${meta.id}_${meta.genotype}_lofreq.vcf.gz.tbi"),
          emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Run LoFreq parallel variant calling.
    #
    # --pp-threads N  : parallel calling threads — bind to task.cpus.
    #                   NOT --threads (that controls pre-processing).
    # --call-indels   : enable indel calling (requires BI/BD tags from
    #                   lofreq indelqual --dindel; mandatory per D9).
    # --min-mq        : minimum mapping quality filter.
    # --min-bq        : minimum base quality filter.
    # --sig           : Bonferroni significance threshold.
    # --min-cov 1     : accept any site with ≥1 read; downstream
    #                   VARIANT_FILTER applies the final DP gate.
    # -f              : reference FASTA (per-genotype consensus).
    # -o              : output VCF (plain text, compressed below).
    # ----------------------------------------------------------------
    lofreq call-parallel \\
        --pp-threads \$(( ${task.cpus} < 8 ? ${task.cpus} : 8 )) \\
        --call-indels \\
        --min-mq ${params.min_mq} \\
        --min-bq ${params.min_bq} \\
        --sig ${params.lofreq_sig} \\
        --min-cov 1 \\
        -f ${ref_fasta} \\
        -o ${meta.id}_${meta.genotype}_lofreq.vcf \\
        ${bam}

    # ----------------------------------------------------------------
    # Compress and index — produces the canonical .vcf.gz + .tbi pair
    # consumed by VARIANT_FILTER and optionally DEVIDER.
    # ----------------------------------------------------------------
    bgzip ${meta.id}_${meta.genotype}_lofreq.vcf
    tabix -p vcf ${meta.id}_${meta.genotype}_lofreq.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(lofreq version 2>&1 | head -1 | sed 's/lofreq version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.genotype}_lofreq.vcf.gz
    touch ${meta.id}_${meta.genotype}_lofreq.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: "2.1.5"
    END_VERSIONS
    """
}
