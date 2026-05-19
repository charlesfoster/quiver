/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOFREQ_CALL — Low-frequency variant calling with LoFreq.

    Purpose:
        Implements Step 5.15 of the data flow specification.  Calls SNVs and indels
        in the preprocessed, depth-capped BAM using LoFreq.  Parallel calling is
        enabled only when params.lofreq_pp_threads > 1.

    Scientific rationale (CLAUDE.md D9, docs/architecture_reasoning.md §7):
        LoFreq is the primary variant caller for this pipeline.  It was designed for
        low-allele-frequency detection with rigorous Bonferroni/FDR control.  The
        ONT preprocessing chain from LOFREQ_PREPROCESS adds indel qualities
        so LoFreq can use its indel model on ONT reads.
        See D9 for the full justification.

    Key flags:
        lofreq call     Serial caller. Used when params.lofreq_pp_threads <= 1,
                        including the docker_mac profile.
        call-parallel   Multi-process caller. Used when params.lofreq_pp_threads > 1.
        --pp-threads N  Capped at min(task.cpus, params.lofreq_pp_threads). Note:
                        --threads is a different flag controlling pre-processing;
                        --pp-threads controls the parallel calling pool.
        --call-indels   Enables indel calling.  Requires BI/BD tags in the BAM
                        (provided by lofreq indelqual --dindel in LOFREQ_PREPROCESS).
                        Must be explicit — LoFreq defaults to SNV-only.
        --min-mq        Minimum mapping quality (params.min_mq, default 20).
        --min-bq        Minimum base quality for all bases (params.min_bq, default 7).
        --min-alt-bq    Minimum base quality for alternate-allele bases
                        (params.min_alt_bq, default 7).  LoFreq's built-in default
                        is 6, which is below our --min-bq of 7 and produces a fatal
                        error ("min base-call quality for all bases larger than min
                        base-call quality for alternate bases").  Must satisfy
                        min_bq <= min_alt_bq; pipeline startup validates this.
        --sig           Bonferroni p-value significance threshold (params.lofreq_sig,
                        default 0.01).
        --min-cov 1     Minimum site coverage to attempt calling.  The downstream
                        VARIANT_FILTER applies the final depth gate (params.min_variant_depth,
                        default 20), so the raw calls are unrestricted here to allow
                        the TSV to report all attempted sites.
        -f              Reference FASTA (per-genotype consensus).
        -o              Output VCF (plain text; bgzip + tabix follow).

    Docker memory note (conf/modules.config):
        This process has `memory = null` in modules.config, which prevents Nextflow
        from passing --memory to Docker.  Without this override, Nextflow emits
        --memory 32768m (from process_high), which sets a cgroup limit exceeding the
        Docker Desktop VM's total RAM on macOS Apple Silicon.  The cgroup OOM killer
        then sends SIGKILL to lofreq worker subprocesses (exit 137).  The fix is to
        remove the --memory cap so the workers compete fairly for available VM memory.
        On native Linux HPC nodes the memory directive poses no problem; HPC profiles
        can re-add memory accounting via their own withName overrides if needed.

    Post-call compression:
        bgzip + tabix produce the standard .vcf.gz + .tbi pair consumed by bcftools
        downstream (VARIANT_FILTER, Prompt 17) and optionally DEVIDER (Prompt 20).

    Inputs:
        meta        — val map with `id` and `genotype` fields
        bam         — preprocessed BAM from LOFREQ_PREPROCESS (indelqual)
        bai         — BAM index
        ref_fasta   — per-genotype consensus FASTA

    Outputs:
        vcf         — [meta, "*.vcf.gz", "*.vcf.gz.tbi"]
        versions    — versions.yml

    Container: quay.io/biocontainers/lofreq:2.1.5--py310h4966b78_15
    Label: process_high (16 CPU by default, memory = null via modules.config, 12 h).
           docker_mac overrides this to 1 CPU / maxForks 1.

    Output published to:
        ${params.outdir}/${meta.id}/variants/${meta.genotype}/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process LOFREQ_CALL {

    label 'process_high'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py310h4966b78_15' :
        'quay.io/biocontainers/lofreq:2.1.5--py310h4966b78_15' }"
    conda "${moduleDir}/environment.yml"

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
    # Copy reference so .fai lives alongside the FASTA in the work directory.
    cp ${ref_fasta} ref.fasta
    samtools faidx ref.fasta
    REF_ABS=\$(pwd)/ref.fasta
    BAM_ABS=\$(readlink -f ${bam})

    LOFREQ_PP_THREADS=\$(( ${params.lofreq_pp_threads} < ${task.cpus} ? ${params.lofreq_pp_threads} : ${task.cpus} ))

    BAQ_FLAG=""
    ${ params.lofreq_no_baq ? 'BAQ_FLAG="-B"' : '' }

    if [ "\${LOFREQ_PP_THREADS}" -gt 1 ]; then
        lofreq call-parallel \\
            --pp-threads "\${LOFREQ_PP_THREADS}" \\
            --call-indels \\
            \${BAQ_FLAG} \\
            --min-mq ${params.min_mq} \\
            --min-bq ${params.min_bq} \\
            --min-alt-bq ${params.min_alt_bq} \\
            --sig ${params.lofreq_sig} \\
            --min-cov 1 \\
            -f "\${REF_ABS}" \\
            -o ${meta.id}_${meta.genotype}_lofreq.vcf \\
            "\${BAM_ABS}"
    else
        lofreq call \\
            --call-indels \\
            \${BAQ_FLAG} \\
            --min-mq ${params.min_mq} \\
            --min-bq ${params.min_bq} \\
            --min-alt-bq ${params.min_alt_bq} \\
            --sig ${params.lofreq_sig} \\
            --min-cov 1 \\
            -f "\${REF_ABS}" \\
            -o ${meta.id}_${meta.genotype}_lofreq.vcf \\
            "\${BAM_ABS}"
    fi

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
