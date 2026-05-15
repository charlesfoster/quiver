/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MOSDEPTH — Coverage quality control on the full-depth Round 2 BAM.

    Purpose:
        Implements Step 5.13 of the data flow specification.  Runs mosdepth on the
        full-depth Round 2 BAM (before any rasusa subsampling) to produce per-window
        and summary coverage statistics.  Emits a LOW_COVERAGE sentinel if mean
        coverage falls below params.min_mean_coverage (default 100×).

    Why full-depth BAM, not the subsampled LoFreq BAM:
        Coverage QC should reflect the true sequencing depth, not the
        artificially capped depth.  The rasusa-subsampled BAMs are used only for
        variant calling and haplotype reconstruction; the coverage report should
        describe what was actually sequenced.  See docs/architecture_reasoning.md §8
        and CLAUDE.md D10.

    mosdepth flags:
        --fast-mode     Uses HTS pileup instead of a per-base array; ~2× faster.
        --no-per-base   Suppresses the per-base depth BED (large; not needed here).
        --by 100        Emits per-window coverage in 100 bp bins; used by the
                        per-sample report to plot regional coverage variation.
        -t              Thread count (process_medium: 4 threads).

    LOW_COVERAGE sentinel:
        The last line of the mosdepth summary file (*.mosdepth.summary.txt) is the
        "total_region" row for the whole-genome mean.  Column 4 (0-indexed col 3)
        is the mean coverage.  If mean < params.min_mean_coverage, a sentinel file
        is created.  The process exits 0 — the sentinel is routed to the per-sample
        report and DEVIDER is skipped for this branch (Prompt 24).

    Inputs:
        meta  — val map with `id` and `genotype` fields
        bam   — sorted, indexed Round 2 BAM (full-depth)
        bai   — BAM index

    Outputs:
        summary      — [meta, "*.mosdepth.summary.txt"]
        regions      — [meta, "*.regions.bed.gz"]
        low_cov_flag — [meta, "*.LOW_COVERAGE"]  optional
        versions     — versions.yml

    Container: quay.io/biocontainers/mosdepth:0.3.10--h4e814b3_1
    Label: process_medium (4 CPU, 4 GB, 10 min — Step 5.13 resource spec).

    Output published to:
        ${params.outdir}/${meta.id}/mapping/${meta.genotype}/coverage/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MOSDEPTH {

    label 'process_medium'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.10--h4e814b3_1' :
        'quay.io/biocontainers/mosdepth:0.3.10--h4e814b3_1' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/mapping/${meta.genotype}/coverage/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}.mosdepth.summary.txt"),
          emit: summary
    tuple val(meta),
          path("${meta.id}_${meta.genotype}.regions.bed.gz"),
          emit: regions
    tuple val(meta),
          path("${meta.id}_${meta.genotype}.LOW_COVERAGE"),
          optional: true,
          emit: low_cov_flag
    path "versions.yml", emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Run mosdepth for coverage QC.
    #
    # --fast-mode     : HTS pileup, ~2× faster than the default.
    # --no-per-base   : Skip the per-base depth BED; saves disk space.
    # --by 100        : Emit per-window coverage in 100 bp bins.
    # -t              : Thread count.
    #
    # mosdepth names outputs by the prefix supplied as the positional arg.
    # Using "${meta.id}_${meta.genotype}" gives stable, unique output names
    # that match the expected patterns in the output declarations above.
    # ----------------------------------------------------------------
    mosdepth \\
        --fast-mode \\
        --no-per-base \\
        --by 100 \\
        -t ${task.cpus} \\
        ${meta.id}_${meta.genotype} \\
        ${bam}

    # ----------------------------------------------------------------
    # LOW_COVERAGE sentinel.
    #
    # The last non-empty line of the summary file is:
    #   total_region  <chrom>  <len>  <mean_cov>  <min>  <max>
    # Column 4 (1-indexed) is the mean coverage over all windows.
    #
    # A Python one-liner compares the float against params.min_mean_coverage.
    # If below threshold: create the sentinel and log a warning.
    # The process exits 0 regardless.
    # ----------------------------------------------------------------
    mean_cov=\$(grep "total_region" ${meta.id}_${meta.genotype}.mosdepth.summary.txt \\
               | awk '{print \$4}')

    if python3 -c "import sys; sys.exit(0 if float('\${mean_cov}') < ${params.min_mean_coverage} else 1)" 2>/dev/null; then
        echo "WARNING: Mean coverage \${mean_cov}× is below threshold ${params.min_mean_coverage}× for sample ${meta.id} genotype ${meta.genotype}. Emitting LOW_COVERAGE sentinel." >&2
        touch ${meta.id}_${meta.genotype}.LOW_COVERAGE
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """

    stub:
    """
    printf "chrom\\tstart\\tend\\tmean_cov\\tmin_cov\\tmax_cov\\n" \\
        > ${meta.id}_${meta.genotype}.mosdepth.summary.txt
    printf "total_region\\t0\\t9646\\t500.0\\t0\\t10000\\n" \\
        >> ${meta.id}_${meta.genotype}.mosdepth.summary.txt
    touch ${meta.id}_${meta.genotype}.regions.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: "0.3.10"
    END_VERSIONS
    """
}
