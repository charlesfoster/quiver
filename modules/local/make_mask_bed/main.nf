/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAKE_MASK_BED — Coverage mask generation using mosdepth quantize mode.

    Runs mosdepth in quantize mode to produce a BED of regions below
    params.min_consensus_cov.  These positions will be written as N
    in the final consensus FASTA (design decision D8).

    mosdepth --quantize 0:<threshold>:
        Produces a file with rows like:
            chr1  0    100  0:10      <- 0–10 coverage (below threshold)
            chr1  100  200  10:500    <- 10–500 coverage (above threshold)
        We keep rows where the label is "0:<threshold>" (the below-threshold bin).
        With threshold=10 (default) the below-threshold label is "0:10".

    Part of the BUILD_CONSENSUS subworkflow (Step 5.10).

    Container: mosdepth 0.3.10
    Label: process_low
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MAKE_MASK_BED {

    label 'process_low'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.10--h4e814b3_1' :
        'quay.io/biocontainers/mosdepth:0.3.10--h4e814b3_1' }"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("mask.bed"), emit: mask_bed
    path "versions.yml",               emit: versions

    script:
    // Build the quantize string: "0:<threshold>:" — three breakpoints that
    // create three bins: [0,1), [1,threshold), [threshold,inf).
    // We want to mask any position in the first two bins.
    def threshold = params.min_consensus_cov
    """
    # ----------------------------------------------------------------
    # Run mosdepth in quantize mode.
    # --no-per-base suppresses the per-base .bed.gz (not needed here).
    # The quantize string "0:${threshold}:" creates bins:
    #   [0, ${threshold})  — low coverage (mask these)
    #   [${threshold}, ∞)  — adequate coverage
    # Note: mosdepth interprets the string as right-exclusive breakpoints
    # so "0:${threshold}:" means: 0 ≤ cov < ${threshold} and ${threshold} ≤ cov.
    # ----------------------------------------------------------------
    mosdepth \\
        --quantize 0:${threshold}: \\
        --no-per-base \\
        ${meta.id}_${meta.genotype}_cov \\
        ${bam}

    # ----------------------------------------------------------------
    # Extract positions in the below-threshold bin.
    # The quantized BED has column 4 = bin label, e.g. "0:10" or "10:500".
    # We keep rows whose label indicates coverage < threshold.
    # With "0:<threshold>:" the below-threshold label is "0:<threshold>".
    # ----------------------------------------------------------------
    LOW_BIN="0:${threshold}"

    zcat ${meta.id}_${meta.genotype}_cov.quantized.bed.gz \\
        | awk -v low_bin="\${LOW_BIN}" '\$4 == low_bin { print \$1"\\t"\$2"\\t"\$3 }' \\
        > mask.bed

    # mask.bed may be empty — that is fine; bcftools consensus accepts an empty -m.
    echo "Mask BED generated: \$(wc -l < mask.bed) regions below coverage ${threshold}." >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """

    stub:
    """
    touch mask.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: "0.3.10"
    END_VERSIONS
    """
}
