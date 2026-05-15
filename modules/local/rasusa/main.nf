/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RASUSA — Depth-capped random subsampling of FASTQ reads.

    Purpose:
        Generic, reusable subsampling module.  Used by PREP_LOFREQ_INPUT (Prompt 15)
        with params.lofreq_max_depth and params.rasusa_seed_lofreq, and by
        PREP_DEVIDER_INPUT (Prompt 19) with params.devider_max_depth and
        params.rasusa_seed_devider.

    Scientific rationale (docs/architecture_reasoning.md §8, CLAUDE.md D10):
        LoFreq sensitivity plateaus at ~1,000–5,000×; beyond ~10,000× error
        stratification increases the false-positive rate.  Cap: 5,000× (default).
        DEVIDER memory scales super-linearly with depth; benchmarks use 100–2,000×.
        Cap: 1,000× (default).
        The seed parameter ensures reproducibility across pipeline re-runs.

    rasusa flag notes:
        --coverage      Target coverage depth.
        --genome-size   Reference genome size in bp.
        --seed / -s     Random seed for reproducibility.
        -o              Output FASTQ (gzipped when extension is .gz).
        <FILE>          Input FASTQ is a positional argument (no -i flag).

    If the input depth is already below the cap, rasusa passes all reads through
    unchanged (it does not up-sample).  The output is always written to a new file,
    so the full-depth FASTQ is preserved.

    This module is parameterised via the `coverage`, `genome_size`, and `seed`
    input channels so callers can pass any cap/seed combination without duplicating
    the process definition.

    Inputs:
        meta         — val map with `id` and `genotype` fields
        reads        — per-genotype FASTQ to subsample
        coverage     — integer coverage cap (e.g. 5000 or 1000)
        genome_size  — integer reference genome size in bp (default HCV: 9646)
        seed         — integer random seed

    Outputs:
        reads        — [meta, "*_subsampled.fastq.gz"]
        versions     — versions.yml

    Container: quay.io/biocontainers/rasusa:2.1.0--hc1c3326_1
    Label: process_medium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RASUSA {

    label 'process_low'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rasusa:2.1.0--hc1c3326_1' :
        'quay.io/biocontainers/rasusa:2.1.0--hc1c3326_1' }"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads), val(coverage), val(genome_size), val(seed)

    output:
    tuple val(meta), path("${meta.id}_${meta.genotype}_subsampled.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Subsample reads to at most <coverage>× of a <genome_size> bp genome.
    #
    # Flag rationale:
    #   --coverage     target depth (rasusa will down-sample to this; if
    #                  input is already below, all reads pass through)
    #   --genome-size  length of the reference in bp (9646 for HCV)
    #   --seed         random seed for reproducibility
    #   -o             output file (.gz extension → gzip compressed)
    #   <FILE>         input is positional (no -i flag in rasusa 2.x)
    # ----------------------------------------------------------------
    rasusa reads \\
        --coverage ${coverage} \\
        --genome-size ${genome_size} \\
        --seed ${seed} \\
        -o ${meta.id}_${meta.genotype}_subsampled.fastq.gz \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: \$(rasusa --version 2>&1 | sed 's/rasusa //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.genotype}_subsampled.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: "2.1.0"
    END_VERSIONS
    """
}
