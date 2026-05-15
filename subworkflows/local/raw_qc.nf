/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RAW_QC — Run NanoPlot and nanoq in parallel on each raw sample.

    Both tools receive the same input channel and are executed concurrently by
    Nextflow's task scheduler. Outputs are collected and emitted jointly so the
    calling workflow can pass them to MultiQC.

    Input:
        reads — channel of [meta, fastq] tuples

    Emit:
        nanoplot_dirs — channel of [meta, nanoplot_dir] tuples
        nanoq_jsons   — channel of [meta, nanoq_json] tuples
        versions      — channel of versions.yml files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NANOPLOT } from '../../modules/local/nanoplot/main'
include { NANOQ    } from '../../modules/local/nanoq/main'

workflow RAW_QC {

    take:
    reads  // channel: [meta, fastq]

    main:
    // Run both tools in parallel on the same input channel.
    // Nextflow forks the channel automatically when it is consumed by two
    // different processes — no explicit branch needed.

    NANOPLOT ( reads )
    NANOQ    ( reads )

    // Collect version files from all tools
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(NANOPLOT.out.versions)
    ch_versions = ch_versions.mix(NANOQ.out.versions)

    emit:
    nanoplot_dirs = NANOPLOT.out.nanoplot_dir  // [meta, dir]
    nanoq_jsons   = NANOQ.out.json             // [meta, json]
    versions      = ch_versions
}
