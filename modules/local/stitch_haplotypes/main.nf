/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STITCH_HAPLOTYPES — Post-hoc stitching of DEVIDER per-window haplotypes.

    Purpose:
        Implements Step 5.20 of the data flow specification.  Wraps
        `bin/stitch_haplotypes.py` to join DEVIDER's per-window haplotype calls
        into longer-range haplotypes using read-spanning evidence from the
        haplotype-tagged BAM emitted by `devider --output-reads`.

    Scientific rationale (CLAUDE.md D11, docs/architecture_reasoning.md §9):
        DEVIDER v0.0.1 has no `--merge-windows` flag.  It emits one haplotype
        FASTA per genomic window and a single haplotype-tagged BAM.  To get
        useful longer-range haplotypes for downstream phylogenetics and
        drug-resistance analysis, we need a post-hoc stitch.  The
        stitch_haplotypes.py script handles that:

            1. Discover DEVIDER outputs via glob (filenames vary by release).
            2. Parse window coordinates and haplotype indices from FASTA
               headers.
            3. Walk the haplotype-tagged BAM, finding reads that span pairs
               of adjacent windows and recording the (hap_in_W_i, hap_in_W_{i+1})
               pairing each read implies.
            4. Threshold by `params.stitch_min_reads`: a link is accepted
               when at least that many reads support it.
            5. Greedily extend chains left-to-right, enumerating all
               supported branches rather than arbitrating between them.

        The algorithm is heuristic — see the module-level docstring of
        `bin/stitch_haplotypes.py` for the limitations.

    Graceful failure:
        The script never raises an unhandled exception.  If DEVIDER failed
        (devider.failed marker or empty output directory) it writes an empty
        FASTA and a JSON with `fallback_used: true`.  The process exits 0 in
        all cases.

    Inputs:
        meta        — val map with `id` and `genotype` fields
        devider_dir — DEVIDER output directory from the DEVIDER process
                      (emit: outdir).  May contain `devider.failed` if
                      DEVIDER exited non-zero.

    Outputs:
        haplotypes — [meta, "*_merged_haplotypes.fasta"]
                     One FASTA record per stitched chain.  Empty if DEVIDER
                     failed or no chains could be built.  Always present so
                     downstream channel-collection (sample report) is reliable.
        report     — [meta, "*_stitch_report.json"]
                     Machine-readable record of windows, junctions, links,
                     and chains.  Schema documented in the script.
        versions   — versions.yml

    Container:
        Reuse the LoFreq biocontainer (pysam + Python 3 inside) — the same
        image used by GENOTYPE_CLASSIFY and PARTITION_READS.  Avoids pulling
        an extra image for the same dependency set.

    Label: process_medium (4 CPU, 8 GB, 30 min — Step 5.20 resource spec).

    Output published to:
        ${params.outdir}/${meta.id}/haplotypes/${meta.genotype}/

    See also:
        docs/data_flow.md              Step 5.20
        docs/implementation_prompts.md Prompt 21
        docs/architecture_reasoning.md §9
        CLAUDE.md                      D11, params.stitch_min_reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process STITCH_HAPLOTYPES {

    label 'process_low'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py310h4a09ff2_0' :
        'quay.io/biocontainers/pysam:0.24.0--py310h4a09ff2_0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/haplotypes/${meta.genotype}/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(devider_dir)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_merged_haplotypes.fasta"),
          emit: haplotypes
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_stitch_report.json"),
          emit: report
    path "versions.yml", emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Stitch DEVIDER per-window haplotypes using the haplotype-tagged
    # BAM produced by `devider --output-reads`.
    #
    # The script discovers DEVIDER outputs by glob (no hard-coded
    # filenames) and writes empty FASTA + a `fallback_used: true`
    # JSON when DEVIDER failed.  It always exits 0 — see the module
    # header for the graceful-failure contract.
    # ----------------------------------------------------------------
    python3 ${projectDir}/bin/stitch_haplotypes.py \\
        --devider-dir ${devider_dir} \\
        --output-fasta ${meta.id}_${meta.genotype}_merged_haplotypes.fasta \\
        --output-json ${meta.id}_${meta.genotype}_stitch_report.json \\
        --sample-id ${meta.id} \\
        --genotype ${meta.genotype} \\
        --min-reads ${params.stitch_min_reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """

    stub:
    """
    # Empty stub outputs — keeps Nextflow output globs happy when this
    # process is executed in -stub-run mode.
    touch ${meta.id}_${meta.genotype}_merged_haplotypes.fasta

    cat <<-PYEOF > ${meta.id}_${meta.genotype}_stitch_report.json
    {
      "sample_id": "${meta.id}",
      "genotype":  "${meta.genotype}",
      "windows_found": 0,
      "haplotypes_per_window": [],
      "chains_stitched": 0,
      "chains": [],
      "unlinked_windows": [],
      "fallback_used": true,
      "reason": "stub run"
    }
    PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.11"
        pysam: "0.24.0"
    END_VERSIONS
    """
}
