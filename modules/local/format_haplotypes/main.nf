/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FORMAT_HAPLOTYPES — Sort, annotate, and publish DEVIDER haplotype output.

    Purpose:
        Wraps bin/format_haplotypes.py.  Receives the DEVIDER output directory,
        re-sorts haplotypes by abundance (highest first), assigns clean sequential
        IDs, and writes enriched FASTA headers:

            >S01_1_haplotype_0 abund:31.98 depth:1392.79 length:9286

        Also writes a header mapping TSV so downstream tools can cross-reference
        new IDs back to DEVIDER's original headers.

    Why ALL-ALL is the normal case:
        In this pipeline, DEVIDER always receives reads >= 2,000 bp and a
        LoFreq-filtered VCF with ~26 SNPs across the 9,286 bp HCV genome.  With
        reads spanning the full genome and a sparse SNP set, DEVIDER reconstructs
        all haplotypes in a single global window (Range:ALL-ALL).  The script
        retains its window-stitching logic as a safe fallback for unusual samples.

    Graceful failure:
        If DEVIDER emitted devider.failed or produced no FASTAs, the script writes
        an empty FASTA, an empty mapping TSV, and a JSON with fallback_used: true.
        The process always exits 0.

    Inputs:
        meta        — val map with `id` and `genotype` fields
        devider_dir — DEVIDER output directory (from DEVIDER process emit: outdir)

    Outputs:
        haplotypes  — [meta, "*_haplotypes.fasta"]     — one record per haplotype
        mapping     — [meta, "*_haplotype_map.tsv"]    — new ID → original header map
        report      — [meta, "*_haplotype_report.json"]— machine-readable summary
        versions    — versions.yml

    Container:
        Reuse the pysam biocontainer (pysam + Python 3).

    Label: process_low (header rewriting; no compute-intensive work)

    Output published to:
        ${params.outdir}/${meta.id}/haplotypes/${meta.genotype}/

    See also:
        docs/data_flow.md              Step 5.20
        CLAUDE.md                      D11
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FORMAT_HAPLOTYPES {

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
          path("${meta.id}_${meta.genotype}${meta.devider_reference == 'round1' ? '_round1_comparison' : ''}_haplotypes.fasta"),
          emit: haplotypes
    tuple val(meta),
          path("${meta.id}_${meta.genotype}${meta.devider_reference == 'round1' ? '_round1_comparison' : ''}_haplotype_map.tsv"),
          emit: mapping
    tuple val(meta),
          path("${meta.id}_${meta.genotype}${meta.devider_reference == 'round1' ? '_round1_comparison' : ''}_haplotype_report.json"),
          emit: report
    path "versions.yml", emit: versions

    script:
    def output_prefix = "${meta.id}_${meta.genotype}${meta.devider_reference == 'round1' ? '_round1_comparison' : ''}"
    def reference_arm = meta.devider_reference ?: 'round2_simple'
    def reference_label = reference_arm == 'round1' ? 'Round 1 consensus' : 'Round 2 simple consensus'
    def primary_arg = reference_arm == 'round2_simple' ? '--primary' : ''
    def id_suffix = reference_arm == 'round1' ? '_round1_comparison' : ''
    """
    python3 ${projectDir}/bin/format_haplotypes.py \\
        --devider-dir ${devider_dir} \\
        --output-fasta ${output_prefix}_haplotypes.fasta \\
        --output-mapping ${output_prefix}_haplotype_map.tsv \\
        --output-json ${output_prefix}_haplotype_report.json \\
        --sample-id ${meta.id} \\
        --genotype ${meta.genotype} \\
        --reference-arm ${reference_arm} \\
        --reference-label '${reference_label}' \\
        --id-suffix '${id_suffix}' \\
        ${primary_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """

    stub:
    def output_prefix_s = "${meta.id}_${meta.genotype}${meta.devider_reference == 'round1' ? '_round1_comparison' : ''}"
    def reference_arm_s = meta.devider_reference ?: 'round2_simple'
    def reference_label_s = reference_arm_s == 'round1' ? 'Round 1 consensus' : 'Round 2 simple consensus'
    def primary_s = reference_arm_s == 'round2_simple' ? 'true' : 'false'
    """
    touch ${output_prefix_s}_haplotypes.fasta
    echo -e "new_id\tabundance_pct\tdepth\tlength_bp\toriginal_devider_header" \\
        > ${output_prefix_s}_haplotype_map.tsv

    cat <<-PYEOF > ${output_prefix_s}_haplotype_report.json
    {
      "sample_id": "${meta.id}",
      "genotype":  "${meta.genotype}",
      "reference_arm": "${reference_arm_s}",
      "reference_label": "${reference_label_s}",
      "primary": ${primary_s},
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
