/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOTYPE_CLASSIFY — Genotype + mixed-infection classification from
    Round 1 competitive-mapping BAM.

    Wraps `bin/classify_genotype.py`. See:
        docs/data_flow.md            Step 5.8
        docs/architecture_reasoning.md Section 5
        docs/configuration.md        (genotype_summary.json schema)

    Inputs:
        tuple val(meta), path(bam), path(bai)
            meta.id is propagated into the JSON `sample_id` field and used as
            the output prefix.

    Outputs:
        summary       — [meta, "${meta.id}.genotype_summary.json"]
                        consumed by the genotype-branching subworkflow
                        (Prompt 10) and the per-sample report (Prompt 22).
        assignments   — [meta, "${meta.id}.read_assignments.tsv"]
                        one row per primary mapped read; useful for QC and
                        for downstream PARTITION_READS (Prompt 9).
        versions      — versions.yml

    Container:
        The LoFreq biocontainer is the only image in the tool inventory
        (CLAUDE.md §4) that ships both Python 3 and pysam.  Re-using it here
        avoids introducing another image whose only job is to host pysam.

    Label: process_medium (2 CPU, 4 GB, 15 min — Step 5.8 resource spec).

    Output published to:
        ${params.outdir}/${meta.id}/genotyping/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GENOTYPE_CLASSIFY {

    label 'process_medium'

    tag "${meta.id}"

    // The LoFreq biocontainer bundles Python 3 + pysam; reuse it here to keep
    // the container manifest small.  Conda fallback installs pysam directly.
    container 'quay.io/biocontainers/lofreq:2.1.5--py310h0dbaff4_3'
    conda 'bioconda::pysam=0.22.1 conda-forge::python=3.11'

    publishDir (
        path: { "${params.outdir}/${meta.id}/genotyping/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.genotype_summary.json"),  emit: summary
    tuple val(meta), path("${meta.id}.read_assignments.tsv"),   emit: assignments
    path "versions.yml",                                        emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Classify reads, detect mixed infection, write summary + TSV.
    # The script writes structured warnings to stderr; exit code is 2
    # if the BAM is genuinely empty (no primary alignments at all).
    # ----------------------------------------------------------------
    python3 ${projectDir}/bin/classify_genotype.py \\
        --bam ${bam} \\
        --sample-id ${meta.id} \\
        --min-secondary-fraction ${params.min_secondary_fraction} \\
        --ambiguous-delta-as ${params.ambiguous_delta_as} \\
        --out-tsv ${meta.id}.read_assignments.tsv \\
        --out-json ${meta.id}.genotype_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """

    stub:
    """
    # Minimal valid JSON matching the schema in docs/configuration.md
    cat > ${meta.id}.genotype_summary.json <<'JSON'
    {
      "sample_id": "${meta.id}",
      "total_mapped_reads": 0,
      "ambiguous_reads": 0,
      "ambiguous_fraction": 0.0,
      "genotypes": [],
      "is_mixed": false,
      "primary_genotype": null,
      "secondary_genotypes": [],
      "branches_to_run": []
    }
    JSON
    printf 'read_name\\treference\\tsubtype\\tgenotype\\tAS\\tXS\\tis_ambiguous\\n' \\
        > ${meta.id}.read_assignments.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.11"
        pysam: "0.22.1"
    END_VERSIONS
    """
}
