/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARTITION_READS — Split Round-1 BAM reads into per-genotype FASTQ files for
    mixed-genotype samples.

    Wraps `bin/partition_reads.py`. See:
        docs/data_flow.md            Step 5.9
        docs/architecture_reasoning.md Section 5
        docs/implementation_prompts.md Prompt 9

    This process is ONLY invoked when GENOTYPE_CLASSIFY flagged the sample as
    mixed (is_mixed = true in genotype_summary.json). Single-genotype samples
    skip this process entirely and pass the host-depleted FASTQ to consensus
    building directly. The conditional gate lives in the calling subworkflow
    (genotype_branch.nf, Prompt 10), not in this module.

    Inputs:
        tuple val(meta), path(bam), path(bai), path(assignments)
            meta.id is used as the output filename prefix and propagated into
            the partition_summary.json. The assignments TSV is the read-level
            output of GENOTYPE_CLASSIFY; the BAM is the Round-1 BAM that fed
            classification (we extract SEQ/QUAL directly from the BAM rather
            than re-opening the original FASTQ).
        val genotypes
            Comma-separated string of genotype labels to emit explicitly
            (e.g. "1,3"). When provided, only these genotypes get dedicated
            output FASTQs and any read assigned to a genotype outside this
            set is routed to the ambiguous pool. Pass an empty string to let
            the script derive the genotypes from the TSV.

    Outputs:
        reads     — [meta, "*.fastq.gz"]  glob-collected per-genotype + ambiguous
                                          FASTQs. Downstream subworkflow expands
                                          via flatMap so each genotype becomes
                                          its own channel item.
        summary   — [meta, "partition_summary.json"]  per-genotype counts.
        versions  — versions.yml

    Container:
        Reuse the LoFreq biocontainer (already used by GENOTYPE_CLASSIFY) — it
        bundles Python 3 + pysam, which is everything partition_reads.py needs.
        Conda fallback installs pysam directly.

    Label: process_medium (2 CPU, 4 GB, 20 min — Step 5.9 resource spec).

    Output published to:
        ${params.outdir}/${meta.id}/genotyping/partitioned/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PARTITION_READS {

    label 'process_medium'

    tag "${meta.id}"

    // Reuse the LoFreq biocontainer (Python 3 + pysam already inside).
    container 'quay.io/biocontainers/lofreq:2.1.5--py310h0dbaff4_3'
    conda 'bioconda::pysam=0.22.1 conda-forge::python=3.11'

    publishDir (
        path: { "${params.outdir}/${meta.id}/genotyping/partitioned/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(bam), path(bai), path(assignments)
    val genotypes   // comma-separated genotype list, e.g. "1,3"; "" = auto

    output:
    tuple val(meta), path("${meta.id}.*.fastq.gz"),        emit: reads
    tuple val(meta), path("partition_summary.json"),       emit: summary
    path "versions.yml",                                   emit: versions

    script:
    // Build the optional --genotypes argument only when the caller supplied
    // a non-empty list. Empty string → let the script auto-derive from the
    // TSV (safe but less deterministic across mixed-vs-non-mixed routing).
    def gt_arg = (genotypes && genotypes.toString().trim()) \
        ? "--genotypes ${genotypes.toString().trim()}" : ''
    """
    # ----------------------------------------------------------------
    # Partition Round-1 BAM reads into per-genotype FASTQs.
    # SEQ and QUAL come from the BAM (no need to re-open the original
    # host-depleted FASTQ). The script keeps BAM orientation as-is — minimap2
    # handles strand on the next mapping.
    # ----------------------------------------------------------------
    python3 ${projectDir}/bin/partition_reads.py \\
        --bam ${bam} \\
        --assignments ${assignments} \\
        --sample-id ${meta.id} \\
        --output-dir . \\
        ${gt_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """

    stub:
    """
    # Minimal stubs — one empty FASTQ per genotype label in the list,
    # plus the ambiguous pool and a zero-count summary JSON.
    gt_list="${genotypes ?: ''}"
    if [ -z "\${gt_list}" ]; then
        gt_list="1a"   # single placeholder so the output channel isn't empty
    fi

    : > genotype_array.tmp
    for gt in \$(echo "\${gt_list}" | tr ',' ' '); do
        echo -n "" | gzip > ${meta.id}.\${gt}.fastq.gz
    done
    echo -n "" | gzip > ${meta.id}.ambiguous.fastq.gz

    python3 - <<PYEOF > partition_summary.json
import json, os
genos = [g for g in "\${gt_list}".split(",") if g]
out = {
    "sample_id":       "${meta.id}",
    "genotypes": [
        {"genotype": g, "reads_written": 0, "output_file": f"${meta.id}.{g}.fastq.gz"}
        for g in genos
    ],
    "ambiguous_reads": 0,
    "ambiguous_file":  "${meta.id}.ambiguous.fastq.gz",
    "skipped": {"no_assignment": 0, "no_seq": 0, "quality_mismatch": 0},
}
print(json.dumps(out, indent=2))
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.11"
        pysam: "0.22.1"
    END_VERSIONS
    """
}
