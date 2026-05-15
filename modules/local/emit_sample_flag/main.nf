/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EMIT_SAMPLE_FLAG — Collect pipeline sentinel files and write a summary flags.txt.

    Sentinels emitted by other processes:
        ${meta.id}.EMPTY_INPUT        — FASTQ was 0 bytes at input (INPUT_CHECK)
        ${meta.id}.ALL_READS_FILTERED — all reads removed by chopper (CHOPPER)
        ${meta.id}.NO_VIRAL_READS_LIKELY — ≥99.5% host reads (HOST_DEPLETE_*)

    This process collects all sentinels for a single sample and writes them
    to a human-readable flags.txt under the sample QC directory.

    Output published to:
        ${params.outdir}/${meta.id}/qc/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EMIT_SAMPLE_FLAG {

    label 'process_low'

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://python:3.11' :
        'python:3.11' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/qc/" },
        mode: 'copy'
    )

    input:
    // Each sentinel is an optional file; use path with optional: true at call site.
    // We receive a list of sentinel paths collected for this sample.
    tuple val(meta), path(sentinels)

    output:
    tuple val(meta), path("flags.txt"), emit: flags

    script:
    """
    python3 - <<'PYEOF'
import os, glob

sentinel_files = glob.glob("*.EMPTY_INPUT") + \\
                 glob.glob("*.ALL_READS_FILTERED") + \\
                 glob.glob("*.NO_VIRAL_READS_LIKELY")

sentinel_meaning = {
    "EMPTY_INPUT":           "Input FASTQ was empty (0 bytes) — sample skipped.",
    "ALL_READS_FILTERED":    "All reads removed by chopper quality/length filter — sample skipped.",
    "NO_VIRAL_READS_LIKELY": "≥99.5% of reads mapped to host genome — very few viral reads remain.",
}

lines = [f"Sample: ${meta.id}", "=" * 60]
if not sentinel_files:
    lines.append("No pipeline flags raised for this sample.")
else:
    lines.append(f"{len(sentinel_files)} flag(s) raised:")
    for sf in sorted(sentinel_files):
        flag_name = sf.split(".", 1)[1] if "." in sf else sf
        meaning   = sentinel_meaning.get(flag_name, "Unknown flag.")
        lines.append(f"  [{flag_name}] {meaning}")

lines.append("")

with open("flags.txt", "w") as fh:
    fh.write("\\n".join(lines))
PYEOF
    """

    stub:
    """
    echo "Sample: ${meta.id}" > flags.txt
    echo "No pipeline flags raised (stub)." >> flags.txt
    """
}
