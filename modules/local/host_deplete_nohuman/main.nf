/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HOST_DEPLETE_NOHUMAN — Remove host (human) reads using nohuman.

    nohuman uses a pre-built Kraken2 database of human k-mers to classify reads.
    Unclassified reads (non-human) are written to output. Faster and lower-memory
    than alignment-based depletion; well-suited for ONT clinical samples where the
    goal is contamination removal rather than sensitive host-genome analysis.

    Database: auto-downloaded by DOWNLOAD_NOHUMAN_DB (params.nohuman_db) and passed
    as input to this process, so it is shared across all samples in a run.

    Stats JSON schema (same as HOST_DEPLETE_MINIMAP2 / HOST_DEPLETE_HOSTILE):
        {
          "sample_id":          "<id>",
          "total_reads":        <int>,
          "host_reads_removed": <int>,
          "kept_reads":         <int>,
          "host_fraction":      <float>
        }

    Warning (not failure):
        If host_fraction >= 0.995 a warning is written to stderr and the sentinel
        file NO_VIRAL_READS_LIKELY is created. The process exits 0 — the pipeline
        continues flagged per design decision D4.

    Container: nohuman biocontainer (verify tag at quay.io/biocontainers/nohuman)
    Conda:     bioconda::nohuman
    Label:     process_high
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HOST_DEPLETE_NOHUMAN {

    label 'process_high'

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nohuman:0.5.1--hfa8f182_1' :
        'quay.io/biocontainers/nohuman:0.5.1--hfa8f182_1' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/reads/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(reads)
    path  nohuman_db   // directory containing Kraken2 database files (hash.k2d, opts.k2d, taxo.k2d)

    output:
    tuple val(meta), path("${meta.id}.hostdep.fastq.gz"),    emit: reads
    tuple val(meta), path("${meta.id}.host_stats.json"),     emit: stats_json
    tuple val(meta), path("${meta.id}.NO_VIRAL_READS_LIKELY"),
          optional: true,                                    emit: sentinel
    path "versions.yml",                                     emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Count input reads before depletion
    # ----------------------------------------------------------------
    total=\$(zcat ${reads} | awk 'NR%4==1' | wc -l | tr -d ' ')

    # ----------------------------------------------------------------
    # Run nohuman — outputs gzip-compressed clean reads directly
    # ----------------------------------------------------------------
    nohuman \\
        --db ${nohuman_db} \\
        -t ${task.cpus} \\
        -F g \\
        -o ${meta.id}.hostdep.fastq.gz \\
        ${reads}

    # ----------------------------------------------------------------
    # Count kept reads and emit stats JSON
    # ----------------------------------------------------------------
    kept=\$(zcat ${meta.id}.hostdep.fastq.gz | awk 'NR%4==1' | wc -l | tr -d ' ')

    python3 - <<PYEOF
import json, sys
total = int("\${total}")
kept  = int("\${kept}")
host  = total - kept
host_frac = round(host / total, 4) if total > 0 else 0.0
stats = {
    "sample_id":           "${meta.id}",
    "total_reads":         total,
    "host_reads_removed":  host,
    "kept_reads":          kept,
    "host_fraction":       host_frac,
}
with open("${meta.id}.host_stats.json", "w") as fh:
    json.dump(stats, fh, indent=2)
if host_frac >= 0.995:
    print(
        f"WARNING: {host_frac*100:.1f}% of reads classified as human for sample ${meta.id}. "
        "Very few viral reads remain. Pipeline continues — check results carefully.",
        file=sys.stderr
    )
    open("${meta.id}.NO_VIRAL_READS_LIKELY", "w").close()
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nohuman: \$(nohuman --version 2>&1 | sed 's/nohuman //')
    END_VERSIONS
    """

    stub:
    """
    echo -n "" | gzip > ${meta.id}.hostdep.fastq.gz
    echo '{"sample_id":"${meta.id}","total_reads":0,"host_reads_removed":0,"kept_reads":0,"host_fraction":0.0}' \\
        > ${meta.id}.host_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nohuman: "0.5.1"
    END_VERSIONS
    """
}
