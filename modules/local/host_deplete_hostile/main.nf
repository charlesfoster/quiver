/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HOST_DEPLETE_HOSTILE — Remove host (human) reads using hostile.

    Alternative host depletion strategy (opt-in: params.use_hostile = true).

    hostile strategy:
        hostile clean --index <index> --fastq <reads> --out-dir hostile_out/
        Hostile renames outputs; we rename to canonical names.

    Stats JSON schema:
        {
          "sample_id": "<id>",
          "total_reads": <int>,
          "host_reads_removed": <int>,
          "kept_reads": <int>,
          "host_fraction": <float>
        }

    Warning (not failure):
        If host_fraction >= 0.995 a warning is written to stderr and a flag
        NO_VIRAL_READS_LIKELY is written alongside the outputs. The process
        still exits 0 — the pipeline continues flagged as per design decision D4.

    Container: hostile 1.1.0
    Label: process_high_memory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HOST_DEPLETE_HOSTILE {

    label 'process_high_memory'

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hostile:1.1.0--pyhdfd78af_0' :
        'quay.io/biocontainers/hostile:1.1.0--pyhdfd78af_0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/reads/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(reads)
    path  host_index_dir   // directory containing hostile-compatible index files

    output:
    tuple val(meta), path("${meta.id}.hostdep.fastq.gz"),    emit: reads
    tuple val(meta), path("${meta.id}.host_stats.json"),     emit: stats_json
    tuple val(meta), path("${meta.id}.NO_VIRAL_READS_LIKELY"),
          optional: true,                                    emit: sentinel
    path "versions.yml",                                     emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Run hostile clean
    # ----------------------------------------------------------------
    mkdir -p hostile_out

    hostile clean \\
        --fastq ${reads} \\
        --index ${host_index_dir} \\
        --out-dir hostile_out \\
        --threads ${task.cpus}

    # hostile names the cleaned file <input_stem>.clean.fastq.gz;
    # rename to our canonical output name.
    cleaned=\$(ls hostile_out/*.clean.fastq.gz | head -1)
    mv "\${cleaned}" ${meta.id}.hostdep.fastq.gz

    # ----------------------------------------------------------------
    # Extract stats from hostile JSON log (if present) or recount
    # ----------------------------------------------------------------
    hostile_json=\$(ls hostile_out/*.json 2>/dev/null | head -1 || true)
    if [ -n "\${hostile_json}" ]; then
        python3 - <<PYEOF
import json, sys
with open("\${hostile_json}") as fh:
    d = json.load(fh)
# hostile JSON structure: list with one item per sample
entry = d[0] if isinstance(d, list) else d
total  = entry.get("reads_in", 0)
kept   = entry.get("reads_out", 0)
host   = total - kept
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
        f"WARNING: {host_frac*100:.1f}% of reads mapped to the host reference for sample ${meta.id}. "
        "Very few viral reads remain. Pipeline continues — check results carefully.",
        file=sys.stderr
    )
    open("${meta.id}.NO_VIRAL_READS_LIKELY", "w").close()
PYEOF
    else
        # Fall back to counting reads directly
        total=\$(zcat ${reads} | awk 'NR%4==1' | wc -l | tr -d ' ')
        kept=\$(zcat  ${meta.id}.hostdep.fastq.gz | awk 'NR%4==1' | wc -l | tr -d ' ')
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
        f"WARNING: {host_frac*100:.1f}% of reads for sample ${meta.id} mapped to the host. "
        "Very few viral reads remain.",
        file=sys.stderr
    )
    open("${meta.id}.NO_VIRAL_READS_LIKELY", "w").close()
PYEOF
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version 2>&1 | sed 's/hostile //')
    END_VERSIONS
    """

    stub:
    """
    echo -n "" | gzip > ${meta.id}.hostdep.fastq.gz
    echo '{"sample_id":"${meta.id}","total_reads":0,"host_reads_removed":0,"kept_reads":0,"host_fraction":0.0}' \\
        > ${meta.id}.host_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: "1.1.0"
    END_VERSIONS
    """
}
