/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HOST_DEPLETE — Remove host (human) reads by mapping to GRCh38 no-alt.

    Two process variants are defined:
      1. HOST_DEPLETE_MINIMAP2  — default; uses minimap2 + samtools (label: process_high_memory)
      2. HOST_DEPLETE_HOSTILE   — alternative; uses hostile (label: process_high_memory)

    The calling workflow selects between them via:
        if (params.use_hostile) { HOST_DEPLETE_HOSTILE(...) }
        else                    { HOST_DEPLETE_MINIMAP2(...) }

    Both processes emit identical output channels so downstream modules are
    agnostic to which deplete strategy was used:
        reads      — [meta, hostdep.fastq.gz]
        stats_json — [meta, host_stats.json]
        versions   — versions.yml

    minimap2 strategy:
        minimap2 → samtools view -f 4 → samtools fastq → pigz
        The -f 4 flag retains only UNMAPPED reads, i.e., non-host.

    hostile strategy:
        hostile clean --index <index> --fastq <reads> --out-dir hostile_out/
        Hostile renames outputs; we rename to canonical names.

    Stats JSON schema (both processes):
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

    minimap2 containers: minimap2 and samtools are in separate biocontainers.
    Because the pipeline command chains both tools in a single script block,
    we use the mulled container that includes both, or fall back to a conda
    env with both packages.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ---------------------------------------------------------------------------
// Process 1 — minimap2 + samtools host depletion (default)
// ---------------------------------------------------------------------------

process HOST_DEPLETE_MINIMAP2 {

    label 'process_high_memory'

    tag "${meta.id}"

    // The samtools biocontainer bundles samtools only; minimap2 is separate.
    // Use a mulled container that provides both tools.
    // mulled-v2 image for minimap2=2.28 + samtools=1.21:
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a969b254c884786eea2ca247ced:3161f532a5ea6f1ade5f7b9af6e853a844a2d2a3-0'
    conda 'bioconda::minimap2=2.28 bioconda::samtools=1.21 conda-forge::pigz'

    publishDir (
        path: { "${params.outdir}/${meta.id}/reads/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(reads)
    path  host_index   // pre-built .mmi file from INDEX_HOST process

    output:
    tuple val(meta), path("${meta.id}.hostdep.fastq.gz"),    emit: reads
    tuple val(meta), path("${meta.id}.host_stats.json"),     emit: stats_json
    tuple val(meta), path("${meta.id}.NO_VIRAL_READS_LIKELY"),
          optional: true,                                    emit: sentinel
    path "versions.yml",                                     emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Map reads to host, extract unmapped (non-host) reads
    # ----------------------------------------------------------------
    minimap2 \\
        -ax map-ont \\
        -t ${task.cpus} \\
        ${host_index} \\
        ${reads} \\
    | samtools view \\
        -@ ${task.cpus} \\
        -b \\
        -f 4 \\
        - \\
    | samtools fastq \\
        -@ ${task.cpus} \\
        - \\
    | pigz -p ${task.cpus} \\
    > ${meta.id}.hostdep.fastq.gz

    # ----------------------------------------------------------------
    # Compute statistics
    # ----------------------------------------------------------------
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
    import sys
    print(
        f"WARNING: {host_frac*100:.1f}% of reads mapped to the host reference for sample ${meta.id}. "
        "Very few viral reads remain. Pipeline continues — check results carefully.",
        file=sys.stderr
    )
    open("${meta.id}.NO_VIRAL_READS_LIKELY", "w").close()
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //')
    END_VERSIONS
    """

    stub:
    """
    echo -n "" | gzip > ${meta.id}.hostdep.fastq.gz
    echo '{"sample_id":"${meta.id}","total_reads":0,"host_reads_removed":0,"kept_reads":0,"host_fraction":0.0}' \\
        > ${meta.id}.host_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "2.28"
        samtools: "1.21"
        pigz: "2.8"
    END_VERSIONS
    """
}

// ---------------------------------------------------------------------------
// Process 2 — hostile host depletion (opt-in: params.use_hostile = true)
// ---------------------------------------------------------------------------

process HOST_DEPLETE_HOSTILE {

    label 'process_high_memory'

    tag "${meta.id}"

    container 'quay.io/biocontainers/hostile:1.1.0--pyhdfd78af_0'
    conda 'bioconda::hostile=1.1.0'

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
