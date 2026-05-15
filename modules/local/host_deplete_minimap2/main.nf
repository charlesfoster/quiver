/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HOST_DEPLETE_MINIMAP2 — Remove host (human) reads by mapping to GRCh38 no-alt.

    Default host depletion strategy: minimap2 + samtools.

    minimap2 strategy:
        minimap2 → samtools view -f 4 → samtools fastq → gzip
        The -f 4 flag retains only UNMAPPED reads, i.e., non-host.

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

    Container: mulled minimap2 2.28 + samtools 1.21
    Label: process_high_memory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HOST_DEPLETE_MINIMAP2 {

    label 'process_high_memory'

    tag "${meta.id}"

    // The samtools biocontainer bundles samtools only; minimap2 is separate.
    // Use a mulled container that provides both tools.
    // mulled-v2 image for minimap2=2.28 + samtools=1.21:
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"
    conda "${moduleDir}/environment.yml"

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
    | gzip -c \\
    > ${meta.id}.hostdep.fastq.gz

    # ----------------------------------------------------------------
    # Compute statistics (pure shell — mulled container has no python3)
    # ----------------------------------------------------------------
    total=\$(zcat ${reads} | awk 'NR%4==1' | wc -l | tr -d ' ')
    kept=\$(zcat  ${meta.id}.hostdep.fastq.gz | awk 'NR%4==1' | wc -l | tr -d ' ')
    host=\$(( total - kept ))
    # host_fraction rounded to 4 dp via awk
    host_frac=\$(awk -v t="\${total}" -v h="\${host}" \\
        'BEGIN { if (t > 0) printf "%.4f", h/t; else print "0.0000" }')

    # Emit JSON stats
    cat > ${meta.id}.host_stats.json <<JSONEOF
{
  "sample_id": "${meta.id}",
  "total_reads": \${total},
  "host_reads_removed": \${host},
  "kept_reads": \${kept},
  "host_fraction": \${host_frac}
}
JSONEOF

    # Warn and flag if >99.5% reads mapped to host
    if awk -v f="\${host_frac}" 'BEGIN { exit (f + 0 >= 0.995) ? 0 : 1 }'; then
        echo "WARNING: \${host_frac} of reads mapped to host for sample ${meta.id}. Very few viral reads remain." >&2
        touch ${meta.id}.NO_VIRAL_READS_LIKELY
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
        gzip: \$(gzip --version 2>&1 | head -1 | sed 's/gzip //')
    END_VERSIONS
    """

    stub:
    """
    echo -n "" | gzip -c > ${meta.id}.hostdep.fastq.gz
    echo '{"sample_id":"${meta.id}","total_reads":0,"host_reads_removed":0,"kept_reads":0,"host_fraction":0.0}' \\
        > ${meta.id}.host_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "2.28"
        samtools: "1.21"
        gzip: "1.12"
    END_VERSIONS
    """
}
