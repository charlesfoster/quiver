/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DOWNLOAD_NOHUMAN_DB — Download the nohuman Kraken2 human-genome database.

    Triggered automatically when:
        --use_nohuman is active (or is the default)  AND
        the database directory does not exist in params.nohuman_db  OR
        --force_nohuman_db_download is set.

    The database (~4 GB) is published to params.nohuman_db so it persists
    across pipeline runs. On subsequent runs the workflow detects the cache
    and this process is skipped entirely.

    Download source:
        nohuman fetches a curated human HPRC Kraken2 database via the --download
        flag. Use --db-version to pin a specific release (default: latest HPRC.r2).
        Run `nohuman --list-db-versions` inside the container for available versions.

    Outputs (published to the parent directory of params.nohuman_db):
        <db_name>/hash.k2d  — Kraken2 k-mer database
        <db_name>/opts.k2d  — Kraken2 options
        <db_name>/taxo.k2d  — Kraken2 taxonomy

    Container: nohuman biocontainer (verify tag at quay.io/biocontainers/nohuman)
    Conda:     bioconda::nohuman
    Label:     process_high  (~4 GB download + decompression)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DOWNLOAD_NOHUMAN_DB {

    label 'process_high'

    tag "nohuman_db"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nohuman:0.5.1--hfa8f182_1' :
        'quay.io/biocontainers/nohuman:0.5.1--hfa8f182_1' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: file(params.nohuman_db).parent,
        mode: 'copy',
        overwrite: true
    )

    output:
    path "${file(params.nohuman_db).name}", emit: db
    path "versions.yml",                    emit: versions

    script:
    def db_name = file(params.nohuman_db).name
    """
    echo "[download_nohuman_db] Downloading nohuman HPRC Kraken2 database (~4 GB)..." >&2
    nohuman --download --db ${db_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nohuman: \$(nohuman --version 2>&1 | sed 's/nohuman //')
    END_VERSIONS
    """

    stub:
    def db_name = file(params.nohuman_db).name
    """
    mkdir -p ${db_name}
    touch ${db_name}/hash.k2d
    touch ${db_name}/opts.k2d
    touch ${db_name}/taxo.k2d

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nohuman: "0.5.1"
    END_VERSIONS
    """
}
