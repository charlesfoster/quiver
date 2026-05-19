/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DOWNLOAD_HOST_REFERENCE — Download GRCh38 no-alt and build minimap2 index.

    Triggered automatically when:
        --host_reference is not supplied  AND
        --skip_host_depletion is false    AND
        the cached .mmi does not exist in params.host_genome_cache  OR
        --force_host_genome_download is set.

    The finished .mmi (and compressed FASTA) are published to
    params.host_genome_cache so they persist across pipeline runs.
    On subsequent runs the workflow reads directly from that path and
    this process is skipped entirely.

    Download source:
        GRCh38 no-alt analysis set (NCBI, hg38)
        https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/
            GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/
            GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
        ~1 GB compressed; minimap2 indexing requires ~8 GB RAM and ~15 min.

    Outputs (in work dir, then published to host_genome_cache):
        GRCh38_no_alt.fna.gz  — compressed FASTA (kept so re-indexing is fast)
        GRCh38_no_alt.mmi     — minimap2 index (-x map-ont)

    Container: minimap2 biocontainer (has minimap2 + curl via conda-forge base)
    Conda:     bioconda::minimap2  +  conda-forge::curl
    Label:     process_high  (indexing needs 8+ GB RAM, ~15 min, 8 threads)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DOWNLOAD_HOST_REFERENCE {

    label 'process_high'

    tag "GRCh38_no_alt"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0' :
        'quay.io/biocontainers/minimap2:2.28--he4a0461_0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path:      params.host_genome_cache,
        mode:      'copy',
        overwrite: true
    )

    output:
    path "GRCh38_no_alt.mmi",    emit: mmi
    path "GRCh38_no_alt.fna.gz", emit: fasta
    path "versions.yml",          emit: versions

    script:
    def grch38_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    """
    # ----------------------------------------------------------------
    # Download GRCh38 no-alt analysis set
    # ----------------------------------------------------------------
    echo "[download_host_reference] Downloading GRCh38 no-alt (~1 GB)..." >&2
    curl --retry 3 --retry-delay 5 -L -o GRCh38_no_alt.fna.gz \\
        "${grch38_url}"

    # Validate gzip integrity
    gzip -t GRCh38_no_alt.fna.gz || {
        echo "ERROR: Downloaded file is corrupt." >&2
        exit 1
    }

    # ----------------------------------------------------------------
    # Build minimap2 index (~8 GB RAM, ~15 min)
    # ----------------------------------------------------------------
    echo "[download_host_reference] Building minimap2 index..." >&2
    minimap2 \\
        -x map-ont \\
        -t ${task.cpus} \\
        -d GRCh38_no_alt.mmi \\
        GRCh38_no_alt.fna.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        curl: \$(curl --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch GRCh38_no_alt.fna.gz
    touch GRCh38_no_alt.mmi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "2.28"
        curl: "8.0.0"
    END_VERSIONS
    """
}
