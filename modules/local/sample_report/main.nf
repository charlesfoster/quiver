/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SAMPLE_REPORT — Per-sample HTML and JSON summary report.

    Purpose:
        Implements Step 5.21 of the data flow specification.  Wraps
        `bin/render_sample_report.py` to produce a self-contained
        `${SID}_summary.html` and `${SID}_summary.json` from all per-sample
        per-genotype outputs collected by the calling workflow.

    Self-contained HTML guarantee:
        The Jinja2 template (`assets/templates/sample_report.html.j2`) uses only
        inline CSS and embeds NanoPlot PNGs as base64 data-URIs.  The resulting
        HTML file opens in a browser without any network requests.

    Optional inputs:
        Every input beyond `genotype_summary` is optional.  The Python script
        handles missing files gracefully — it produces a partial report rather
        than failing.  This means:
          - Samples with LOW_COVERAGE get a report without coverage/variant sections.
          - Samples where DEVIDER failed get a report with a "failed" note.
          - Samples without NanoPlot dirs get a report without embedded images.

    Sentinel flags:
        Sentinel files (*.LOW_COVERAGE, *.NO_HCV_DETECTED, etc.) emitted by
        upstream processes are collected into `flag_files` and passed to the
        script via --flags.  The script reads each file's extension to determine
        the flag name and renders colour-coded badges.

    Container:
        quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0 — reuses the MultiQC
        container, which ships Jinja2 as a dependency and includes procps/ps
        required by Nextflow process monitoring under Docker/WSL.

    Label: process_low (1 CPU, 2 GB, 5 min — Step 5.21 resource spec).

    Inputs:
        meta            — val map: at minimum { id: <sample_id> }
        genotype_summary — path to *_genotype_summary.json (required)
        nanoq_raw       — path to raw *.nanoq.json (optional)
        nanoq_filtered  — path to post-filter *.nanoq.json (optional)
        host_stats      — path to *.host_stats.json (optional)
        mosdepth_files  — list of *.mosdepth.summary.txt from all branches (optional)
        variant_tsvs    — list of *_variants.tsv from all branches (optional)
        flagstat_files  — list of *_round2.flagstat from all branches (optional)
        haplotype_reports — list of *_haplotype_report.json from all branches (optional)
        nanoplot_dirs   — list of NanoPlot output directories (optional)
        flag_files      — list of sentinel flag files (optional)

    Outputs:
        html     — [meta, "${meta.id}_summary.html"]
        json     — [meta, "${meta.id}_summary.json"]
        versions — versions.yml

    PublishDir: ${params.outdir}/${meta.id}/reports/

    See also:
        docs/data_flow.md              Step 5.21
        docs/implementation_prompts.md Prompt 22
        bin/render_sample_report.py    Report generation logic
        assets/templates/              Jinja2 HTML template
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SAMPLE_REPORT {

    label 'process_low'

    tag "${meta.id}"

    // MultiQC container — has Jinja2 pre-installed (MultiQC depends on it)
    // and includes procps/ps required by Nextflow process monitoring.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/reports/" },
        mode: 'copy'
    )

    input:
    tuple val(meta),
          path(genotype_summary),   // required: *_genotype_summary.json
          path(nanoq_raw),          // optional: raw nanoq JSON
          path(nanoq_filtered),     // optional: post-filter nanoq JSON
          path(host_stats),         // optional: host_stats.json
          path(mosdepth_files),     // optional: list of *.mosdepth.summary.txt
          path(mosdepth_beds),      // optional: list of *.regions.bed.gz (per-window coverage)
          path(variant_tsvs),       // optional: list of *_variants.tsv
          path(flagstat_files),     // optional: list of *_round2.flagstat
          path(haplotype_reports),   // optional: list of *_haplotype_report.json
          path(nanoplot_dirs),      // optional: list of NanoPlot directories
          path(flag_files)          // optional: list of *.FLAGNAME sentinel files

    output:
    tuple val(meta),
          path("${meta.id}_summary.html"),
          emit: html
    tuple val(meta),
          path("${meta.id}_summary.json"),
          emit: json
    path "versions.yml", emit: versions

    script:
    // Build CLI argument fragments.
    //
    // Each optional block uses a ternary: if the variable is a non-empty
    // path/list, emit the flag; otherwise emit an empty string.
    //
    // Nextflow stages optional path inputs as empty strings when not supplied;
    // we check for actual files with -f / -d guards.

    def nanoq_raw_arg      = nanoq_raw      ? "--nanoq-raw ${nanoq_raw}"           : ""
    def nanoq_filt_arg     = nanoq_filtered ? "--nanoq-filtered ${nanoq_filtered}"  : ""
    def host_stats_arg     = host_stats     ? "--host-stats ${host_stats}"          : ""

    def _ri_b64 = groovy.json.JsonOutput.toJson([
        pipeline_version: (workflow.manifest.version ?: 'QuIVER'),
        nextflow_version: workflow.nextflow.version.toString(),
        params          : params.findAll { true },
    ]).bytes.encodeBase64().toString()

    // For list inputs (mosdepth_files, etc.) Nextflow stages them as space-separated
    // filenames in the work directory.  We pass them as-is; the Python script
    // iterates them with nargs="*".
    """
    python3 -c "import base64,json; open('run_info.json','w').write(json.dumps(json.loads(base64.b64decode('${_ri_b64}')),indent=2))"

    # ----------------------------------------------------------------
    # Build the optional argument lists.
    #
    # We expand each optional path collection into its constituent
    # filenames using shell globs so the Python argparse nargs="*"
    # receives individual file paths.
    #
    # Sentinel files are collected by matching all flag-named extensions.
    # ----------------------------------------------------------------

    mosdepth_args=""
    if compgen -G "*.mosdepth.summary.txt" > /dev/null 2>&1; then
        mosdepth_args="--mosdepth-summaries \$(ls *.mosdepth.summary.txt | tr '\\n' ' ')"
    fi

    mosdepth_bed_args=""
    if compgen -G "*.regions.bed.gz" > /dev/null 2>&1; then
        mosdepth_bed_args="--mosdepth-beds \$(ls *.regions.bed.gz | tr '\\n' ' ')"
    fi

    variant_args=""
    if compgen -G "*_variants.tsv" > /dev/null 2>&1; then
        variant_args="--variant-tsvs \$(ls *_variants.tsv | tr '\\n' ' ')"
    fi

    flagstat_args=""
    if compgen -G "*.flagstat" > /dev/null 2>&1; then
        flagstat_args="--flagstats \$(ls *.flagstat | tr '\\n' ' ')"
    fi

    haplotype_args=""
    if compgen -G "*_haplotype_report.json" > /dev/null 2>&1; then
        haplotype_args="--haplotype-reports \$(ls *_haplotype_report.json | tr '\\n' ' ')"
    fi

    nanoplot_args=""
    if ls -d *_nanoplot/ > /dev/null 2>&1; then
        nanoplot_args="--nanoplot-dirs \$(ls -d *_nanoplot/ | tr '\\n' ' ')"
    fi

    flag_args=""
    # Collect all recognised sentinel files (any file matching *.UPPER_CASE extension)
    sentinel_files=\$(find . -maxdepth 1 -name "*.NO_HCV_DETECTED" \\
                                          -o -name "*.LOW_COVERAGE" \\
                                          -o -name "*.LOW_MAPPING_RATE" \\
                                          -o -name "*.NO_VIRAL_READS_LIKELY" \\
                                          -o -name "*.ALL_READS_FILTERED" \\
                                          -o -name "*.EMPTY_INPUT" \\
                                          -o -name "*.LOW_COVERAGE_CONSENSUS" \\
                                          -o -name "*.DEVIDER_FAILED" \\
                                          2>/dev/null | tr '\\n' ' ')
    if [ -n "\${sentinel_files}" ]; then
        flag_args="--flags \${sentinel_files}"
    fi

    # ----------------------------------------------------------------
    # Run the report generator.
    # ----------------------------------------------------------------
    python3 ${projectDir}/bin/render_sample_report.py \\
        --sample-id ${meta.id} \\
        --genotype-summary ${genotype_summary} \\
        ${nanoq_raw_arg} \\
        ${nanoq_filt_arg} \\
        ${host_stats_arg} \\
        --run-info run_info.json \\
        --template ${projectDir}/assets/templates/sample_report.html.j2 \\
        --pipeline-version "${workflow.manifest.version ?: 'QuIVER'}" \\
        --output-html ${meta.id}_summary.html \\
        --output-json ${meta.id}_summary.json \\
        \${mosdepth_args} \\
        \${mosdepth_bed_args} \\
        \${variant_args} \\
        \${flagstat_args} \\
        \${haplotype_args} \\
        \${nanoplot_args} \\
        \${flag_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        jinja2: \$(python3 -c 'import jinja2; print(jinja2.__version__)')
    END_VERSIONS
    """

    stub:
    """
    cat <<-HTMLEOF > ${meta.id}_summary.html
    <!DOCTYPE html><html><body><h1>Stub report for ${meta.id}</h1></body></html>
    HTMLEOF

    cat <<-JSONEOF > ${meta.id}_summary.json
    {
      "sample_id":        "${meta.id}",
      "overall_status":   "PASS",
      "is_mixed":         false,
      "primary_genotype": null,
      "flags":            [],
      "_stub":            true
    }
    JSONEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.11"
        jinja2: "3.1"
    END_VERSIONS
    """
}
