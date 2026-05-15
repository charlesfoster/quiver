/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MULTIQC — Cross-sample QC aggregation using MultiQC.

    Purpose:
        Implements Step 5.22 of the data flow specification.  Collects all per-sample
        QC files from all pipeline stages (NanoPlot, mosdepth, samtools flagstat,
        bcftools stats) and produces a single comprehensive MultiQC report.

    Parsers used (MultiQC 1.25.1):
        nanoplot      NanoStats.txt (tsv_stats output from NANOPLOT process)
        mosdepth      *.mosdepth.summary.txt (from MOSDEPTH process)
        samtools      *.flagstat (from MINIMAP2_ROUND2 process)
        bcftools      *.bcftools_stats.txt (from bcftools stats post-processing)

    MultiQC input strategy:
        All QC files from all samples are collected with `.collect()` in the calling
        workflow and passed as a flat channel.  MultiQC's recursive search then finds
        each tool's output by filename pattern as configured in
        `assets/multiqc_config.yml`.

    Configuration file:
        `assets/multiqc_config.yml` sets:
          - Report title and subtitle
          - Module list (nanoplot, mosdepth, samtools, bcftools)
          - Table column visibility and names
          - Search pattern overrides

    Container:
        quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0
        (CLAUDE.md Section 4, pinned to MultiQC 1.25.1)

    Label: process_low (2 CPU, 4 GB, 10 min — Step 5.22 resource spec).

    Outputs:
        report  — multiqc_report/multiqc_report.html (self-contained HTML)
        data    — multiqc_report/multiqc_data/ (JSON data for downstream tools)
        versions — versions.yml

    PublishDir: ${params.outdir}/reports/

    See also:
        docs/data_flow.md              Step 5.22
        docs/implementation_prompts.md Prompt 23
        assets/multiqc_config.yml      MultiQC configuration
        CLAUDE.md                      D15 — QC aggregation decision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MULTIQC {

    label 'process_low'

    tag "multiqc"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: "${params.outdir}/reports/",
        mode: 'copy'
    )

    input:
    // All QC files from all samples passed as a flat collection.
    // The calling workflow uses:
    //   ch_multiqc_files.collect()
    // to gather NanoPlot NanoStats.txt, mosdepth summary, flagstat, and
    // bcftools stats files from all samples before feeding this process.
    path(multiqc_files, stageAs: "multiqc_input/*")

    // MultiQC configuration file (assets/multiqc_config.yml)
    path multiqc_config

    output:
    path "multiqc_report/multiqc_report.html", emit: report
    path "multiqc_report/multiqc_data/",       emit: data
    path "versions.yml",                       emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Run MultiQC across all staged QC files.
    #
    # --config  : use the pipeline's multiqc_config.yml for title,
    #             module selection, table customisation.
    # -f        : force overwrite of any existing multiqc_report/.
    # --outdir  : explicit output directory.
    # --title   : CLI title overrides config title in some MultiQC
    #             versions; we pass both for robustness.
    # -q        : quiet mode (suppress per-file log lines in NF log).
    #
    # MultiQC's recursive file search is started from the staged
    # input directory (multiqc_input/).  It finds tool outputs by
    # filename pattern as defined in multiqc_config.yml and in
    # MultiQC's built-in search patterns.
    # ----------------------------------------------------------------
    multiqc \\
        --config ${multiqc_config} \\
        --title "hcv-quasi pipeline QC report" \\
        --force \\
        --outdir multiqc_report \\
        -q \\
        multiqc_input/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | sed 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p multiqc_report/multiqc_data
    cat <<-HTMLEOF > multiqc_report/multiqc_report.html
    <!DOCTYPE html><html><body>
    <h1>MultiQC stub report</h1>
    <p>hcv-quasi pipeline QC report (stub run)</p>
    </body></html>
    HTMLEOF

    echo '{"report_general_stats_data": []}' > multiqc_report/multiqc_data/multiqc_general_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: "1.25.1"
    END_VERSIONS
    """
}
