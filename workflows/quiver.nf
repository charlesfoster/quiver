/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QuIVER — Top-level workflow.

    Wiring overview (see CLAUDE.md Section 2 for the ASCII diagram and
    docs/data_flow.md for the per-step rationale):

        INPUT_CHECK
            -> RAW_QC (NanoPlot + nanoq) ─── parallel branch
            -> CHOPPER  (read filter)
            -> NANOQ_FILT (post-filter QC)
            -> HOST_DEPLETE_{MINIMAP2|HOSTILE}
            -> NANOQ_POSTHOST
            -> INDEX_PANEL (singleton)
            -> MINIMAP2_ROUND1
            -> GENOTYPE_CLASSIFY
            -> GENOTYPE_BRANCH ─── fan-out one item per [sample × genotype]
                |
                |   (per-branch, branch_meta = meta + [genotype: gt])
                v
                BUILD_CONSENSUS
                    -> MINIMAP2_ROUND2 (full-depth)
                    -> MOSDEPTH
                    -> PREP_LOFREQ_INPUT
                        -> LOFREQ_CALL -> VARIANT_FILTER
                        -> (optional) CLAIR3 + concordance
                    -> PREP_DEVIDER_INPUT
                        -> DEVIDER (gated by LOW_COVERAGE)
                        -> FORMAT_HAPLOTYPES
                |
                v
            groupTuple(by: 0) on meta.id  ── per-sample collapse
                -> SAMPLE_REPORT
            collect QC across all samples + branches
                -> MULTIQC
                -> RENDER_RUN_SUMMARY

    Failure-mode routing (architecture_reasoning.md §12):
        EMPTY_INPUT             — INPUT_CHECK errors before workflow proceeds.
        ALL_READS_FILTERED      — CHOPPER sentinel; sample short-circuits to flag.
        NO_VIRAL_READS_LIKELY   — HOST_DEPLETE sentinel; warning, pipeline continues.
        NO_HCV_DETECTED         — MINIMAP2_ROUND1 / GENOTYPE_BRANCH; sample routed
                                  to no_hcv channel, no downstream processing.
        LOW_COVERAGE            — MOSDEPTH sentinel; skip DEVIDER for that branch
                                  but keep LoFreq.  Flag passed to SAMPLE_REPORT.
        LOW_COVERAGE_CONSENSUS  — APPLY_CONSENSUS sentinel; downstream continues,
                                  flag passed to SAMPLE_REPORT.
        devider.failed          — DEVIDER process; FORMAT_HAPLOTYPES degrades
                                  gracefully (empty FASTA + fallback JSON).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ---------------------------------------------------------------------------
// Module / subworkflow includes
// ---------------------------------------------------------------------------

include { INPUT_CHECK              } from '../modules/local/input_check/main'
include { RAW_QC                   } from '../subworkflows/local/raw_qc'
include { CHOPPER                  } from '../modules/local/chopper/main'
include { NANOQ as NANOQ_FILT      } from '../modules/local/nanoq/main'
include { NANOQ as NANOQ_POSTHOST  } from '../modules/local/nanoq/main'
include { DOWNLOAD_NOHUMAN_DB      } from '../modules/local/download_nohuman_db/main'
include { HOST_DEPLETE_NOHUMAN     } from '../modules/local/host_deplete_nohuman/main'
include { DOWNLOAD_HOST_REFERENCE  } from '../modules/local/download_host_reference/main'
include { INDEX_HOST               } from '../modules/local/index_host/main'
include { HOST_DEPLETE_MINIMAP2    } from '../modules/local/host_deplete_minimap2/main'
include { HOST_DEPLETE_HOSTILE     } from '../modules/local/host_deplete_hostile/main'
include { INDEX_PANEL              } from '../modules/local/index_panel/main'
include { MINIMAP2_ROUND1          } from '../modules/local/minimap2_round1/main'
include { GENOTYPE_CLASSIFY        } from '../modules/local/genotype_classify/main'
include { GENOTYPE_BRANCH          } from '../subworkflows/local/genotype_branch'
include { BUILD_CONSENSUS          } from '../subworkflows/local/build_consensus'
include { MINIMAP2_ROUND2          } from '../modules/local/minimap2_round2/main'
include { MOSDEPTH                 } from '../modules/local/mosdepth/main'
include { PREP_LOFREQ_INPUT        } from '../subworkflows/local/prep_lofreq_input'
include { LOFREQ_CALL              } from '../modules/local/lofreq_call/main'
include { VARIANT_FILTER           } from '../modules/local/variant_filter/main'
include { CLAIR3                   } from '../modules/local/clair3/main'
include { PREP_DEVIDER_INPUT       } from '../subworkflows/local/prep_devider_input'
include { FILTER_VCF_FOR_DEVIDER   } from '../modules/local/filter_vcf_for_devider/main'
include { DEVIDER                  } from '../modules/local/devider/main'
include { FORMAT_HAPLOTYPES        } from '../modules/local/format_haplotypes/main'
include { MAKE_ROUND2_MASK         } from '../modules/local/make_round2_mask/main'
include { BUILD_ROUND2_CONSENSUS   } from '../modules/local/build_round2_consensus/main'
include { SAMPLE_REPORT            } from '../modules/local/sample_report/main'
include { MULTIQC                  } from '../modules/local/multiqc/main'


// ---------------------------------------------------------------------------
// Local helper processes — scoped to this workflow
// ---------------------------------------------------------------------------

/*
    DUMP_SOFTWARE_VERSIONS — collect all per-process versions.yml fragments
    into a single software_versions.yml under the run pipeline_info directory.

    Each module emits a YAML fragment of the form:
        "<task.process>":
            tool: version
    We concatenate them verbatim (no deduplication needed; Nextflow already
    caches and resumes by content hash, so re-runs of the same process emit
    the same fragment).
*/
process DUMP_SOFTWARE_VERSIONS {

    label 'process_low'

    tag 'software_versions'

    container 'python:3.11'
    conda 'conda-forge::python=3.11'

    publishDir (
        path: "${params.outdir}/pipeline_info/",
        mode: 'copy'
    )

    input:
    // Stage each input file under a unique name to avoid the
    // "input file name collision" error — every module emits
    // a file literally called `versions.yml`.
    path versions_yamls, stageAs: "versions_*/versions.yml"

    output:
    path "software_versions.yml"

    script:
    """
    {
        echo "# QuIVER software versions"
        echo "# pipeline_version: ${workflow.manifest.version ?: 'unknown'}"
        echo "# nextflow_version: ${nextflow.version}"
        echo "# run_date: \$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo
        for f in versions_*/versions.yml; do
            cat "\$f"
            echo
        done
    } > software_versions.yml
    """

    stub:
    """
    echo "# stub software versions" > software_versions.yml
    """
}

/*
    RENDER_RUN_SUMMARY — run-level HTML + JSON aggregation of per-sample summaries.

    Wraps bin/render_run_summary.py.  Receives the collected list of per-sample
    *_summary.json files; emits run_summary.html and run_summary.json under
    ${params.outdir}/reports/.
*/
process RENDER_RUN_SUMMARY {

    label 'process_low'

    tag 'run_summary'

    container 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'
    conda 'conda-forge::python=3.11 conda-forge::jinja2'

    publishDir (
        path: "${params.outdir}/reports/",
        mode: 'copy'
    )

    input:
    path sample_jsons

    output:
    path "run_summary.html"
    path "run_summary.json"

    script:
    // Read samplesheet to capture the original input order of sample IDs.
    // This lets the run summary table match the user's samplesheet ordering.
    def _sample_order = []
    try {
        def _ss = file(params.input)
        if (_ss.exists()) {
            _ss.readLines().drop(1).each { line ->
                def sid = line.split(',')[0].trim().replaceAll('"', '').replaceAll("'", '')
                if (sid) _sample_order << sid
            }
        }
    } catch (Exception _e) { /* ignore — table falls back to status sort */ }

    def _ri_b64 = groovy.json.JsonOutput.toJson([
        pipeline_version: (workflow.manifest.version ?: 'QuIVER'),
        nextflow_version: workflow.nextflow.version.toString(),
        params          : params.findAll { true },
        sample_order    : _sample_order,
    ]).bytes.encodeBase64().toString()
    """
    python3 -c "import base64,json; open('run_info.json','w').write(json.dumps(json.loads(base64.b64decode('${_ri_b64}')),indent=2))"

    python3 ${projectDir}/bin/render_run_summary.py \\
        --sample-jsons ${sample_jsons} \\
        --run-info run_info.json \\
        --output-html run_summary.html \\
        --output-json run_summary.json \\
        --pipeline-version "${workflow.manifest.version ?: 'QuIVER'}"
    """

    stub:
    """
    echo "<html><body>Stub run summary</body></html>" > run_summary.html
    echo '{"run":"stub"}' > run_summary.json
    """
}

/*
    EMIT_FAILURE_FLAG — produce a published flag file for samples that hit a
    pipeline-terminating sentinel (NO_HCV_DETECTED, ALL_READS_FILTERED, etc.).

    The flag file lives under ${params.outdir}/${meta.id}/ so analysts can see
    at a glance why a sample was short-circuited.
*/
process EMIT_FAILURE_FLAG {

    label 'process_low'

    tag "${meta.id}"

    container 'python:3.11'
    conda 'conda-forge::python=3.11'

    publishDir (
        path: { "${params.outdir}/${meta.id}/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), val(reason)

    output:
    tuple val(meta), path("${meta.id}.PIPELINE_FLAG.txt"), emit: flag

    script:
    """
    {
        echo "Pipeline status flag: ${reason}"
        echo "Sample: ${meta.id}"
        echo "Timestamp: \$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo
        echo "This sample was short-circuited by the QuIVER pipeline."
        echo "See ${params.outdir}/${meta.id}/ for the upstream sentinel file."
    } > ${meta.id}.PIPELINE_FLAG.txt
    """

    stub:
    """
    echo "Sample: ${meta.id}" > ${meta.id}.PIPELINE_FLAG.txt
    echo "Reason: ${reason}" >> ${meta.id}.PIPELINE_FLAG.txt
    """
}


// ===========================================================================
//                              MAIN WORKFLOW
// ===========================================================================

workflow QUIVER {

    // ---------------------------------------------------------------------
    // Master version channel — every module appends its versions.yml here.
    // Aggregated at the end via DUMP_SOFTWARE_VERSIONS.
    // ---------------------------------------------------------------------
    ch_versions = Channel.empty()

    // =====================================================================
    // Parameter validation
    // =====================================================================
    if (!params.input) {
        error "Required parameter --input (samplesheet CSV) is not set."
    }
    if (params.min_bq > params.min_alt_bq) {
        error "Parameter error: --min_bq (${params.min_bq}) must be <= --min_alt_bq (${params.min_alt_bq}). " +
              "LoFreq requires that the general base-quality floor (Q) does not exceed " +
              "the alternate-allele base-quality floor (q)."
    }
    if (params.min_length >= params.max_length) {
        error "Parameter error: --min_length (${params.min_length}) must be < --max_length (${params.max_length})."
    }
    if (params.min_qual < 0) {
        error "Parameter error: --min_qual (${params.min_qual}) must be >= 0."
    }
    if (params.lofreq_sig <= 0 || params.lofreq_sig > 1) {
        error "Parameter error: --lofreq_sig (${params.lofreq_sig}) must be in the range (0, 1]."
    }
    if ((params.lofreq_pp_threads as int) < 1) {
        error "Parameter error: --lofreq_pp_threads (${params.lofreq_pp_threads}) must be >= 1."
    }
    if (params.devider_min_abund <= 0 || params.devider_min_abund >= 1) {
        error "Parameter error: --devider_min_abund (${params.devider_min_abund}) must be in the range (0, 1)."
    }
    if (params.devider_min_af <= 0 || params.devider_min_af > 1) {
        error "Parameter error: --devider_min_af (${params.devider_min_af}) must be in the range (0, 1]."
    }
    if (params.devider_min_af < params.min_report_af) {
        error "Parameter error: --devider_min_af (${params.devider_min_af}) must be >= --min_report_af (${params.min_report_af}). " +
              "The DEVIDER VCF filter cannot be looser than the analysis VCF filter."
    }
    if (params.min_secondary_fraction < 0 || params.min_secondary_fraction >= 1) {
        error "Parameter error: --min_secondary_fraction (${params.min_secondary_fraction}) must be in the range [0, 1)."
    }


    // =====================================================================
    // Step 5.1 — INPUT_CHECK
    // =====================================================================
    ch_samplesheet = Channel.fromPath(params.input, checkIfExists: true)
    INPUT_CHECK(ch_samplesheet)

    // The python validator writes a JSON array of {sample_id, fastq, metadata}
    // records.  Parse it back into a per-sample channel of [meta, fastq].
    ch_reads = INPUT_CHECK.out.json
        .flatMap { json_path ->
            // json_path is a java.nio Path; JsonSlurper.parse accepts an
            // InputStream or a Reader.  Use newReader() (a Path method) to
            // hand it a Reader without needing a File cast.
            def records = new groovy.json.JsonSlurper().parse(json_path.newReader())
            records.collect { rec ->
                def meta = [
                    id:       rec.id as String,
                    metadata: rec.metadata ?: [:],
                ]
                tuple(meta, file(rec.fastq as String, checkIfExists: true))
            }
        }


    // =====================================================================
    // Step 5.2 — RAW_QC (NanoPlot + nanoq, parallel)
    // =====================================================================
    RAW_QC(ch_reads)
    ch_versions = ch_versions.mix(RAW_QC.out.versions)


    // =====================================================================
    // Step 5.3 — CHOPPER read filtering
    // =====================================================================
    CHOPPER(ch_reads)
    ch_versions = ch_versions.mix(CHOPPER.out.versions)

    // Samples whose filter sentinel fired short-circuit here.  We gate on the
    // absence of an ALL_READS_FILTERED file by joining the reads channel with
    // the sentinel channel (left-outer join semantics via `.join(..., remainder: true)`)
    // and routing accordingly.
    ch_filtered_routed = CHOPPER.out.reads
        .join(CHOPPER.out.sentinel, by: 0, remainder: true)
        .branch { meta, reads, sentinel ->
            // sentinel is null when the optional emit did not produce a file
            filtered:    sentinel == null
            all_filtered: true
        }

    ch_post_filter_reads = ch_filtered_routed.filtered.map { meta, reads, sentinel -> tuple(meta, reads) }

    // Emit a published PIPELINE_FLAG for ALL_READS_FILTERED samples.
    ch_all_filtered_failed = ch_filtered_routed.all_filtered.map { meta, reads, sentinel ->
        tuple(meta, 'ALL_READS_FILTERED')
    }


    // =====================================================================
    // Step 5.3b — post-filter nanoq (re-uses the NANOQ process, distinct name)
    // =====================================================================
    NANOQ_FILT(ch_post_filter_reads)
    ch_versions = ch_versions.mix(NANOQ_FILT.out.versions)


    // =====================================================================
    // Step 5.4 — HOST_DEPLETE
    //
    // Three methods — specify at most one flag; nohuman is the default:
    //   --use_nohuman   (default) Kraken2-based; DB auto-downloaded to params.nohuman_db
    //   --use_hostile             alignment-based; requires --host_reference (hostile index dir)
    //   --use_minimap2            alignment vs GRCh38; resolves reference:
    //       --host_reference <.mmi>  → use pre-built index directly
    //       --host_reference <FASTA> → build index via INDEX_HOST
    //       (none)                   → auto-download GRCh38 no-alt + build .mmi;
    //                                  cached in params.host_genome_cache for future runs
    // =====================================================================
    def _method_count = [params.use_nohuman, params.use_hostile, params.use_minimap2].count { it }
    if (_method_count > 1) {
        error "Specify at most one of --use_nohuman, --use_hostile, --use_minimap2."
    }
    def host_method = params.use_hostile  ? 'hostile'
                    : params.use_minimap2 ? 'minimap2'
                    : 'nohuman'

    if (params.skip_host_depletion) {
        log.warn "Host depletion skipped (--skip_host_depletion). Not recommended for clinical samples."
        ch_hostdep_reads = ch_post_filter_reads
        ch_host_stats    = Channel.empty()

    } else if (host_method == 'hostile') {
        if (!params.host_reference) {
            error "hostile (--use_hostile) requires --host_reference pointing to a hostile index directory."
        }
        def ch_hostile_idx = Channel.fromPath(params.host_reference, checkIfExists: true).first()
        HOST_DEPLETE_HOSTILE(ch_post_filter_reads, ch_hostile_idx)
        ch_versions      = ch_versions.mix(HOST_DEPLETE_HOSTILE.out.versions)
        ch_hostdep_reads = HOST_DEPLETE_HOSTILE.out.reads
        ch_host_stats    = HOST_DEPLETE_HOSTILE.out.stats_json

    } else if (host_method == 'nohuman') {
        def ch_nohuman_db
        def cached_db = file(params.nohuman_db)
        def db_is_cached = cached_db.isDirectory() && cached_db.list().size() > 0
        if (db_is_cached && !params.force_nohuman_db_download) {
            log.info "Using cached nohuman database: ${params.nohuman_db}"
            ch_nohuman_db = Channel.value(cached_db)
        } else {
            def reason = params.force_nohuman_db_download
                ? "Re-downloading nohuman database (--force_nohuman_db_download)"
                : "No nohuman database found in cache (${params.nohuman_db}); downloading (~4 GB)"
            log.warn reason
            DOWNLOAD_NOHUMAN_DB()
            ch_versions   = ch_versions.mix(DOWNLOAD_NOHUMAN_DB.out.versions)
            ch_nohuman_db = DOWNLOAD_NOHUMAN_DB.out.db.first()
        }
        HOST_DEPLETE_NOHUMAN(ch_post_filter_reads, ch_nohuman_db)
        ch_versions      = ch_versions.mix(HOST_DEPLETE_NOHUMAN.out.versions)
        ch_hostdep_reads = HOST_DEPLETE_NOHUMAN.out.reads
        ch_host_stats    = HOST_DEPLETE_NOHUMAN.out.stats_json

    } else {
        // minimap2 alignment-based depletion
        def ch_host_mmi
        if (params.host_reference) {
            def href = params.host_reference.toString()
            if (href.endsWith('.mmi')) {
                ch_host_mmi = Channel.fromPath(params.host_reference, checkIfExists: true).first()
            } else {
                def ch_host_fasta = Channel.fromPath(params.host_reference, checkIfExists: true)
                INDEX_HOST(ch_host_fasta)
                ch_versions = ch_versions.mix(INDEX_HOST.out.versions)
                ch_host_mmi = INDEX_HOST.out.index.first()
            }
        } else {
            def cached_mmi = file("${params.host_genome_cache}/GRCh38_no_alt.mmi")
            if (cached_mmi.exists() && !params.force_host_genome_download) {
                log.info "Using cached host reference: ${cached_mmi}"
                ch_host_mmi = Channel.value(cached_mmi)
            } else {
                def reason = cached_mmi.exists()
                    ? "Re-downloading host reference (--force_host_genome_download)"
                    : "No host reference found in cache (${params.host_genome_cache}); downloading GRCh38 no-alt (~1 GB, ~15 min)"
                log.warn reason
                DOWNLOAD_HOST_REFERENCE()
                ch_versions = ch_versions.mix(DOWNLOAD_HOST_REFERENCE.out.versions)
                ch_host_mmi = DOWNLOAD_HOST_REFERENCE.out.mmi.first()
            }
        }
        HOST_DEPLETE_MINIMAP2(ch_post_filter_reads, ch_host_mmi)
        ch_versions      = ch_versions.mix(HOST_DEPLETE_MINIMAP2.out.versions)
        ch_hostdep_reads = HOST_DEPLETE_MINIMAP2.out.reads
        ch_host_stats    = HOST_DEPLETE_MINIMAP2.out.stats_json
    }


    // =====================================================================
    // Step 5.5 — Post-host nanoq
    // =====================================================================
    NANOQ_POSTHOST(ch_hostdep_reads)
    ch_versions = ch_versions.mix(NANOQ_POSTHOST.out.versions)


    // =====================================================================
    // Step 5.6 — INDEX_PANEL (singleton; cached per run)
    // =====================================================================
    ch_panel_fasta = Channel.fromPath(params.reference_panel, checkIfExists: true)
    INDEX_PANEL(ch_panel_fasta)
    ch_versions = ch_versions.mix(INDEX_PANEL.out.versions)

    // Materialise as singleton so it can be `.combine`d / `.first`d into
    // per-sample and per-branch streams without duplicate emission.
    ch_panel = INDEX_PANEL.out.index.first()


    // =====================================================================
    // Step 5.7 — MINIMAP2_ROUND1 competitive mapping
    // =====================================================================
    MINIMAP2_ROUND1(ch_hostdep_reads, ch_panel)
    ch_versions = ch_versions.mix(MINIMAP2_ROUND1.out.versions)

    // Route samples that triggered NO_HCV_DETECTED away from downstream
    // processing.  Use a left-outer join with the optional sentinel channel.
    ch_r1_routed = MINIMAP2_ROUND1.out.bam
        .join(MINIMAP2_ROUND1.out.sentinel, by: 0, remainder: true)
        .branch { meta, bam, bai, sentinel ->
            no_hcv:      sentinel != null
            processable: true
        }

    ch_round1_bam  = ch_r1_routed.processable.map { meta, bam, bai, sentinel -> tuple(meta, bam, bai) }
    ch_no_hcv_r1   = ch_r1_routed.no_hcv.map      { meta, bam, bai, sentinel -> tuple(meta, 'NO_HCV_DETECTED') }


    // =====================================================================
    // Step 5.8 — GENOTYPE_CLASSIFY
    // =====================================================================
    GENOTYPE_CLASSIFY(ch_round1_bam)
    ch_versions = ch_versions.mix(GENOTYPE_CLASSIFY.out.versions)


    // =====================================================================
    // Steps 5.9–5.10 — GENOTYPE_BRANCH fan-out
    // =====================================================================
    GENOTYPE_BRANCH(
        ch_hostdep_reads,
        GENOTYPE_CLASSIFY.out.summary,
        ch_round1_bam,
        GENOTYPE_CLASSIFY.out.assignments,
    )
    ch_versions = ch_versions.mix(GENOTYPE_BRANCH.out.versions)

    ch_branches  = GENOTYPE_BRANCH.out.branches    // [branch_meta, reads, dom_ref_id]
    ch_no_hcv_gb = GENOTYPE_BRANCH.out.no_hcv.map  { meta -> tuple(meta, 'NO_HCV_DETECTED') }

    // Union all NO_HCV channels; downstream EMIT_FAILURE_FLAG publishes a
    // status file per sample.
    ch_failed_samples = Channel.empty()
        .mix(ch_no_hcv_r1)
        .mix(ch_no_hcv_gb)
        .mix(ch_all_filtered_failed)


    // =====================================================================
    // Step 5.10 — BUILD_CONSENSUS (per branch)
    // =====================================================================
    BUILD_CONSENSUS(ch_branches, ch_panel)
    ch_versions = ch_versions.mix(BUILD_CONSENSUS.out.versions)


    // =====================================================================
    // Step 5.11 — MINIMAP2_ROUND2 (full-depth per branch)
    //
    // Join the per-branch consensus with the per-branch reads, both keyed by
    // [meta.id, meta.genotype].
    // =====================================================================
    ch_consensus = BUILD_CONSENSUS.out.consensus
        .map { meta, fasta, fai, mmi -> tuple([meta.id, meta.genotype], meta, fasta, fai, mmi) }

    ch_branch_reads_keyed = ch_branches
        .map { meta, reads, dom_ref -> tuple([meta.id, meta.genotype], reads) }

    ch_round2_input = ch_consensus
        .join(ch_branch_reads_keyed, by: 0)
        .map { key, meta, fasta, fai, mmi, reads ->
            tuple(meta, fasta, fai, mmi, reads)
        }

    MINIMAP2_ROUND2(ch_round2_input)
    ch_versions = ch_versions.mix(MINIMAP2_ROUND2.out.versions)


    // =====================================================================
    // Step 5.13 — MOSDEPTH (coverage QC on full-depth round 2)
    // =====================================================================
    MOSDEPTH(MINIMAP2_ROUND2.out.bam)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)


    // =====================================================================
    // Steps 5.14–5.16 — LoFreq prep / call / filter
    //
    // PREP_LOFREQ_INPUT expects: [meta, bam, bai, consensus_fasta]
    // Feed it the full-depth round 2 BAM directly; RASUSA_ALN inside the
    // subworkflow subsamples to params.lofreq_max_depth using `rasusa aln`,
    // which is coverage-accurate and avoids a redundant minimap2 remap.
    // =====================================================================
    ch_consensus_fasta_only = BUILD_CONSENSUS.out.consensus
        .map { meta, fasta, fai, mmi -> tuple([meta.id, meta.genotype], fasta) }

    ch_round2_bam_keyed = MINIMAP2_ROUND2.out.bam
        .map { meta, bam, bai -> tuple([meta.id, meta.genotype], meta, bam, bai) }

    ch_prep_lofreq_input = ch_round2_bam_keyed
        .join(ch_consensus_fasta_only, by: 0)
        .map { key, meta, bam, bai, fasta -> tuple(meta, bam, bai, fasta) }

    PREP_LOFREQ_INPUT(ch_prep_lofreq_input)
    ch_versions = ch_versions.mix(PREP_LOFREQ_INPUT.out.versions)

    // LoFreq call input: [meta, bam, bai, ref_fasta]
    ch_lofreq_bam_keyed = PREP_LOFREQ_INPUT.out.lofreq_bam
        .map { meta, bam, bai -> tuple([meta.id, meta.genotype], meta, bam, bai) }

    ch_lofreq_call_input = ch_lofreq_bam_keyed
        .join(ch_consensus_fasta_only, by: 0)
        .map { key, meta, bam, bai, fasta -> tuple(meta, bam, bai, fasta) }

    LOFREQ_CALL(ch_lofreq_call_input)
    ch_versions = ch_versions.mix(LOFREQ_CALL.out.versions)

    VARIANT_FILTER(LOFREQ_CALL.out.vcf)
    ch_versions = ch_versions.mix(VARIANT_FILTER.out.versions)


    // =====================================================================
    // Step 5.16b — Round 2 consensus (published output only)
    //
    // Uses the full-depth Round 2 BAM for coverage masking and the LoFreq
    // filtered VCF to generate two consensus FASTAs per branch:
    //   *_round2_consensus_simple.fasta — majority-allele (AF >= 0.5)
    //   *_round2_consensus_iupac.fasta  — IUPAC codes at AF < 0.5 sites
    // These are published under consensus/${GT}/ and do not feed back
    // into any downstream pipeline step.
    // =====================================================================
    MAKE_ROUND2_MASK(MINIMAP2_ROUND2.out.bam)
    ch_versions = ch_versions.mix(MAKE_ROUND2_MASK.out.versions)

    // Join: [meta, vcf, tbi] + [meta, fasta, fai] + [meta, mask_bed]
    // keyed by [meta.id, meta.genotype]
    ch_r2cons_vcf  = VARIANT_FILTER.out.vcf
        .map { meta, vcf, tbi -> tuple([meta.id, meta.genotype], meta, vcf, tbi) }
    ch_r2cons_ref  = BUILD_CONSENSUS.out.consensus
        .map { meta, fasta, fai, mmi -> tuple([meta.id, meta.genotype], fasta, fai) }
    ch_r2cons_mask = MAKE_ROUND2_MASK.out.mask_bed
        .map { meta, bed -> tuple([meta.id, meta.genotype], bed) }

    ch_r2cons_input = ch_r2cons_vcf
        .join(ch_r2cons_ref,  by: 0)
        .join(ch_r2cons_mask, by: 0)
        .map { key, meta, vcf, tbi, fasta, fai, bed ->
            tuple(meta, vcf, tbi, fasta, fai, bed)
        }

    BUILD_ROUND2_CONSENSUS(ch_r2cons_input)
    ch_versions = ch_versions.mix(BUILD_ROUND2_CONSENSUS.out.versions)


    // =====================================================================
    // Step 5.17 — Optional CLAIR3 corroboration
    // =====================================================================
    ch_clair3_vcf = Channel.empty()
    if (params.run_clair3) {
        // CLAIR3 expects: [meta, bam, bai, ref_fasta, ref_fai]
        ch_consensus_fasta_fai = BUILD_CONSENSUS.out.consensus
            .map { meta, fasta, fai, mmi -> tuple([meta.id, meta.genotype], fasta, fai) }

        ch_clair3_input = ch_lofreq_bam_keyed
            .join(ch_consensus_fasta_fai, by: 0)
            .map { key, meta, bam, bai, fasta, fai -> tuple(meta, bam, bai, fasta, fai) }

        CLAIR3(ch_clair3_input)
        ch_versions   = ch_versions.mix(CLAIR3.out.versions)
        ch_clair3_vcf = CLAIR3.out.vcf
    }


    // =====================================================================
    // Steps 5.18–5.20 — DEVIDER prep, run, stitch
    //
    // LOW_COVERAGE branches skip DEVIDER but keep LoFreq.
    // =====================================================================
    PREP_DEVIDER_INPUT(ch_prep_lofreq_input)  // [meta, bam, bai, consensus_fasta] — same shape as LoFreq prep
    ch_versions = ch_versions.mix(PREP_DEVIDER_INPUT.out.versions)

    // Apply DEVIDER-specific AF threshold to the analysis VCF.
    // VARIANT_FILTER outputs PASS variants at >= min_report_af (1%); DEVIDER
    // needs a higher floor (params.devider_min_af, default 5%) to avoid graph
    // saturation from low-AF noise variants.
    FILTER_VCF_FOR_DEVIDER(VARIANT_FILTER.out.vcf)
    ch_versions = ch_versions.mix(FILTER_VCF_FOR_DEVIDER.out.versions)

    // Build DEVIDER input: [meta, bam, bai, vcf, tbi, ref_fasta], but only for
    // branches that did NOT trigger LOW_COVERAGE in MOSDEPTH.
    ch_devider_bam_keyed = PREP_DEVIDER_INPUT.out.devider_bam
        .map { meta, bam, bai -> tuple([meta.id, meta.genotype], meta, bam, bai) }

    ch_filtered_vcf_keyed = FILTER_VCF_FOR_DEVIDER.out.vcf
        .map { meta, vcf, tbi -> tuple([meta.id, meta.genotype], vcf, tbi) }

    ch_consensus_for_devider = BUILD_CONSENSUS.out.consensus
        .map { meta, fasta, fai, mmi -> tuple([meta.id, meta.genotype], fasta) }

    // Optional-emit LOW_COVERAGE flags from MOSDEPTH; key by [id, genotype].
    ch_low_cov_keyed = MOSDEPTH.out.low_cov_flag
        .map { meta, flag -> tuple([meta.id, meta.genotype], flag) }

    // Left-outer join with low-cov flag, then branch on its presence.
    ch_devider_assembled = ch_devider_bam_keyed
        .join(ch_filtered_vcf_keyed, by: 0)
        .join(ch_consensus_for_devider, by: 0)
        .join(ch_low_cov_keyed, by: 0, remainder: true)
        .branch { key, meta, bam, bai, vcf, tbi, fasta, low_cov ->
            adequate:    low_cov == null
            low_coverage: true
        }

    ch_devider_run_input = ch_devider_assembled.adequate
        .map { key, meta, bam, bai, vcf, tbi, fasta, low_cov ->
            tuple(meta, bam, bai, vcf, tbi, fasta)
        }

    DEVIDER(ch_devider_run_input)
    ch_versions = ch_versions.mix(DEVIDER.out.versions)

    FORMAT_HAPLOTYPES(DEVIDER.out.outdir)
    ch_versions = ch_versions.mix(FORMAT_HAPLOTYPES.out.versions)


    // =====================================================================
    // Step 5.21 — SAMPLE_REPORT (per sample, collapse all branches)
    //
    // Collapse per-branch outputs back to per-sample tuples keyed by meta.id.
    // The SAMPLE_REPORT process declares all per-branch inputs as `path` —
    // Nextflow will stage them as a flat list in the work directory.  We
    // pass each list-of-paths as a single channel item per sample.
    // =====================================================================

    // Collapse per-branch [meta, file] channels into per-sample
    // [sample_id, [files...]] using groupTuple.  Branches with no output
    // (e.g. FORMAT_HAPLOTYPES for low-cov samples that skipped DEVIDER) are absent
    // from the channel — groupTuple still emits at least one item per
    // sample as long as at least one branch contributed.
    ch_mosdepth_by_sample = MOSDEPTH.out.summary
        .map { meta, f -> tuple(meta.id, f) }
        .groupTuple(by: 0)

    ch_mosdepth_bed_by_sample = MOSDEPTH.out.regions
        .map { meta, f -> tuple(meta.id, f) }
        .groupTuple(by: 0)

    ch_variants_by_sample = VARIANT_FILTER.out.tsv
        .map { meta, f -> tuple(meta.id, f) }
        .groupTuple(by: 0)

    ch_flagstat_by_sample = MINIMAP2_ROUND2.out.flagstat
        .map { meta, f -> tuple(meta.id, f) }
        .groupTuple(by: 0)

    ch_haplotype_report_by_sample = FORMAT_HAPLOTYPES.out.report
        .map { meta, f -> tuple(meta.id, f) }
        .groupTuple(by: 0)

    ch_nanoplot_by_sample = RAW_QC.out.nanoplot_dirs
        .map { meta, f -> tuple(meta.id, f) }
        .groupTuple(by: 0)

    // Per-sample (not per-branch) channels: genotype summary, nanoq, host stats.
    ch_summary_by_sample   = GENOTYPE_CLASSIFY.out.summary.map { meta, j -> tuple(meta.id, j) }
    ch_nanoq_raw_by_sample = RAW_QC.out.nanoq_jsons.map         { meta, j -> tuple(meta.id, j) }
    ch_nanoq_filt_by_sample= NANOQ_FILT.out.json.map            { meta, j -> tuple(meta.id, j) }
    ch_host_stats_by_sample= ch_host_stats.map                  { meta, j -> tuple(meta.id, j) }

    // Recover meta from ch_reads — always populated for every sample regardless
    // of whether downstream processes (genotype classify, etc.) completed.
    // Using ch_reads as the base ensures the report join chain always has a
    // non-empty starting channel, preventing null-collapse in remainder joins.
    ch_meta_keyed = ch_reads.map { meta, fastq -> tuple(meta.id, meta) }

    // Samples short-circuited before genotype classify (NO_HCV / ALL_READS_FILTERED)
    // are NOT joined into the SAMPLE_REPORT input — they have no genotype summary
    // to render against.  Instead they are routed to EMIT_FAILURE_FLAG below,
    // which publishes a PIPELINE_FLAG.txt under ${params.outdir}/${meta.id}/.

    // We use `.join(..., remainder: true)` extensively so that branches that
    // produced no output (e.g. low_cov samples that skipped DEVIDER, or
    // missing host_stats when --host_reference is not set) still appear as
    // null in the report inputs — the SAMPLE_REPORT process handles null
    // (empty) inputs gracefully via its optional-path declarations.
    //
    // Join order matters: start with the meta channel and add columns one at
    // a time, keying always by sample_id.
    ch_report_input = ch_meta_keyed
        .join(ch_summary_by_sample,    by: 0, remainder: true)
        .join(ch_nanoq_raw_by_sample,  by: 0, remainder: true)
        .join(ch_nanoq_filt_by_sample, by: 0, remainder: true)
        .join(ch_host_stats_by_sample, by: 0, remainder: true)
        .join(ch_mosdepth_by_sample,     by: 0, remainder: true)
        .join(ch_mosdepth_bed_by_sample, by: 0, remainder: true)
        .join(ch_variants_by_sample,     by: 0, remainder: true)
        .join(ch_flagstat_by_sample,     by: 0, remainder: true)
        .join(ch_haplotype_report_by_sample, by: 0, remainder: true)
        .join(ch_nanoplot_by_sample,     by: 0, remainder: true)
        .map { sid, meta, summary, nanoq_raw, nanoq_filt, host_stats,
               mosdepth_files, mosdepth_beds, variant_tsvs, flagstat_files,
               haplotype_reports, nanoplot_dirs ->
            tuple(
                meta,
                summary           ?: [],
                nanoq_raw         ?: [],
                nanoq_filt        ?: [],
                host_stats        ?: [],
                mosdepth_files    ?: [],
                mosdepth_beds     ?: [],
                variant_tsvs      ?: [],
                flagstat_files    ?: [],
                haplotype_reports ?: [],
                nanoplot_dirs  ?: [],
                []  // flag_files placeholder — sentinels published separately
            )
        }
        // Drop samples with no genotype summary (NO_HCV_DETECTED short-circuits).
        // Those samples get a PIPELINE_FLAG.txt via EMIT_FAILURE_FLAG instead.
        .filter { meta, summary, _na, _nf, _hs, _md, _mb, _vt, _fs, _hr, _np, _ff ->
            summary != null && !(summary instanceof List && summary.isEmpty())
        }

    SAMPLE_REPORT(ch_report_input)
    ch_versions = ch_versions.mix(SAMPLE_REPORT.out.versions)


    // =====================================================================
    // EMIT_FAILURE_FLAG — publish a one-line status file for short-circuited
    // samples.  These samples still appear in MultiQC via their raw-QC outputs.
    // =====================================================================
    EMIT_FAILURE_FLAG(ch_failed_samples)


    // =====================================================================
    // Step 5.22 — MULTIQC aggregation across all samples + branches
    //
    // Collect every QC file we want MultiQC to parse, into a single flat
    // channel, then `.collect()` to gate the process on all upstream
    // completions.
    // =====================================================================
    ch_multiqc_files = Channel.empty()
        .mix(RAW_QC.out.nanoplot_dirs.map        { meta, d -> d })
        .mix(RAW_QC.out.nanoq_jsons.map          { meta, j -> j })
        .mix(NANOQ_FILT.out.json.map             { meta, j -> j })
        .mix(NANOQ_POSTHOST.out.json.map         { meta, j -> j })
        .mix(MINIMAP2_ROUND1.out.flagstat.map    { meta, f -> f })
        .mix(MINIMAP2_ROUND2.out.flagstat.map    { meta, f -> f })
        .mix(MOSDEPTH.out.summary.map            { meta, f -> f })
        .mix(MOSDEPTH.out.regions.map            { meta, f -> f })
        .mix(SAMPLE_REPORT.out.multiqc_sample.map { meta, f -> f })
        .mix(SAMPLE_REPORT.out.multiqc_branch.map { meta, f -> f })

    ch_multiqc_config = Channel.value(file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true))

    MULTIQC(ch_multiqc_files.collect(), ch_multiqc_config)
    ch_versions = ch_versions.mix(MULTIQC.out.versions)


    // =====================================================================
    // Run summary — collect all per-sample JSON files and render HTML/JSON.
    // =====================================================================
    RENDER_RUN_SUMMARY(SAMPLE_REPORT.out.json.map { meta, j -> j }.collect())


    // =====================================================================
    // Software versions — collect every versions.yml fragment into a single
    // published manifest.
    // =====================================================================
    DUMP_SOFTWARE_VERSIONS(ch_versions.unique().collect())
}
