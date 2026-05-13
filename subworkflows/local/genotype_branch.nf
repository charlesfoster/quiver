/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOTYPE_BRANCH — Per-genotype fan-out subworkflow.

    This subworkflow is the routing hub of the pipeline.  Given the genotype
    summary JSON produced by `classify_genotype.py`, plus the host-depleted
    FASTQ and the Round-1 BAM, it emits one downstream branch per genotype
    that needs to be processed.  Each downstream module (consensus build,
    Round-2 mapping, variant calling, haplotype reconstruction) then operates
    on one branch at a time.

    See:
        docs/data_flow.md            Steps 5.8 (routing note), 5.9, 5.10
        docs/architecture_reasoning.md Section 5 (mixed detection),
                                       Section 12 (failure modes)
        docs/implementation_prompts.md Prompt 10
        CLAUDE.md                    Section 3 D12

    -----------------------------------------------------------------
    Routing logic
    -----------------------------------------------------------------
    1. Parse the JSON to extract `is_mixed`, `branches_to_run`, and the
       per-genotype `top_reference`.
    2. Samples with an empty `branches_to_run` (i.e. NO_HCV_DETECTED or the
       all-ambiguous edge case) are routed to the `no_hcv` emit channel.
       They do not enter any downstream processing.
    3. Single-genotype samples (is_mixed = false, one branch in
       `branches_to_run`) bypass PARTITION_READS — their host-depleted
       FASTQ is passed through verbatim.
    4. Mixed samples invoke PARTITION_READS to split the BAM into one
       FASTQ per genotype, then re-pair each per-genotype FASTQ with its
       branch metadata for downstream consumption.

    -----------------------------------------------------------------
    Channel semantics
    -----------------------------------------------------------------
    The output `branches` channel emits one tuple per (sample × genotype):

        [branch_meta, reads_fastq, dominant_ref_id]

    where:
        branch_meta       = meta + [genotype: <gt>]
                            All original meta fields are preserved (id,
                            metadata, etc.); `genotype` is added so downstream
                            modules can pin per-branch outputs by
                            `meta.id`/`meta.genotype`.
        reads_fastq       = per-genotype FASTQ (partitioned, OR the full
                            host-depleted FASTQ for single-genotype samples).
        dominant_ref_id   = the `top_reference` value from the JSON
                            (e.g. "1a_M62321.1").  The consensus subworkflow
                            (Prompt 11) extracts this reference from the panel
                            FASTA via `samtools faidx`.

    The join between per-branch metadata and per-genotype FASTQs uses
    `[meta.id, genotype]` as the structured key, so multiple samples and
    multiple genotypes within a sample never collide.

    No `.collect()` is used anywhere — every operator streams.  This preserves
    per-sample parallelism on the downstream side.

    Emits:
        branches  — [branch_meta, reads_fastq, dominant_ref_id]  per branch
        no_hcv    — [meta]  for samples with no HCV signal (routed to
                            EMIT_FAILURE_REPORT in Prompt 24)
        versions  — version files from included modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARTITION_READS } from '../../modules/local/partition_reads'


workflow GENOTYPE_BRANCH {

    take:
    ch_reads        // channel: [meta, fastq]                 — host-depleted reads
    ch_summary      // channel: [meta, genotype_summary_json] — classifier output
    ch_round1_bam   // channel: [meta, bam, bai]              — Round-1 BAM
    ch_assignments  // channel: [meta, read_assignments_tsv]  — per-read TSV

    main:

    ch_versions = Channel.empty()

    // ----------------------------------------------------------------
    // 1. Parse the JSON in a closure.  This produces, per sample, a
    //    rich parsed record we can fan out from.
    //
    //    The JsonSlurper().parse(file) form accepts a File / Path
    //    directly — no `.text` round-trip needed.
    // ----------------------------------------------------------------
    ch_parsed = ch_summary.map { meta, json_path ->
        def parsed = new groovy.json.JsonSlurper().parse(json_path as File)
        // Build a fast-lookup map: genotype → top_reference
        def top_ref_by_gt = [:]
        parsed.genotypes.each { entry ->
            top_ref_by_gt[entry.genotype as String] = entry.top_reference as String
        }
        tuple(
            meta,
            [
                is_mixed:        parsed.is_mixed as boolean,
                branches:        (parsed.branches_to_run ?: []).collect { it as String },
                top_ref_by_gt:   top_ref_by_gt,
            ]
        )
    }

    // ----------------------------------------------------------------
    // 2. Route: samples with no branches (NO_HCV or all-ambiguous edge
    //    case) go to the `no_hcv` emit.  Everything else flows on.
    // ----------------------------------------------------------------
    ch_parsed
        .branch { meta, parsed ->
            no_hcv:      parsed.branches.isEmpty()
            processable: true
        }
        .set { ch_routed }

    // ----------------------------------------------------------------
    // 3. Split processable samples by is_mixed.  Mixed samples need
    //    PARTITION_READS; single-genotype samples bypass it.
    //
    //    Note: a sample with is_mixed = false but multiple entries in
    //    branches_to_run is logically impossible per the classifier
    //    (branches_to_run only adds entries that crossed the secondary
    //    threshold, which by definition flips is_mixed = true).  We
    //    still defensively treat the single-branch case as the only
    //    valid non-mixed shape; defensive branching below covers the
    //    pathological case by routing through the mixed path.
    // ----------------------------------------------------------------
    ch_routed.processable
        .branch { meta, parsed ->
            mixed:  parsed.is_mixed || parsed.branches.size() > 1
            single: true
        }
        .set { ch_split }

    // ----------------------------------------------------------------
    // 4a. Single-genotype path.
    //
    //     Emit one branch per sample.  We carry the host-depleted FASTQ
    //     through unchanged.  Join with ch_reads on meta.id so that
    //     samples whose ordering differs across upstream channels still
    //     pair correctly.
    // ----------------------------------------------------------------
    ch_single_keyed = ch_split.single.map { meta, parsed ->
        def gt = parsed.branches[0]
        def dom_ref = parsed.top_ref_by_gt[gt]
        tuple(meta.id, meta, gt, dom_ref)
    }

    ch_reads_keyed = ch_reads.map { meta, fastq -> tuple(meta.id, fastq) }

    ch_single_branches = ch_single_keyed
        .combine(ch_reads_keyed, by: 0)
        .map { sid, meta, gt, dom_ref, fastq ->
            def branch_meta = meta + [genotype: gt]
            tuple(branch_meta, fastq, dom_ref)
        }

    // ----------------------------------------------------------------
    // 4b. Mixed-genotype path.
    //
    //     Build the PARTITION_READS input tuple by joining the Round-1
    //     BAM with the read-assignments TSV on meta.  We also need the
    //     comma-separated list of genotypes (so the script emits exactly
    //     the branches the classifier sentenced).
    //
    //     The PARTITION_READS output channel is then fanned out: one
    //     FASTQ per genotype + the ambiguous pool.  The ambiguous pool
    //     is dropped here (it is published via publishDir but does not
    //     feed downstream processing).
    // ----------------------------------------------------------------
    ch_mixed_meta = ch_split.mixed.map { meta, parsed ->
        // Sort the genotype list (non-destructive copy) to keep the CLI string
        // deterministic — cache-friendly across resumes.
        def sorted_branches = (parsed.branches as List).toSorted()
        tuple(meta.id, meta, sorted_branches, parsed.top_ref_by_gt)
    }

    // Join the BAM and the assignments TSV on meta.id; both upstream
    // modules emit [meta, ...] so we key by meta.id for a clean join.
    ch_bam_keyed = ch_round1_bam.map { meta, bam, bai -> tuple(meta.id, bam, bai) }
    ch_assn_keyed = ch_assignments.map { meta, tsv -> tuple(meta.id, tsv) }

    ch_partition_input = ch_mixed_meta
        .combine(ch_bam_keyed, by: 0)
        .combine(ch_assn_keyed, by: 0)
        .map { sid, meta, branches, top_ref_by_gt, bam, bai, tsv ->
            // Tuple expected by PARTITION_READS:
            //   tuple(meta, bam, bai, tsv)  +  val(genotypes_csv)
            // Pack the routing context into a sidecar map (`bundle`) we
            // can recover after the process returns.
            def bundle = [
                meta:           meta,
                branches:       branches,
                top_ref_by_gt:  top_ref_by_gt,
                genotypes_csv:  branches.join(','),
            ]
            tuple(bundle, tuple(meta, bam, bai, tsv))
        }

    // Run PARTITION_READS.  Build the two parallel channels Nextflow
    // expects from a two-input process by `multiMap`-ing the joined
    // tuple stream.
    ch_partition_input
        .multiMap { bundle, partition_tuple ->
            inputs:    partition_tuple
            genotypes: bundle.genotypes_csv
            bundles:   bundle
        }
        .set { ch_partition_mm }

    PARTITION_READS (
        ch_partition_mm.inputs,
        ch_partition_mm.genotypes
    )

    ch_versions = ch_versions.mix(PARTITION_READS.out.versions)

    // ----------------------------------------------------------------
    // 5. Fan out PARTITION_READS output.
    //
    //    PARTITION_READS.out.reads emits [meta, [fastq, fastq, ...]] —
    //    one tuple per sample with the glob-collected per-genotype and
    //    ambiguous FASTQs.  We need one tuple per (sample × genotype).
    //
    //    Strategy: re-pair each output tuple with its bundle (via
    //    meta.id), then flatMap over the file list inferring the
    //    genotype from the filename (`${sid}.${gt}.fastq.gz`).  The
    //    ambiguous file (`${sid}.ambiguous.fastq.gz`) is dropped.
    // ----------------------------------------------------------------
    ch_bundles_keyed = ch_partition_mm.bundles.map { bundle -> tuple(bundle.meta.id, bundle) }

    ch_partition_keyed = PARTITION_READS.out.reads.map { meta, fastqs ->
        // `fastqs` is a list when the glob matches multiple files; pysam-style
        // singleton when only one.  Normalise to a list.
        def file_list = fastqs instanceof List ? fastqs : [fastqs]
        tuple(meta.id, file_list)
    }

    ch_mixed_branches = ch_partition_keyed
        .combine(ch_bundles_keyed, by: 0)
        .flatMap { sid, file_list, bundle ->
            // For each requested genotype, locate its FASTQ in the list
            // (the partition script always writes `${sid}.${gt}.fastq.gz`).
            bundle.branches.collect { gt ->
                def expected_name = "${bundle.meta.id}.${gt}.fastq.gz".toString()
                def hit = file_list.find { f -> f.getName() == expected_name }
                if (hit == null) {
                    // Defensive: PARTITION_READS always produces this file
                    // (possibly empty) when --genotypes is passed.  If it's
                    // missing we have a bug or a renamed output; fail loud
                    // so the operator notices rather than silently dropping
                    // the branch.
                    def got_names = file_list.collect { it.getName() }
                    throw new RuntimeException(
                        "PARTITION_READS for sample ${bundle.meta.id} did not produce expected output ${expected_name}. Got: ${got_names}"
                    )
                }
                def branch_meta = bundle.meta + [genotype: gt]
                def dom_ref = bundle.top_ref_by_gt[gt]
                tuple(branch_meta, hit, dom_ref)
            }
        }

    // ----------------------------------------------------------------
    // 6. Merge single + mixed branches into a single output channel.
    //    Order doesn't matter — every downstream consumer operates on
    //    one branch at a time and is keyed by (meta.id, meta.genotype).
    // ----------------------------------------------------------------
    ch_branches = ch_single_branches.mix(ch_mixed_branches)

    // ----------------------------------------------------------------
    // 7. NO_HCV channel — just the meta (downstream report module will
    //    read the genotype_summary.json itself).
    // ----------------------------------------------------------------
    ch_no_hcv = ch_routed.no_hcv.map { meta, parsed -> meta }

    emit:
    branches = ch_branches   // [branch_meta, reads_fastq, dominant_ref_id]
    no_hcv   = ch_no_hcv     // [meta]
    versions = ch_versions
}
