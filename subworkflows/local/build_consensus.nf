/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BUILD_CONSENSUS — Per-genotype sample-specific consensus building subworkflow.

    Purpose:
        Implements Step 5.10 of the data flow specification.  Given per-genotype
        reads and the full HCV panel, this subworkflow produces a polished
        sample-specific consensus FASTA, a samtools .fai index, and a minimap2
        .mmi index for use in Round 2 mapping.

    Scientific rationale (docs/architecture_reasoning.md §1–2, CLAUDE.md D1, D8):
        Round 2 mapping aligns per-genotype reads against a sample-specific
        consensus rather than the panel reference.  This reduces reference-
        divergence noise in variant calling and haplotype reconstruction.  The
        consensus is built from majority-allele (AF >= 50%) variants called
        against the single dominant panel reference, with low-coverage positions
        (< params.min_consensus_cov, default 10×) masked as N.

    Modules chained (in order):
        1. EXTRACT_REF                  — extract dominant reference from the panel
        2. MINIMAP2_CONSENSUS_MAP       — map per-genotype reads to dominant reference
        3. BCFTOOLS_CONSENSUS_CALL      — call majority-allele variants
        4. MAKE_MASK_BED                — generate low-coverage mask BED via mosdepth
        5. APPLY_CONSENSUS              — apply variants + mask, rename header
        6. INDEX_CONSENSUS              — build samtools .fai and minimap2 .mmi indices
        7. CHECK_CONSENSUS_DIVERGENCE   — pairwise identity vs panel ref; novel-subtype flag

    Input channels:
        ch_branch  — [branch_meta, reads_fastq, dominant_ref_id]
                     Emitted by GENOTYPE_BRANCH (Prompt 10).
                     branch_meta carries both `id` (sample) and `genotype` fields.
        ch_panel   — [panel_mmi, panel_fasta, panel_fai]
                     Singleton emitted by INDEX_PANEL (Prompt 6).

    Output channels:
        consensus          — [meta, consensus_fasta, consensus_fai, consensus_mmi]
                             Per-branch; consumed by MINIMAP2_ROUND2 (Prompt 12).
        low_cov_sentinel   — [meta, LOW_COVERAGE_CONSENSUS file]  optional
        divergent_sentinel — [meta, DIVERGENT_CONSENSUS file]  optional
        divergence_stats   — [meta, *_divergence.json]
        versions           — version files from all included processes

    Channel join strategy:
        ch_branch carries `dominant_ref_id` (a val) — it is combined with
        ch_panel (a singleton) so every branch gets the full panel without
        requiring a `groupTuple`.

    No `.collect()` is used — the subworkflow stays streaming to preserve
    per-sample parallelism.

    See also:
        docs/data_flow.md              Step 5.10
        docs/implementation_prompts.md Prompt 11
        CLAUDE.md                      D1, D8
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXTRACT_REF                  } from '../../modules/local/extract_ref/main'
include { MINIMAP2_CONSENSUS_MAP       } from '../../modules/local/minimap2_consensus_map/main'
include { BCFTOOLS_CONSENSUS_CALL      } from '../../modules/local/bcftools_consensus_call/main'
include { MAKE_MASK_BED                } from '../../modules/local/make_mask_bed/main'
include { APPLY_CONSENSUS              } from '../../modules/local/apply_consensus/main'
include { INDEX_CONSENSUS              } from '../../modules/local/index_consensus/main'
include { CHECK_CONSENSUS_DIVERGENCE   } from '../../modules/local/check_consensus_divergence/main'


workflow BUILD_CONSENSUS {

    take:
    ch_branch   // channel: [branch_meta, reads_fastq, dominant_ref_id]
    ch_panel    // channel: [panel_mmi, panel_fasta, panel_fai]  (singleton)

    main:

    ch_versions = Channel.empty()

    // ----------------------------------------------------------------
    // Step 1: Extract the dominant reference from the panel FASTA.
    //
    // ch_panel is a singleton value-channel (emitted by INDEX_PANEL via
    // `.first()` in the main workflow).  We combine it with ch_branch so
    // that each branch gets its own copy of the panel file paths.
    //
    // Input to EXTRACT_REF:
    //   tuple val(meta), path(panel_fasta), path(panel_fai), val(ref_id)
    // ----------------------------------------------------------------
    ch_extract_input = ch_branch
        .combine(ch_panel)
        .map { branch_meta, reads, dominant_ref_id, panel_mmi, panel_fasta, panel_fai ->
            // Drop panel_mmi — EXTRACT_REF only needs the FASTA + FAI.
            tuple(branch_meta, panel_fasta, panel_fai, dominant_ref_id)
        }

    EXTRACT_REF(ch_extract_input)
    ch_versions = ch_versions.mix(EXTRACT_REF.out.versions)

    // ----------------------------------------------------------------
    // Step 2: Map per-genotype reads to the dominant reference.
    //
    // Re-join the reads from ch_branch with the extracted ref FASTA.
    // Key by [meta.id, meta.genotype] to avoid cross-sample collisions.
    //
    // Input to MINIMAP2_CONSENSUS_MAP:
    //   tuple val(meta), path(ref_fasta), path(reads)
    // ----------------------------------------------------------------
    ch_reads_keyed = ch_branch.map { meta, reads, dominant_ref_id ->
        tuple([meta.id, meta.genotype], reads)
    }

    ch_ref_keyed = EXTRACT_REF.out.ref_fasta.map { meta, ref_fasta ->
        tuple([meta.id, meta.genotype], meta, ref_fasta)
    }

    ch_map_input = ch_ref_keyed
        .join(ch_reads_keyed, by: 0)
        .map { key, meta, ref_fasta, reads ->
            tuple(meta, ref_fasta, reads)
        }

    MINIMAP2_CONSENSUS_MAP(ch_map_input)
    ch_versions = ch_versions.mix(MINIMAP2_CONSENSUS_MAP.out.versions)

    // ----------------------------------------------------------------
    // Step 3: Call high-confidence variants for consensus building.
    //
    // Join the extracted ref FASTA with the consensus-map BAM.
    // Key by [meta.id, meta.genotype].
    //
    // Input to BCFTOOLS_CONSENSUS_CALL:
    //   tuple val(meta), path(ref_fasta), path(bam), path(bai)
    // ----------------------------------------------------------------
    ch_bam_keyed = MINIMAP2_CONSENSUS_MAP.out.bam.map { meta, bam, bai ->
        tuple([meta.id, meta.genotype], bam, bai)
    }

    ch_call_input = ch_ref_keyed
        .join(ch_bam_keyed, by: 0)
        .map { key, meta, ref_fasta, bam, bai ->
            tuple(meta, ref_fasta, bam, bai)
        }

    BCFTOOLS_CONSENSUS_CALL(ch_call_input)
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS_CALL.out.versions)

    // ----------------------------------------------------------------
    // Step 4: Generate the low-coverage mask BED.
    //
    // Input to MAKE_MASK_BED:
    //   tuple val(meta), path(bam), path(bai)
    // ----------------------------------------------------------------
    MAKE_MASK_BED(MINIMAP2_CONSENSUS_MAP.out.bam)
    ch_versions = ch_versions.mix(MAKE_MASK_BED.out.versions)

    // ----------------------------------------------------------------
    // Step 5: Apply variants + mask, rename header.
    //
    // Join: ref_fasta + vcf (filtered) + vcf_csi + mask_bed.
    // All keyed by [meta.id, meta.genotype].
    //
    // Input to APPLY_CONSENSUS:
    //   tuple val(meta), path(ref_fasta), path(vcf), path(vcf_csi), path(mask_bed)
    // ----------------------------------------------------------------
    ch_vcf_keyed = BCFTOOLS_CONSENSUS_CALL.out.vcf.map { meta, vcf, vcf_csi ->
        tuple([meta.id, meta.genotype], vcf, vcf_csi)
    }

    ch_mask_keyed = MAKE_MASK_BED.out.mask_bed.map { meta, mask_bed ->
        tuple([meta.id, meta.genotype], mask_bed)
    }

    ch_apply_input = ch_ref_keyed
        .join(ch_vcf_keyed,  by: 0)
        .join(ch_mask_keyed, by: 0)
        .map { key, meta, ref_fasta, vcf, vcf_csi, mask_bed ->
            tuple(meta, ref_fasta, vcf, vcf_csi, mask_bed)
        }

    APPLY_CONSENSUS(ch_apply_input)
    ch_versions = ch_versions.mix(APPLY_CONSENSUS.out.versions)

    // ----------------------------------------------------------------
    // Step 6: Build samtools .fai and minimap2 .mmi indices.
    //
    // Input to INDEX_CONSENSUS:
    //   tuple val(meta), path(consensus_fasta)
    //
    // Output: [meta, fasta, fai, mmi]
    // ----------------------------------------------------------------
    INDEX_CONSENSUS(APPLY_CONSENSUS.out.consensus)
    ch_versions = ch_versions.mix(INDEX_CONSENSUS.out.versions)

    // ----------------------------------------------------------------
    // Step 7: Check pairwise identity between sample consensus and the
    //         dominant panel reference.
    //
    // minimap2 asm5 alignment → identity = matches / alignment_length.
    // Emits DIVERGENT_CONSENSUS sentinel when identity < params.min_consensus_identity.
    //
    // Join: consensus FASTA (from APPLY_CONSENSUS) + ref FASTA (from
    // EXTRACT_REF), keyed by [meta.id, meta.genotype].
    // ----------------------------------------------------------------
    ch_consensus_for_div = APPLY_CONSENSUS.out.consensus.map { meta, fasta ->
        tuple([meta.id, meta.genotype], meta, fasta)
    }

    ch_ref_for_div = EXTRACT_REF.out.ref_fasta.map { meta, ref_fasta ->
        tuple([meta.id, meta.genotype], ref_fasta)
    }

    ch_divergence_input = ch_consensus_for_div
        .join(ch_ref_for_div, by: 0)
        .map { key, meta, consensus_fasta, ref_fasta ->
            tuple(meta, consensus_fasta, ref_fasta)
        }

    CHECK_CONSENSUS_DIVERGENCE(ch_divergence_input)
    ch_versions = ch_versions.mix(CHECK_CONSENSUS_DIVERGENCE.out.versions)

    emit:
    consensus           = INDEX_CONSENSUS.out.consensus                    // [meta, fasta, fai, mmi]
    low_cov_sentinel    = APPLY_CONSENSUS.out.low_cov_sentinel             // [meta, sentinel] optional
    divergent_sentinel  = CHECK_CONSENSUS_DIVERGENCE.out.sentinel          // [meta, sentinel] optional
    divergence_stats    = CHECK_CONSENSUS_DIVERGENCE.out.stats             // [meta, json]
    versions            = ch_versions
}
