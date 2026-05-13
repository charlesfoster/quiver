/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREP_LOFREQ_INPUT — Subsample, remap, and preprocess reads for LoFreq calling.

    Purpose:
        Implements Step 5.14 of the data flow specification.  Chains three modules
        to produce a LoFreq-ready BAM at ≤ params.lofreq_max_depth (default 5,000×):

            1. RASUSA          — depth-capped random subsampling of the per-genotype
                                 reads to params.lofreq_max_depth.
            2. MINIMAP2_ROUND2 — remap the subsampled reads to the per-genotype
                                 consensus (reuses the Round 2 mapping module).
            3. LOFREQ_PREPROCESS — run lofreq indelqual --dindel + lofreq alnqual on
                                   the subsampled BAM.

    Why subsample before variant calling (CLAUDE.md D10, docs/architecture_reasoning.md §8):
        LoFreq sensitivity plateaus at ~1,000–5,000×.  Above ~10,000×, error
        stratification increases the false-positive rate.  The cap is applied
        BEFORE remapping so that the BAM used for calling reflects the capped depth
        exactly.  The full-depth BAM from MINIMAP2_ROUND2 (called from the main
        workflow) is preserved for coverage QC and reporting.

    The MINIMAP2_ROUND2 module is re-used here with subsampled reads.  Its output
    BAM name pattern (*_round2.bam) is unique per meta.id + meta.genotype, so there
    is no filename collision with the full-depth BAM (which lives in the main workflow's
    work directory, not this subworkflow's).

    Input channel (`ch_input`) expected shape:
        tuple val(meta),
              path(reads),
              path(consensus_fasta),
              path(consensus_fai),
              path(consensus_mmi)

    This is formed in the calling workflow by combining the branch reads channel with
    the consensus channel from BUILD_CONSENSUS.

    Output channel (`lofreq_bam`):
        tuple val(meta),
              path("*_preprocessed.bam"),
              path("*_preprocessed.bam.bai")

    See also:
        docs/data_flow.md              Step 5.14
        docs/implementation_prompts.md Prompt 15
        CLAUDE.md                      D10
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RASUSA             } from '../../modules/local/rasusa'
include { MINIMAP2_ROUND2    } from '../../modules/local/minimap2_round2'
include { LOFREQ_PREPROCESS  } from '../../modules/local/lofreq_preprocess'


workflow PREP_LOFREQ_INPUT {

    take:
    ch_input    // channel: [meta, reads, consensus_fasta, consensus_fai, consensus_mmi]

    main:

    ch_versions = Channel.empty()

    // ----------------------------------------------------------------
    // Step 1: Subsample reads to params.lofreq_max_depth.
    //
    // RASUSA expects:
    //   tuple val(meta), path(reads), val(coverage), val(genome_size), val(seed)
    //
    // We inject the coverage cap, genome size, and seed from params here.
    // The genome_size defaults to 9646 (HCV genome length per the data flow spec).
    // Using a fixed 9646 avoids a circular dependency on a per-sample reference length
    // — and the HCV genome is consistently ~9,646 bp across all subtypes.
    // ----------------------------------------------------------------
    ch_rasusa_input = ch_input.map { meta, reads, consensus_fasta, consensus_fai, consensus_mmi ->
        tuple(
            meta,
            reads,
            params.lofreq_max_depth,
            9646,
            params.rasusa_seed_lofreq
        )
    }

    RASUSA(ch_rasusa_input)
    ch_versions = ch_versions.mix(RASUSA.out.versions)

    // ----------------------------------------------------------------
    // Step 2: Remap the subsampled reads to the per-genotype consensus.
    //
    // MINIMAP2_ROUND2 expects:
    //   tuple val(meta), path(consensus_fasta), path(consensus_fai),
    //         path(consensus_mmi), path(reads)
    //
    // Join the subsampled reads (keyed by [meta.id, meta.genotype]) with
    // the consensus files from the input channel.
    //
    // Both streams carry unique [meta.id, meta.genotype] combinations,
    // so a simple key-based join is safe.
    // ----------------------------------------------------------------
    ch_consensus = ch_input.map { meta, reads, consensus_fasta, consensus_fai, consensus_mmi ->
        tuple([meta.id, meta.genotype], consensus_fasta, consensus_fai, consensus_mmi)
    }

    ch_subsampled_keyed = RASUSA.out.reads.map { meta, subsampled_reads ->
        tuple([meta.id, meta.genotype], meta, subsampled_reads)
    }

    ch_round2_input = ch_subsampled_keyed
        .join(ch_consensus, by: 0)
        .map { key, meta, subsampled_reads, consensus_fasta, consensus_fai, consensus_mmi ->
            tuple(meta, consensus_fasta, consensus_fai, consensus_mmi, subsampled_reads)
        }

    MINIMAP2_ROUND2(ch_round2_input)
    ch_versions = ch_versions.mix(MINIMAP2_ROUND2.out.versions)

    // ----------------------------------------------------------------
    // Step 3: Run LoFreq preprocessing on the subsampled BAM.
    //
    // LOFREQ_PREPROCESS expects:
    //   tuple val(meta), path(bam), path(bai), path(ref_fasta)
    //
    // Join the remapped BAM with the consensus FASTA (needed as the
    // -f reference for lofreq indelqual and alnqual).
    // ----------------------------------------------------------------
    ch_consensus_fasta = ch_input.map { meta, reads, consensus_fasta, consensus_fai, consensus_mmi ->
        tuple([meta.id, meta.genotype], consensus_fasta)
    }

    ch_bam_keyed = MINIMAP2_ROUND2.out.bam.map { meta, bam, bai ->
        tuple([meta.id, meta.genotype], meta, bam, bai)
    }

    ch_preprocess_input = ch_bam_keyed
        .join(ch_consensus_fasta, by: 0)
        .map { key, meta, bam, bai, consensus_fasta ->
            tuple(meta, bam, bai, consensus_fasta)
        }

    LOFREQ_PREPROCESS(ch_preprocess_input)
    ch_versions = ch_versions.mix(LOFREQ_PREPROCESS.out.versions)

    emit:
    lofreq_bam = LOFREQ_PREPROCESS.out.bam   // [meta, *_preprocessed.bam, *_preprocessed.bam.bai]
    versions   = ch_versions
}
