/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREP_LOFREQ_INPUT — Subsample and preprocess the full-depth round 2 BAM for
    LoFreq variant calling.

    Purpose:
        Implements Step 5.14 of the data flow specification.  Takes the full-depth
        round 2 BAM produced by MINIMAP2_ROUND2 and chains two modules to produce a
        LoFreq-ready BAM at ≤ params.lofreq_max_depth (default 5,000×):

            1. RASUSA_ALN      — depth-accurate subsampling of the full-depth BAM
                                 using `rasusa aln`, which calculates coverage from
                                 actual per-position alignment depth rather than
                                 estimating from read count × genome size.
            2. LOFREQ_PREPROCESS — run lofreq indelqual --dindel + alnqual on the
                                   subsampled BAM.

    Why rasusa aln over rasusa reads (see also prep_devider_input.nf):
        - Coverage is measured from real alignment depth, not estimated.
        - The full-depth BAM is already available from MINIMAP2_ROUND2; subsampling
          it directly avoids a redundant minimap2 remapping step.
        - No genome-size approximation needed (the HCV consensus length varies
          slightly per sample; 9646 is only an estimate).

    Why subsample before variant calling (CLAUDE.md D10, docs/architecture_reasoning.md §8):
        LoFreq sensitivity plateaus at ~1,000–5,000×.  Above ~10,000×, error
        stratification increases the false-positive rate.  The full-depth BAM is
        preserved in the main workflow for coverage QC and reporting.

    Input channel (`ch_input`) expected shape:
        tuple val(meta),
              path(bam),
              path(bai),
              path(consensus_fasta)

    This is formed in the calling workflow from MINIMAP2_ROUND2.out.bam joined
    with the consensus FASTA from BUILD_CONSENSUS.

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

include { RASUSA_ALN         } from '../../modules/local/rasusa_aln/main'
include { LOFREQ_PREPROCESS  } from '../../modules/local/lofreq_preprocess/main'


workflow PREP_LOFREQ_INPUT {

    take:
    ch_input    // channel: [meta, bam, bai, consensus_fasta]

    main:

    ch_versions = Channel.empty()

    // ----------------------------------------------------------------
    // Step 1: Subsample the full-depth BAM to params.lofreq_max_depth.
    //
    // RASUSA_ALN expects:
    //   tuple val(meta), path(bam), path(bai), val(coverage), val(seed)
    // ----------------------------------------------------------------
    ch_rasusa_input = ch_input.map { meta, bam, bai, consensus_fasta ->
        tuple(meta, bam, bai, params.lofreq_max_depth, params.rasusa_seed_lofreq)
    }

    RASUSA_ALN(ch_rasusa_input)
    ch_versions = ch_versions.mix(RASUSA_ALN.out.versions)

    // ----------------------------------------------------------------
    // Step 2: Run LoFreq preprocessing on the subsampled BAM.
    //
    // LOFREQ_PREPROCESS expects:
    //   tuple val(meta), path(bam), path(bai), path(ref_fasta)
    //
    // Join the subsampled BAM with the consensus FASTA from the input channel.
    // ----------------------------------------------------------------
    ch_consensus_fasta = ch_input.map { meta, bam, bai, consensus_fasta ->
        tuple([meta.id, meta.genotype], consensus_fasta)
    }

    ch_subsampled_keyed = RASUSA_ALN.out.bam.map { meta, bam, bai ->
        tuple([meta.id, meta.genotype], meta, bam, bai)
    }

    ch_preprocess_input = ch_subsampled_keyed
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
