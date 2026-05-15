/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREP_DEVIDER_INPUT — Subsample, remap, and preprocess reads for DEVIDER haplotype
    reconstruction.

    Purpose:
        Implements Step 5.18 of the data flow specification.  Chains three modules
        to produce a DEVIDER-ready BAM at ≤ params.devider_max_depth (default 1,000×):

            1. RASUSA          — depth-capped random subsampling of the per-genotype
                                 reads to params.devider_max_depth.
            2. MINIMAP2_ROUND2 — remap the subsampled reads to the per-genotype
                                 consensus (reuses the Round 2 mapping module).
            3. LOFREQ_PREPROCESS — run lofreq indelqual --dindel + lofreq alnqual on
                                   the subsampled BAM.

    Why subsample before DEVIDER (CLAUDE.md D10, D11, docs/architecture_reasoning.md §9):
        DEVIDER's memory usage scales super-linearly with depth.  Benchmarks show
        stable behaviour at 100–2,000×; above that, memory can exhaust even 64 GB
        nodes.  The 1,000× cap is deliberately lower than the LoFreq cap (5,000×)
        because DEVIDER's windowed phasing algorithm does not gain sensitivity from
        additional depth once windows are well-covered.
        The cap is applied BEFORE remapping so that the BAM used for phasing reflects
        the capped depth exactly.  The full-depth BAM from MINIMAP2_ROUND2 (called from
        the main workflow) is preserved for coverage QC and reporting.

    Why run LOFREQ_PREPROCESS for DEVIDER (docs/architecture_reasoning.md §9):
        DEVIDER phases over existing SNPs supplied via the input VCF; it does not
        call variants itself.  However, quality-calibrated BAMs improve DEVIDER's
        internal read-to-haplotype assignment because the tool uses base quality
        information when computing assignment likelihoods.  Using the same indelqual +
        alnqual preprocessing chain as LoFreq ensures the BAM qualities are consistent
        with those used to produce the filtered VCF that DEVIDER consumes.

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

    Output channel (`devider_bam`):
        tuple val(meta),
              path("*_preprocessed.bam"),
              path("*_preprocessed.bam.bai")

    See also:
        docs/data_flow.md              Step 5.18
        docs/implementation_prompts.md Prompt 19
        CLAUDE.md                      D10, D11
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RASUSA             } from '../../modules/local/rasusa/main'
include { MINIMAP2_ROUND2    } from '../../modules/local/minimap2_round2/main'
include { LOFREQ_PREPROCESS  } from '../../modules/local/lofreq_preprocess/main'


workflow PREP_DEVIDER_INPUT {

    take:
    ch_input    // channel: [meta, reads, consensus_fasta, consensus_fai, consensus_mmi]

    main:

    ch_versions = Channel.empty()

    // ----------------------------------------------------------------
    // Step 1: Subsample reads to params.devider_max_depth.
    //
    // RASUSA expects:
    //   tuple val(meta), path(reads), val(coverage), val(genome_size), val(seed)
    //
    // We inject the DEVIDER-specific coverage cap, genome size, and seed
    // from params here.  The genome_size is fixed at 9646 (HCV genome
    // length per the data flow spec) — the same value used for LoFreq
    // subsampling.  Using a different seed (params.rasusa_seed_devider,
    // default 43) than the LoFreq seed (default 42) ensures the two
    // subsampled sets are statistically independent.
    // ----------------------------------------------------------------
    ch_rasusa_input = ch_input.map { meta, reads, consensus_fasta, consensus_fai, consensus_mmi ->
        tuple(
            meta,
            reads,
            params.devider_max_depth,
            9646,
            params.rasusa_seed_devider
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
    devider_bam = LOFREQ_PREPROCESS.out.bam   // [meta, *_preprocessed.bam, *_preprocessed.bam.bai]
    versions    = ch_versions
}
