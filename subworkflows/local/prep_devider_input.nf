/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREP_DEVIDER_INPUT — Subsample and preprocess the full-depth round 2 BAM for
    DEVIDER haplotype reconstruction.

    Purpose:
        Implements Step 5.18 of the data flow specification.  Takes the full-depth
        round 2 BAM produced by MINIMAP2_ROUND2 and chains two modules to produce a
        DEVIDER-ready BAM at ≤ params.devider_max_depth (default 1,000×):

            1. RASUSA_ALN      — depth-accurate subsampling of the full-depth BAM
                                 using `rasusa aln`, which calculates coverage from
                                 actual per-position alignment depth rather than
                                 estimating from read count × genome size.
            2. LOFREQ_PREPROCESS — run lofreq indelqual --dindel on the subsampled
                                   BAM.

    Why rasusa aln over rasusa reads (see also prep_lofreq_input.nf):
        - Coverage is measured from real alignment depth, not estimated.
        - The full-depth BAM is already available from MINIMAP2_ROUND2; subsampling
          it directly avoids a redundant minimap2 remapping step.
        - No genome-size approximation needed (the HCV consensus length varies
          slightly per sample; 9646 is only an estimate).

    Why subsample before DEVIDER (CLAUDE.md D10, D11, docs/architecture_reasoning.md §9):
        DEVIDER's memory usage scales super-linearly with depth.  Benchmarks show
        stable behaviour at 100–2,000×; above that, memory can exhaust even 64 GB
        nodes.  The 1,000× cap is deliberately lower than the LoFreq cap (5,000×)
        because DEVIDER's windowed phasing algorithm does not gain sensitivity from
        additional depth once windows are well-covered.
        The full-depth BAM is preserved in the main workflow for coverage QC and
        reporting.

    Why run LOFREQ_PREPROCESS for DEVIDER (docs/architecture_reasoning.md §9):
        DEVIDER phases over existing SNPs supplied via the input VCF; it does not
        call variants itself.  However, quality-calibrated BAMs improve DEVIDER's
        internal read-to-haplotype assignment because the tool uses base quality
        information when computing assignment likelihoods.  Using the same indelqual
        preprocessing chain as LoFreq ensures the BAM qualities are consistent
        with those used to produce the filtered VCF that DEVIDER consumes.

    Input channel (`ch_input`) expected shape:
        tuple val(meta),
              path(bam),
              path(bai),
              path(consensus_fasta)

    This is formed in the calling workflow from MINIMAP2_ROUND2.out.bam joined
    with the consensus FASTA from BUILD_CONSENSUS.

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

include { SAMTOOLS_LENGTH_FILTER } from '../../modules/local/samtools_length_filter/main'
include { RASUSA_ALN             } from '../../modules/local/rasusa_aln/main'
include { LOFREQ_PREPROCESS      } from '../../modules/local/lofreq_preprocess/main'


workflow PREP_DEVIDER_INPUT {

    take:
    ch_input    // channel: [meta, bam, bai, consensus_fasta]

    main:

    ch_versions = Channel.empty()

    // ----------------------------------------------------------------
    // Step 1: Filter reads shorter than params.devider_min_read_length.
    //
    // Short reads span too few SNPs to contribute phasing information to
    // DEVIDER's de Bruijn graph.  Removing them before subsampling ensures
    // the depth cap (Step 2) is drawn exclusively from reads long enough
    // to link adjacent SNP positions.
    //
    // Default: 2,000 bp — spans ~6 SNPs given HCV's ~357 bp inter-SNP
    // spacing, producing well-connected graph edges.
    // ----------------------------------------------------------------
    ch_length_filter_input = ch_input.map { meta, bam, bai, consensus_fasta ->
        tuple(meta, bam, bai, params.devider_min_read_length)
    }

    SAMTOOLS_LENGTH_FILTER(ch_length_filter_input)
    ch_versions = ch_versions.mix(SAMTOOLS_LENGTH_FILTER.out.versions)

    // ----------------------------------------------------------------
    // Step 2: Subsample the length-filtered BAM to params.devider_max_depth.
    //
    // RASUSA_ALN expects:
    //   tuple val(meta), path(bam), path(bai), val(coverage), val(seed)
    //
    // A different seed from the LoFreq subsample (params.rasusa_seed_devider,
    // default 43 vs. LoFreq default 42) ensures the two subsampled sets are
    // statistically independent.
    // ----------------------------------------------------------------
    ch_rasusa_input = SAMTOOLS_LENGTH_FILTER.out.bam.map { meta, bam, bai ->
        tuple(meta, bam, bai, params.devider_max_depth, params.rasusa_seed_devider)
    }

    RASUSA_ALN(ch_rasusa_input)
    ch_versions = ch_versions.mix(RASUSA_ALN.out.versions)

    // ----------------------------------------------------------------
    // Step 3: Run LoFreq preprocessing on the subsampled BAM.
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
    devider_bam = LOFREQ_PREPROCESS.out.bam   // [meta, *_preprocessed.bam, *_preprocessed.bam.bai]
    versions    = ch_versions
}
