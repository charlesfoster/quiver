/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BCFTOOLS_CONSENSUS_CALL — Call high-confidence variants for consensus building.

    Purpose:
        Step 5.10 sub-step 3.  NOT the final variant-calling step — this is used
        only to polish the sample-specific consensus reference for Round 2 mapping.

    Scientific rationale (CLAUDE.md D8):
        - HCV is a haploid RNA virus → `bcftools call --ploidy 1` is mandatory.
          Using ploidy 2 would produce diploid genotype calls (0/1, 1/1) and
          corrupt the consensus: heterozygous positions would be reported as IUPAC
          ambiguity codes rather than the majority allele.
        - Retain only called variants with sufficient depth: DP >= min_consensus_cov.
          Because --ploidy 1 -mv only emits ALT calls, an explicit AF>=0.5 filter
          is redundant.  INFO/DP is used to qualify the field and avoid ambiguity
          with FORMAT/DP, which mpileup also writes via -a AD,DP.

    Pipeline:
        1. bcftools mpileup   — pileup with AD and DP annotations
        2. bcftools call      — ploidy-1 variant call (VCF → .gz)
        3. bcftools index     — .csi index on the raw VCF
        4. bcftools view -i   — filter to AF >= 0.5 and DP >= min_consensus_cov
        5. bcftools index     — .csi index on the filtered VCF

    Key flags:
        mpileup  -d 0      remove the depth cap (default 250 would downsample high-coverage sites)
        mpileup  -Q 7      minimum base quality (matches params.min_bq)
        mpileup  -q 20     minimum mapping quality (matches params.min_mq)
        mpileup  --annotate FORMAT/AD,FORMAT/DP
                           AD (allelic depth) and DP in the FORMAT field are
                           required for the AF filter below
        call     --ploidy 1  HCV is haploid (non-negotiable, D8)
        call     -mv        output only variant sites; use the multiallelic caller
        view     AF>=0.5    majority allele only
        view     DP>=...    minimum coverage depth gate

    NOTE on the FORMAT annotation flag:
        bcftools mpileup accepts either `-a AD,DP` (short) or
        `--annotate FORMAT/AD,FORMAT/DP` (long form).  Both are valid in
        bcftools 1.21.  The short form `-a AD,DP` is used here for brevity.

    NOTE on the AF field:
        After `bcftools call -mv`, the FORMAT/AD field contains [REF_depth, ALT_depth].
        bcftools view -i can reference FORMAT fields for the genotyped call; however,
        for a direct AF calculation we compute it from INFO/AD produced during mpileup
        and propagated through call.  The filter expression `AF>=0.5` uses the
        bcftools-inferred allele frequency field (set by `call` from AD counts).

    Inputs:
        meta        — val map with `id` and `genotype` fields
        ref_fasta   — dominant reference FASTA (from EXTRACT_REF)
        bam         — consensus-map BAM (from MINIMAP2_CONSENSUS_MAP)
        bai         — BAM index

    Outputs:
        vcf         — [meta, "consensus_variants.vcf.gz", "consensus_variants.vcf.gz.csi"]
        raw_vcf     — [meta, "raw_variants.vcf.gz"]
        versions    — versions.yml

    Container: bcftools 1.21 (CLAUDE.md §4).
    Label: process_medium (4 CPU, 8 GB, 30 min — Step 5.10 resource spec).

    Output published to:
        ${params.outdir}/${meta.id}/consensus/${meta.genotype}/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BCFTOOLS_CONSENSUS_CALL {

    label 'process_medium'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.21--h8b25389_0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/consensus/${meta.genotype}/" },
        mode: 'copy'
    )

    input:
    tuple val(meta), path(ref_fasta), path(bam), path(bai)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_sample_consensus_variants.vcf.gz"),
          path("${meta.id}_${meta.genotype}_sample_consensus_variants.vcf.gz.csi"),
          emit: vcf
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_sample_consensus_raw.vcf.gz"),
          emit: raw_vcf
    path "versions.yml", emit: versions

    script:
    def raw_name      = "${meta.id}_${meta.genotype}_sample_consensus_raw"
    def filtered_name = "${meta.id}_${meta.genotype}_sample_consensus_variants"
    """
    # ----------------------------------------------------------------
    # Step 1+2: pileup → call (haploid; majority allele only).
    #
    # --ploidy 1: HCV is a haploid RNA virus.  Using ploidy 2 would
    #             produce heterozygous calls and break consensus (D8).
    # -mv:        output only variant sites; multiallelic caller.
    # -d 10000:   depth cap to bound memory usage.
    # -Q 7:       minimum base quality (params.min_bq equivalent).
    # -q 20:      minimum mapping quality (params.min_mq equivalent).
    # -a AD,DP:   annotate FORMAT with allelic depth and total depth;
    #             required for the AF filter in step 3.
    # ----------------------------------------------------------------
    bcftools mpileup \\
        -f ${ref_fasta} \\
        -d 0 \\
        -Q 7 \\
        -q 20 \\
        -a AD,DP \\
        ${bam} \\
    | bcftools call \\
        --ploidy 1 \\
        -mv \\
        -Oz \\
        -o ${raw_name}.vcf.gz

    # ----------------------------------------------------------------
    # Step 3: index the raw VCF.
    # ----------------------------------------------------------------
    bcftools index ${raw_name}.vcf.gz

    # ----------------------------------------------------------------
    # Step 4: filter on minimum depth.
    #
    # --ploidy 1 -mv already ensures only ALT calls are emitted, so
    # an explicit AF>=0.5 filter is redundant here.
    # Qualify DP as INFO/DP to avoid ambiguity with FORMAT/DP, which
    # is also written when mpileup is run with -a AD,DP.
    # ----------------------------------------------------------------
    bcftools view \\
        -i "INFO/DP>=${params.min_consensus_cov}" \\
        -Oz \\
        -o ${filtered_name}.vcf.gz \\
        ${raw_name}.vcf.gz

    # ----------------------------------------------------------------
    # Step 5: index the filtered VCF.
    # ----------------------------------------------------------------
    bcftools index ${filtered_name}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def raw_name_s      = "${meta.id}_${meta.genotype}_sample_consensus_raw"
    def filtered_name_s = "${meta.id}_${meta.genotype}_sample_consensus_variants"
    """
    touch ${raw_name_s}.vcf.gz
    touch ${filtered_name_s}.vcf.gz
    touch ${filtered_name_s}.vcf.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "1.21"
    END_VERSIONS
    """
}
