/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BUILD_ROUND2_CONSENSUS — Apply variants with masking → simple + IUPAC consensus FASTAs.

    Generates two published output FASTAs per branch:
        *_round2_consensus_simple.fasta  — majority-allele consensus
        *_round2_consensus_iupac.fasta   — IUPAC ambiguity consensus

    Neither output feeds back into the pipeline — both are published outputs only.

    Approach for IUPAC codes:
        LoFreq 2.1.5 emits 8-column VCF (no FORMAT/SAMPLE columns).  We add a
        synthetic FORMAT=GT column with GT=1/1 for all called variants, then use
        `bcftools +setGT` to demote sites with INFO/AF < 0.5 to GT=0/1.
        `bcftools consensus -H I` then emits IUPAC codes at heterozygous positions.

    Container: bcftools 1.21
    Label: process_low

    Published to:
        ${params.outdir}/${meta.id}/consensus/${meta.genotype}/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BUILD_ROUND2_CONSENSUS {

    label 'process_low'

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
    tuple val(meta),
          path(filtered_vcf),
          path(filtered_tbi),
          path(ref_fasta),
          path(ref_fai),
          path(mask_bed)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_round2_consensus_simple.fasta"),
          path("${meta.id}_${meta.genotype}_round2_consensus_iupac.fasta"),
          emit: consensus
    path "versions.yml", emit: versions

    script:
    def simple_name = "${meta.id}_${meta.genotype}_round2_consensus_simple"
    def iupac_name  = "${meta.id}_${meta.genotype}_round2_consensus_iupac"
    """
    # ------------------------------------------------------------------
    # Step 1: Add FORMAT/GT column to LoFreq VCF.
    # LoFreq 2.1.5 emits 8-column VCF (no FORMAT/SAMPLE columns).
    # We add GT=1/1 to every variant record so that bcftools +setGT
    # and bcftools consensus can operate on them.
    # The awk guard handles any future LoFreq version that already
    # emits FORMAT/GT by checking the column count.
    # ------------------------------------------------------------------
    zcat ${filtered_vcf} \\
    | awk '
        BEGIN { has_gt = 0 }
        /^##FORMAT=<ID=GT/ { has_gt = 1 }
        /^##/ { print; next }
        /^#CHROM/ {
            if (!has_gt) {
                print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">"
            }
            if (NF == 8) print \$0 "\\tFORMAT\\tSAMPLE"
            else         print \$0
            next
        }
        NF == 8 { print \$0 "\\tGT\\t1/1"; next }
        { print }
    ' | bgzip -c > with_gt.vcf.gz
    tabix -p vcf with_gt.vcf.gz

    # ------------------------------------------------------------------
    # Step 2: Simple consensus VCF — retain only majority-allele sites
    # (AF >= 0.5).  All retained records have GT=1/1; bcftools consensus
    # applies them as ALT alleles, producing a pure majority consensus.
    # ------------------------------------------------------------------
    bcftools view -i 'INFO/AF>=0.5' with_gt.vcf.gz -Oz -o simple.vcf.gz
    tabix -p vcf simple.vcf.gz

    # ------------------------------------------------------------------
    # Step 3: IUPAC consensus VCF — demote minority-frequency sites
    # (AF < 0.5) to GT=0/1.  bcftools consensus -H I will then emit
    # the appropriate IUPAC ambiguity code at each such site.
    # Sites with AF >= 0.5 retain GT=1/1 and are applied as ALT.
    # ------------------------------------------------------------------
    bcftools +setGT \\
        -Oz -o iupac.vcf.gz \\
        with_gt.vcf.gz \\
        -- -t q -i 'INFO/AF<0.5' -n 'c:0/1'
    tabix -p vcf iupac.vcf.gz

    # ------------------------------------------------------------------
    # Step 4a: Build simple majority-allele consensus.
    # ------------------------------------------------------------------
    bcftools consensus \\
        -f ${ref_fasta} \\
        -m ${mask_bed} \\
        -o simple_raw.fasta \\
        simple.vcf.gz

    awk 'NR==1{print ">${simple_name}"; next} 1' simple_raw.fasta > ${simple_name}.fasta

    # ------------------------------------------------------------------
    # Step 4b: Build IUPAC ambiguity consensus.
    # -H I: use IUPAC codes at heterozygous (0/1) positions.
    # ------------------------------------------------------------------
    bcftools consensus \\
        -H I \\
        -f ${ref_fasta} \\
        -m ${mask_bed} \\
        -o iupac_raw.fasta \\
        iupac.vcf.gz

    awk 'NR==1{print ">${iupac_name}"; next} 1' iupac_raw.fasta > ${iupac_name}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def simple_name_s = "${meta.id}_${meta.genotype}_round2_consensus_simple"
    def iupac_name_s  = "${meta.id}_${meta.genotype}_round2_consensus_iupac"
    """
    printf '>${simple_name_s}\\nACGTACGT\\n' > ${simple_name_s}.fasta
    printf '>${iupac_name_s}\\nACGYMRWS\\n' > ${iupac_name_s}.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "1.21"
    END_VERSIONS
    """
}
