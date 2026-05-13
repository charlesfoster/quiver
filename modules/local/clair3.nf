process CLAIR3 {
    tag "${meta.id}_${meta.genotype}"
    label 'process_high'

    container 'quay.io/biocontainers/clair3:1.0.10--py39hf5e1c6e_2'

    publishDir { "${params.outdir}/${meta.id}/variants/${meta.genotype}/clair3/" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(ref_fasta), path(ref_fai)

    output:
    tuple val(meta), path("clair3_output/merge_output.vcf.gz"), path("clair3_output/merge_output.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("clair3_output/"), emit: outdir
    path "versions.yml", emit: versions

    script:
    def model_dir = params.clair3_model_dir ?: '/usr/local/lib/python3.9/dist-packages/clair3/models/r1041_e82_400bps_hac_v520'
    """
    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${ref_fasta} \\
        --threads=${task.cpus} \\
        --platform=ont \\
        --model_path=${model_dir} \\
        --output=clair3_output \\
        --min_mq=${params.min_mq} \\
        --min_coverage=2 \\
        --haploid_precise \\
        --no_phasing_for_fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version 2>&1 | grep -oP '(?<=v)[0-9.]+' || echo "1.0.10")
    END_VERSIONS
    """
}
