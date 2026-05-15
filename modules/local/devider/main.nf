/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEVIDER — Windowed haplotype reconstruction from phased long reads.

    Purpose:
        Implements Step 5.19 of the data flow specification.  Runs DEVIDER on the
        depth-capped, quality-calibrated BAM produced by PREP_DEVIDER_INPUT, using the
        LoFreq filtered VCF as the SNP set over which DEVIDER phases reads.

    Scientific rationale (CLAUDE.md D11, docs/architecture_reasoning.md §9):
        DEVIDER was designed for noisy long reads (ONT, PacBio CLR) and reconstructs
        haplotypes in overlapping windows across the reference, then links windows via
        shared reads.  It does not call variants — it phases over an externally supplied
        SNP set (here: LoFreq filtered VCF).  The windowed mode is DEVIDER's design
        intent; attempting to run a single global window blows up memory and typically
        fails on viral genomes with heterozygous site densities > ~1/500 bp.

    Graceful failure handling:
        DEVIDER can exit non-zero for several legitimate reasons — insufficient depth
        in a window, too few heterozygous sites, or no reads spanning a window.
        These conditions do not constitute pipeline errors; they simply mean no
        haplotypes can be reconstructed for this sample/genotype.

        When DEVIDER exits non-zero, the shell OR-block:
            1. Creates the output directory (so Nextflow output globs still resolve).
            2. Touches `devider.failed` as a named-sentinel consumed by downstream
               modules (Prompt 21 stitching, Prompt 22 sample report).
            3. Writes a human-readable message to stderr.

        The process itself always exits 0.  This keeps the sample alive in the
        pipeline; the per-sample report will note that haplotypes are unavailable.

    DEVIDER v0.0.1 flag notes — verified against `devider --help`:
        -b          Input BAM (short form; v0.0.1 does not accept --bam).
        -r          Reference FASTA (short form; NOT --ref or --reference).
        -v          Input VCF — SNP set for phasing (short form; NOT --vcf).
        -o          Output directory name (short form).
        -O          Overwrite output directory if it exists.
                    Required for Nextflow work-dir re-runs; without -O, DEVIDER
                    exits non-zero if the directory already exists.
        -t          Thread count (short form).
        --preset nanopore-r10
                    R10.4.1 error profile.  Use nanopore-r10, NOT nanopore-r9
                    (R9.4.1 — too lenient for R10) and NOT ont (does not exist).
        --min-cov   Per-window minimum depth floor.  Windows below this are skipped.
        --min-abund Minimum haplotype relative abundance (fraction, e.g. 0.25 = 25%).
        --output-reads
                    Emit haplotype-tagged BAM.  Required by the stitching agent
                    (Step 5.20) which walks reads spanning window junctions.
        --allele-output
                    Write nucleotide alleles (A/C/G/T) rather than 0/1 codes.
                    Required by the stitching script and sample report renderer.
        NO --merge-windows flag in v0.0.1.  Cross-region stitching is Step 5.20.

    Inputs:
        meta            — val map with `id` and `genotype` fields
        bam             — quality-calibrated BAM (from PREP_DEVIDER_INPUT)
        bai             — BAM index
        vcf             — LoFreq filtered VCF (.vcf.gz, tabix-indexed)
        tbi             — tabix index for vcf
        ref_fasta       — per-genotype consensus FASTA (same reference used in mapping)

    Outputs:
        outdir          — [meta, path("devider_output/")] — always emitted; empty if
                          DEVIDER failed (indicated by presence of devider.failed).
        failed_flag     — [meta, path("devider.failed")] — optional; present only when
                          DEVIDER exited non-zero.
        versions        — versions.yml

    Container:  quay.io/biocontainers/devider:0.0.1--ha6fb395_3
    Conda:      bioconda::devider=0.0.1  (includes osx-arm64 build)
    Label:      process_high_memory  (worst case: 16 CPU, 64 GB, 4 hours — Step 5.19)

    Output published to:
        ${params.outdir}/${meta.id}/haplotypes/${meta.genotype}/devider/

    See also:
        docs/data_flow.md              Step 5.19
        docs/implementation_prompts.md Prompt 20
        CLAUDE.md                      D11, Section 4 (DEVIDER version note)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DEVIDER {

    label 'process_high_memory'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/devider:0.0.1--ha6fb395_3' :
        'quay.io/biocontainers/devider:0.0.1--ha6fb395_3' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/haplotypes/${meta.genotype}/devider/" },
        mode: 'copy'
    )

    input:
    tuple val(meta),
          path(bam),
          path(bai),
          path(vcf),
          path(tbi),
          path(ref_fasta)

    output:
    tuple val(meta),
          path("devider_output/"),
          emit: outdir
    tuple val(meta),
          path("devider.failed"),
          optional: true,
          emit: failed_flag
    path "versions.yml", emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Run DEVIDER windowed haplotype reconstruction.
    #
    # All flags are in short form (-b -r -v -o -t) — v0.0.1 does not
    # accept the long-form equivalents.
    #
    # The OR-block ensures the process exits 0 even when DEVIDER fails,
    # by creating the output directory and touching devider.failed.
    # Downstream modules gate on the presence of devider.failed rather
    # than on the process exit code.
    # ----------------------------------------------------------------
    devider \\
        -b ${bam} \\
        -r ${ref_fasta} \\
        -v ${vcf} \\
        -o devider_output \\
        -O \\
        -t ${task.cpus} \\
        --preset nanopore-r10 \\
        --min-cov ${params.devider_min_cov} \\
        --min-abund ${params.devider_min_abund} \\
        --output-reads \\
        --allele-output \\
    || {
        echo "DEVIDER exited non-zero for ${meta.id} ${meta.genotype}" >&2
        echo "No haplotypes will be available for this sample/genotype." >&2
        mkdir -p devider_output
        touch devider.failed
    }

    # ----------------------------------------------------------------
    # Produce a haplotype-tagged BAM using devider's haplotag_bam helper.
    #
    # devider --output-reads writes ids.txt (read → haplotype assignments).
    # haplotag_bam re-tags the original BAM with those assignments so the
    # downstream stitcher (stitch_haplotypes.py) can walk read-spanning
    # junction evidence across DEVIDER windows.
    #
    # Only run when ids.txt was actually produced (devider succeeded).
    # The BAM is placed inside devider_output/ so the stitcher discovers it
    # via its top-level glob — no index needed because pysam uses until_eof.
    # ----------------------------------------------------------------
    if [ -f devider_output/ids.txt ]; then
        haplotag_bam ${bam} -i devider_output/ids.txt \\
            > devider_output/haplotagged.bam \\
        || echo "haplotag_bam failed; stitcher will run without BAM evidence" >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        devider: \$(devider --version 2>&1 | head -1 || echo "0.0.1")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p devider_output
    touch devider_output/stub_haplotypes.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        devider: "0.0.1"
    END_VERSIONS
    """
}
