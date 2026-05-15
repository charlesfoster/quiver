#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    hcv-quasi
    HCV quasispecies pipeline: low-frequency variant calling and haplotype
    reconstruction from ONT PromethION reads.
    Documentation: CLAUDE.md and docs/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HCV_QUASI } from './workflows/hcv_quasi'

// -----------------------------------------------------------------------
// Help text
// -----------------------------------------------------------------------
def helpMessage() {
    log.info """
    =========================================
    hcv-quasi  v${workflow.manifest.version}
    =========================================
    HCV quasispecies pipeline: low-frequency variant calling and haplotype
    reconstruction from ONT PromethION reads. Handles mixed-genotype infections
    via automatic detection and per-genotype branching.

    Usage:
        nextflow run main.nf --input samplesheet.csv [options]
        nextflow run main.nf --input samplesheet.csv -profile docker [options]

    Required:
        --input FILE                 Samplesheet CSV with columns: sample_id,fastq,metadata_json

    Inputs:
        --reference_panel FILE       HCV reference panel FASTA [default: assets/hcv_references.fasta]
        --host_reference FILE        GRCh38 FASTA or .mmi for host depletion [default: null]
        --outdir DIR                 Output directory [default: results]

    Read filtering:
        --min_length INT             Minimum read length in bp [default: 200]
        --max_length INT             Maximum read length in bp [default: 10000]
        --min_qual INT               Minimum mean Phred quality score [default: 8]

    Genotyping:
        --min_secondary_fraction FLOAT  Fraction threshold for mixed infection flag [default: 0.05]
        --ambiguous_delta_as INT        AS score gap to flag a read as ambiguous [default: 20]
        --min_round1_mapped INT         Minimum mapped reads to proceed [default: 100]

    Coverage:
        --min_mean_coverage INT      Mean coverage below this triggers LOW_COVERAGE flag [default: 100]
        --min_consensus_cov INT      Per-position depth below this masks with N [default: 10]
        --min_variant_depth INT      LoFreq DP filter threshold [default: 20]

    Variant calling:
        --min_call_af FLOAT          LoFreq lower bound for raw calls [default: 0.005]
        --min_report_af FLOAT        Reporting AF threshold [default: 0.01]
        --min_mq INT                 Minimum mapping quality [default: 20]
        --min_bq INT                 Minimum base quality [default: 7]
        --lofreq_sig FLOAT           LoFreq significance threshold [default: 0.01]
        --lofreq_pp_threads INT      LoFreq parallel workers; 1 uses serial call [default: 8]

    Depth normalisation:
        --lofreq_max_depth INT       rasusa depth cap for LoFreq [default: 5000]
        --devider_max_depth INT      rasusa depth cap for DEVIDER [default: 1000]
        --rasusa_seed_lofreq INT     Random seed for rasusa (LoFreq) [default: 42]
        --rasusa_seed_devider INT    Random seed for rasusa (DEVIDER) [default: 43]

    Haplotype reconstruction:
        --devider_min_cov INT        DEVIDER --min-cov: per-window depth floor [default: 50]
        --devider_min_abund FLOAT    DEVIDER --min-abund: minimum haplotype abundance [default: 0.25]
        --stitch_min_reads INT       Minimum spanning reads to stitch windows [default: 5]

    Run-mode toggles:
        --run_clair3                 Run Clair3 corroboration calling [default: false]
        --use_hostile                Use hostile for host depletion [default: false]

    Resources:
        --max_cpus INT               Maximum CPUs per process [default: 16]
        --max_memory STR             Maximum memory per process [default: 64.GB]
        --max_time STR               Maximum wall time per process [default: 24.h]

    Profiles:
        -profile docker              Local Docker
        -profile docker_mac          Apple Silicon Docker; serial LoFreq for stability
        -profile conda               Local conda/micromamba
        -profile katana              UNSW Katana HPC, SLURM + Singularity
        -profile gadi                NCI Gadi HPC, SLURM + Singularity (requires --gadi_project)
        -profile test                Synthetic minimal test dataset

    Full specification: CLAUDE.md and docs/
    """.stripIndent()
}

// -----------------------------------------------------------------------
// Entrypoint workflow
// -----------------------------------------------------------------------
workflow {

    // Print help and exit cleanly
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Validate required parameter
    if (!params.input) {
        log.error """
        =========================================
        ERROR: --input is required.
        =========================================
        Please provide a samplesheet CSV via --input.
        Example:
            nextflow run main.nf --input samplesheet.csv -profile docker

        Run with --help for full parameter documentation.
        """.stripIndent()
        exit 1
    }

    // Launch main pipeline workflow
    HCV_QUASI()
}

// -----------------------------------------------------------------------
// Completion handler
// -----------------------------------------------------------------------
workflow.onComplete {
    def ms      = workflow.duration.toMillis()
    def hours   = (ms / 3600000) as int
    def mins    = ((ms % 3600000) / 60000) as int
    def secs    = ((ms % 60000) / 1000) as int
    def elapsed = String.format("%d h %02d min %02d s", hours, mins, secs)

    if (workflow.success) {
        log.info """
        =========================================
        hcv-quasi  v${workflow.manifest.version}
        =========================================
        Analysis successfully completed in ${elapsed}.
        Results located in: ${params.outdir}
        =========================================
        """.stripIndent()
    } else {
        log.error """
        =========================================
        hcv-quasi  v${workflow.manifest.version}
        =========================================
        Pipeline completed with errors after ${elapsed}.
        Check the Nextflow log and work/ directory for details.
        =========================================
        """.stripIndent()
    }
}
