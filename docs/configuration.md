# Configuration — HCV Quasispecies Pipeline

## Parameters (`nextflow.config`)

All parameters live under `params { }`. Defaults tuned for the dataset described in the design brief.

```groovy
params {
    // --- Inputs ---
    input               = null   // samplesheet CSV (required)
    reference_panel     = "${projectDir}/assets/hcv_references.fasta"
    host_reference      = null   // local path to GRCh38 fasta or .mmi
    outdir              = "results"

    // --- Read filtering ---
    min_length          = 200
    max_length          = 10000
    min_qual            = 8

    // --- Genotyping / mixed infection ---
    min_secondary_fraction  = 0.05   // ≥5% of reads → mixed flag
    ambiguous_delta_as      = 20     // AS gap to flag a read as ambiguous
    min_round1_mapped       = 100    // minimum mapped reads to proceed

    // --- Coverage thresholds ---
    min_mean_coverage       = 100    // below → skip DEVIDER, flag LOW_COVERAGE
    min_consensus_cov       = 10     // below → mask with N in consensus
    min_variant_depth       = 20     // LoFreq DP filter

    // --- Variant calling ---
    min_call_af             = 0.005  // LoFreq lower bound for raw calls
    min_report_af           = 0.01   // Reporting threshold (1%)
    min_mq                  = 20
    min_bq                  = 7
    min_alt_bq              = 7
    lofreq_sig              = 0.01
    lofreq_pp_threads       = 8

    // --- Depth normalisation ---
    lofreq_max_depth        = 5000
    devider_max_depth       = 5000
    rasusa_seed_lofreq      = 42
    rasusa_seed_devider     = 43

    // --- Haplotype ---
    devider_min_cov         = 20      // DEVIDER --min-cov (haplotype depth floor)
    devider_min_abund       = 0.25    // DEVIDER --min-abund (% abundance floor)
    stitch_min_reads        = 5       // post-hoc stitcher minimum spanning reads

    // --- Run-mode toggles ---
    run_clair3              = false
    allow_conda_fallback    = false
    use_hostile             = false

    // --- Resources (overridable per-profile) ---
    max_cpus                = 16
    max_memory              = '64.GB'
    max_time                = '24.h'
}
```

## Samplesheet schema

```
sample_id,fastq,metadata_json
P001,/path/to/P001.fastq.gz,
P002,/path/to/P002.fastq.gz,{"collection_date":"2026-01-15"}
```

Constraints: `sample_id` must match `^[A-Za-z0-9._-]+$`; must be unique; FASTQ must be readable.

## `genotype_summary.json` schema (Step 5.8 output)

```json
{
  "sample_id": "P001",
  "total_mapped_reads": 853221,
  "ambiguous_reads": 12000,
  "ambiguous_fraction": 0.014,
  "genotypes": [
    {"genotype": "1", "fraction": 0.92, "reads": 785000, "top_subtype": "1a", "top_reference": "1a_M62321.1"},
    {"genotype": "3", "fraction": 0.07, "reads": 60000, "top_subtype": "3a", "top_reference": "3a_D17763.1"}
  ],
  "is_mixed": true,
  "primary_genotype": "1",
  "secondary_genotypes": ["3"],
  "branches_to_run": ["1", "3"]
}
```

---

## Compute Profiles

Select via `-profile <name>`. Profiles defined in `conf/`.

### `docker` — Local Docker

```groovy
process {
    executor = 'local'
    cpus   = { Math.min(task.cpus, 12) }
    memory = { Math.min(task.memory.toGiga(), 56).GB }
}
docker {
    enabled    = true
    runOptions = '--platform=linux/amd64'  // Rosetta fallback if no arm64 image
}
params.max_cpus   = 12
params.max_memory = '56.GB'
```

> Leave 8 GB headroom for macOS. Biocontainers are increasingly multi-arch; if a tool is x86-only, Docker uses Rosetta 2 (~30–50% slower). For Apple Silicon Docker runs prefer the `docker_mac` profile.

### `docker_mac` — Apple Silicon Docker

The `docker_mac` profile inherits `docker`, sets `params.lofreq_pp_threads = 1`, and constrains `LOFREQ_CALL` to one CPU and one fork. This keeps LoFreq serial under Docker Desktop/Rosetta to avoid `call-parallel` OOM/SIGKILL failures.

### `conda` — Local conda fallback

```groovy
process {
    executor = 'local'
}
conda {
    enabled       = true
    useMicromamba = true
}
```

### `katana` — UNSW Katana (SLURM)

```groovy
process {
    executor       = 'slurm'
    queue          = { task.memory > 64.GB ? 'high_mem' : 'normal' }
    clusterOptions = '--account=oz000'   // user-overridable
    scratch        = '$TMPDIR'
}
singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = "/srv/scratch/${USER}/.singularity_cache"
}
params.max_cpus   = 32
params.max_memory = '256.GB'
params.max_time   = '48.h'
```

### `gadi` — NCI Gadi (SLURM)

```groovy
process {
    executor       = 'slurm'
    queue          = 'normal'
    clusterOptions = { "-P ${params.gadi_project} -l storage=scratch/${params.gadi_project}+gdata/${params.gadi_project}" }
    scratch        = '$PBS_JOBFS'
}
singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = "/scratch/${params.gadi_project}/${USER}/.singularity_cache"
}
params.gadi_project = null   // required at runtime: --gadi_project xyz
params.max_cpus     = 48
params.max_memory   = '192.GB'
params.max_time     = '48.h'
```

### `test` — synthetic minimal test

```groovy
params {
    input           = "${projectDir}/test_data/samplesheet.csv"
    reference_panel = "${projectDir}/test_data/mini_panel.fasta"
    host_reference  = "${projectDir}/test_data/mini_host.fasta"
    max_cpus        = 4
    max_memory      = '8.GB'
    max_time        = '2.h'
}
```
