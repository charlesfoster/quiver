/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT_CHECK — Validate samplesheet CSV and emit sample channel.

    Wraps bin/check_samplesheet.py which validates:
      - Required columns: sample_id, fastq (metadata_json is optional)
      - sample_id character set: ^[A-Za-z0-9._-]+$
      - No duplicate sample_ids
      - FASTQ exists, is readable, is non-empty

    Emits:
      reads  — channel of [meta, file(fastq)] tuples
                 meta = [id: <sample_id>, metadata: <map>]

    On validation failure the Python script exits non-zero and Nextflow
    propagates the error immediately (fail-fast).

    Sentinel: if a FASTQ is 0 bytes the script exits non-zero with a clear
    message — the caller should never reach the sentinel path via this module.
    Empty-FASTQ detection happens at the Python level so the pipeline stops
    early with a descriptive error rather than producing an empty output channel.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process INPUT_CHECK {

    label 'process_low'

    container 'python:3.11-slim'
    conda 'conda-forge::python=3.11'

    input:
    path samplesheet

    output:
    // The process publishes nothing — its only product is the parsed JSON
    // which is read back into Nextflow via a channel transformation in the
    // calling workflow.  We emit a single path so Nextflow tracks the work
    // directory correctly.
    path 'validated_samples.json', emit: json

    script:
    """
    check_samplesheet.py ${samplesheet} > validated_samples.json
    """

    stub:
    """
    echo '[]' > validated_samples.json
    """
}
