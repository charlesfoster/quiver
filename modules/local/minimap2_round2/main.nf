/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MINIMAP2_ROUND2 — Round 2 mapping of per-genotype reads to the sample-specific
    consensus sequence.

    Purpose:
        Implements Step 5.11 of the data flow specification.  Maps per-genotype
        reads against the polished, sample-specific consensus produced by
        BUILD_CONSENSUS (Prompt 11).  This is a single-reference mapping step;
        the competitive flags used in Round 1 are therefore omitted.

    Scientific rationale (docs/architecture_reasoning.md §1–2, §6, CLAUDE.md D1, D3):
        Round 1 used a 238-sequence panel with `-N 5 --secondary=no` for competitive
        assignment.  Round 2 maps against a single per-genotype consensus, so those
        flags would either be no-ops or harmful (secondary suppression has no meaning
        with a single reference).  All other flags are retained:
            -Y      soft-clip supplementary — required by LoFreq indelqual/alnqual
            --MD    MD mismatch tag for IGV and LoFreq
            --eqx   =/X CIGAR ops — required by LoFreq alnqual and bcftools
            -R      read group — required by LoFreq (errors without it)

    Mapping quality check:
        After indexing, samtools flagstat is run and the primary-mapped percentage
        is computed.  If < 90%, a LOW_MAPPING_RATE sentinel file is emitted.
        The process exits 0 — the sentinel is consumed by the reporting subworkflow.

    This module is re-used by PREP_LOFREQ_INPUT (Prompt 15) with downsampled reads
    as input.  The output filename pattern (*_round2.bam) is shared across both
    invocations; the meta.id / meta.genotype combination ensures no filename
    collision.

    Inputs:
        meta            — val map with `id` and `genotype` fields
        consensus_fasta — per-genotype consensus FASTA from APPLY_CONSENSUS
        consensus_fai   — samtools .fai index of the consensus
        consensus_mmi   — minimap2 .mmi index of the consensus
        reads           — per-genotype FASTQ (full-depth or downsampled)

    Outputs:
        bam             — [meta, "*_round2.bam", "*_round2.bam.bai"]
        flagstat        — [meta, "*_round2.flagstat"]
        low_rate_flag   — [meta, "*.LOW_MAPPING_RATE"]  optional
        versions        — versions.yml

    Container: mulled image providing both minimap2 2.28 and samtools 1.21.
    Label: process_high (16 CPU, 16 GB — Step 5.11 resource spec).

    Output published to:
        ${params.outdir}/${meta.id}/mapping/${meta.genotype}/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MINIMAP2_ROUND2 {

    label 'process_high'

    tag "${meta.id}:${meta.genotype}"

    // Mulled container providing minimap2 2.28 + samtools 1.21 in a single image.
    // Same image used by MINIMAP2_ROUND1 and MINIMAP2_CONSENSUS_MAP.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/mapping/${meta.genotype}/" },
        mode: 'copy'
    )

    input:
    tuple val(meta),
          path(consensus_fasta),
          path(consensus_fai),
          path(consensus_mmi),
          path(reads)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_round2.bam"),
          path("${meta.id}_${meta.genotype}_round2.bam.bai"),
          emit: bam
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_round2.flagstat"),
          emit: flagstat
    tuple val(meta),
          path("${meta.id}_${meta.genotype}.LOW_MAPPING_RATE"),
          optional: true,
          emit: low_rate_flag
    path "versions.yml", emit: versions

    script:
    """
    # ----------------------------------------------------------------
    # Round 2 mapping — single-reference, per-genotype consensus.
    #
    # Key differences from Round 1 (D3, docs/architecture_reasoning.md §6):
    #   - No --secondary=no : we are mapping to a single reference, so
    #     secondary suppression is not needed.
    #   - No -N 5           : N-best candidates are only relevant for
    #     competitive multi-reference mapping.
    # Retained flags:
    #   -Y   soft-clip supplementary (required by LoFreq)
    #   --MD mismatch tag (LoFreq alnqual + IGV)
    #   --eqx =/X CIGAR (LoFreq alnqual + bcftools)
    #   -R   read group (LoFreq errors without it)
    #
    # Pipe directly to samtools sort — no intermediate SAM on disk.
    # ----------------------------------------------------------------
    minimap2 \\
        -ax map-ont \\
        -t ${task.cpus} \\
        -Y \\
        --MD \\
        --eqx \\
        -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ONT" \\
        ${consensus_mmi} \\
        ${reads} \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -O bam \\
        -o ${meta.id}_${meta.genotype}_round2.bam \\
        -

    # ----------------------------------------------------------------
    # Index the sorted BAM.
    # ----------------------------------------------------------------
    samtools index -@ ${task.cpus} ${meta.id}_${meta.genotype}_round2.bam

    # ----------------------------------------------------------------
    # Flagstat — written to file for MultiQC and the mapping-rate gate.
    # ----------------------------------------------------------------
    samtools flagstat \\
        -@ ${task.cpus} \\
        ${meta.id}_${meta.genotype}_round2.bam \\
        > ${meta.id}_${meta.genotype}_round2.flagstat

    # ----------------------------------------------------------------
    # Mapping-rate gate.
    #
    # Extract the count of primary-mapped reads and the total primary
    # reads from flagstat.  Compute integer percentage and emit a
    # sentinel if < 90%.  The process exits 0 so the sample continues;
    # the sentinel is consumed by the reporting / branching logic.
    #
    # flagstat lines of interest:
    #   "N + 0 primary mapped ..."  — primary reads that mapped
    #   "N + 0 primary"             — all primary reads (total)
    # ----------------------------------------------------------------
    mapped=\$(grep "primary mapped" ${meta.id}_${meta.genotype}_round2.flagstat | awk '{print \$1}')
    total=\$(grep -E "^[0-9]+ \\+ [0-9]+ primary\$" ${meta.id}_${meta.genotype}_round2.flagstat | awk '{print \$1}')

    if [ "\${total}" -gt 0 ]; then
        pct=\$(python3 -c "print(int(100 * \${mapped} / \${total}))" 2>/dev/null || echo 0)
    else
        pct=0
    fi

    if [ "\${pct}" -lt 90 ]; then
        echo "WARNING: Only \${pct}% of primary reads mapped in Round 2 for sample ${meta.id} genotype ${meta.genotype} (threshold: 90%). Emitting LOW_MAPPING_RATE sentinel." >&2
        touch ${meta.id}_${meta.genotype}.LOW_MAPPING_RATE
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.genotype}_round2.bam
    touch ${meta.id}_${meta.genotype}_round2.bam.bai
    printf "0 + 0 in total (QC-passed reads + QC-failed reads)\\n0 + 0 primary\\n0 + 0 secondary\\n0 + 0 supplementary\\n0 + 0 duplicates\\n0 + 0 primary mapped (0.00%% : N/A)\\n" \\
        > ${meta.id}_${meta.genotype}_round2.flagstat
    touch ${meta.id}_${meta.genotype}.LOW_MAPPING_RATE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "2.28"
        samtools: "1.21"
    END_VERSIONS
    """
}
