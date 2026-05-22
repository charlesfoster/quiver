/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECK_CONSENSUS_DIVERGENCE — Pairwise identity between the sample consensus
    and the dominant panel reference; flags potential novel subtypes.

    Purpose:
        After APPLY_CONSENSUS produces the polished sample-specific consensus,
        this process aligns it back against the extracted dominant panel reference
        using minimap2 in asm5 mode and computes pairwise nucleotide identity.

        If identity falls below params.min_consensus_identity (default 0.90), a
        DIVERGENT_CONSENSUS sentinel is emitted.  The sample continues normally
        downstream — this flag is informational, not a hard stop.

    Why 90% as the default threshold?
        Within-subtype nucleotide variation is typically 2–8%, placing the floor
        of "definitely same subtype" assignments at ~92% identity.  The 90%
        default sits 2% below that floor: it catches sequences that diverge more
        than any normal within-subtype variation would explain, without being so
        tight that it fires for divergent-but-valid strains.  Genotype-6 users
        may wish to lower the threshold to 0.88 via `--min_consensus_identity`.

    Important interaction with LOW_COVERAGE_CONSENSUS:
        If the consensus was heavily N-masked (LOW_COVERAGE_CONSENSUS sentinel
        also fired), computed identity will be artificially low because Ns never
        match the reference.  The SAMPLE_REPORT renderer notes both flags together
        so the user can distinguish "truly divergent" from "just low coverage".

    Algorithm:
        1. minimap2 -x asm5 --cs <ref> <consensus> → PAF
        2. awk: sum column 10 (residue matches) and column 11 (alignment block
           length) across all alignment records.
        3. identity = total_matches / total_alignment_length
        4. coverage = total_alignment_length / query_length
        5. Emit JSON; emit sentinel file when identity < threshold.

    Container: minimap2:2.28 (pure C binary; awk/bash available in container)
    Label:     process_low
    PublishDir: ${params.outdir}/${meta.id}/consensus/${meta.genotype}/

    Input:
        meta             — val map with { id, genotype, ... }
        consensus_fasta  — path to *_consensus.fasta (from APPLY_CONSENSUS)
        ref_fasta        — path to dominant panel reference FASTA (from EXTRACT_REF)

    Output:
        stats    — [meta, "${meta.id}_${meta.genotype}_divergence.json"]
        sentinel — [meta, "${meta.id}_${meta.genotype}.DIVERGENT_CONSENSUS"]
                   optional: present only when identity < threshold
        versions — versions.yml
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CHECK_CONSENSUS_DIVERGENCE {

    label 'process_low'

    tag "${meta.id}:${meta.genotype}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0' :
        'quay.io/biocontainers/minimap2:2.28--he4a0461_0' }"
    conda "${moduleDir}/environment.yml"

    publishDir (
        path: { "${params.outdir}/${meta.id}/consensus/${meta.genotype}/" },
        mode: 'copy',
        pattern: '*_divergence.json'
    )

    input:
    tuple val(meta), path(consensus_fasta), path(ref_fasta)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.genotype}_divergence.json"),
          emit: stats
    tuple val(meta),
          path("${meta.id}_${meta.genotype}.DIVERGENT_CONSENSUS"),
          optional: true,
          emit: sentinel
    path "versions.yml", emit: versions

    script:
    def threshold = params.min_consensus_identity
    """
    # Align consensus to the dominant panel reference using the assembly preset.
    # asm5 is tuned for sequences with < 5% divergence but degrades gracefully
    # above that — it still produces useful alignments up to ~15% divergence.
    minimap2 -x asm5 --cs "${ref_fasta}" "${consensus_fasta}" > aln.paf

    # ---- Parse PAF columns (1-indexed in awk) ----
    # \$2  = query sequence length
    # \$10 = number of residue matches
    # \$11 = alignment block length (matches + mismatches + ref gaps)
    #
    # Sum across all alignment records; a fragmented alignment (multiple
    # records) is handled correctly by accumulating all blocks.
    MATCHES=\$(awk 'NF>=11{m+=\$10} END{print m+0}' aln.paf)
    ALN_LEN=\$(awk 'NF>=11{l+=\$11} END{print l+0}' aln.paf)
    QUERY_LEN=\$(awk 'NR==1&&NF>=2{print \$2; exit} END{if(NR==0) print 0}' aln.paf)

    IDENTITY=\$(awk -v m="\${MATCHES}" -v l="\${ALN_LEN}" \\
        'BEGIN{if(l>0) printf "%.6f", m/l; else print "0.000000"}')
    COVERAGE=\$(awk -v l="\${ALN_LEN}" -v q="\${QUERY_LEN}" \\
        'BEGIN{if(q>0) printf "%.6f", l/q; else print "0.000000"}')
    IS_DIVERGENT=\$(awk -v id="\${IDENTITY}" -v thr="${threshold}" \\
        'BEGIN{print (id+0 < thr+0) ? "true" : "false"}')

    # ---- Emit JSON ----
    cat > ${meta.id}_${meta.genotype}_divergence.json << JSONEOF
{
  "sample_id": "${meta.id}",
  "genotype": "${meta.genotype}",
  "consensus_identity": \${IDENTITY},
  "alignment_coverage": \${COVERAGE},
  "matches": \${MATCHES},
  "alignment_length": \${ALN_LEN},
  "query_length": \${QUERY_LEN},
  "threshold": ${threshold},
  "is_divergent": \${IS_DIVERGENT}
}
JSONEOF

    # ---- Emit sentinel when divergent ----
    if [ "\${IS_DIVERGENT}" = "true" ]; then
        touch "${meta.id}_${meta.genotype}.DIVERGENT_CONSENSUS"
        echo "DIVERGENT_CONSENSUS: ${meta.id} (${meta.genotype}) identity=\${IDENTITY} threshold=${threshold}" >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    cat > ${meta.id}_${meta.genotype}_divergence.json << 'JSONEOF'
{
  "sample_id": "${meta.id}",
  "genotype": "${meta.genotype}",
  "consensus_identity": 0.960000,
  "alignment_coverage": 0.990000,
  "matches": 9200,
  "alignment_length": 9584,
  "query_length": 9646,
  "threshold": 0.9,
  "is_divergent": false
}
JSONEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "2.28"
    END_VERSIONS
    """
}
