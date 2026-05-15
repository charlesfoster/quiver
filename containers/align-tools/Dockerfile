# containers/align-tools/Dockerfile
#
# Single image providing:
#   minimap2 2.28  — ONT alignment
#   samtools 1.21  — BAM sort / index / stats
#   bcftools 1.21  — consensus calling (APPLY_CONSENSUS)
#   python 3.11    — inline scripts in APPLY_CONSENSUS
#
# Used by: INDEX_PANEL, HOST_DEPLETE_MINIMAP2, MINIMAP2_ROUND1, MINIMAP2_ROUND2,
#          MINIMAP2_CONSENSUS_MAP, APPLY_CONSENSUS (build_consensus_fasta).
#
# Build target: linux/amd64 (M-series Macs run amd64 containers via Rosetta).
#
# Local build:
#   docker build -t hcv-quasi/align-tools:1.0 containers/align-tools/
#
# Production: pin by digest before deploying to shared infrastructure.

FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="hcv-quasi align-tools"
LABEL org.opencontainers.image.description="minimap2 2.28 + samtools 1.21 + bcftools 1.21 + python 3.11"

USER root

RUN micromamba install -y -n base \
        -c bioconda \
        -c conda-forge \
        -c defaults \
        minimap2=2.28 \
        samtools=1.21 \
        bcftools=1.21 \
        htslib=1.21 \
        python=3.11 \
    && micromamba clean -afy

ENV PATH="/opt/conda/bin:${PATH}"

ENTRYPOINT []
CMD ["/bin/bash"]
