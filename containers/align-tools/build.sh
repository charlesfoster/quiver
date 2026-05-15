#!/usr/bin/env bash
# containers/align-tools/build.sh
#
# Build the align-tools container image from the Dockerfile in this directory.
#
# Usage:
#   bash containers/align-tools/build.sh
#
# The resulting image is tagged hcv-quasi/align-tools:1.0 — the exact name
# referenced in modules/local/{index_panel,host_deplete,minimap2_round1,
# minimap2_round2,minimap2_consensus_map,build_consensus_fasta}.nf.
#
# Build target: linux/amd64
#   On Apple Silicon (M-series) Docker uses Rosetta to emulate amd64.
#
# Estimated build time: 3–8 minutes (conda solve + package download dominates).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

docker build \
    --platform linux/amd64 \
    -t hcv-quasi/align-tools:1.0 \
    "${SCRIPT_DIR}"

echo "Built hcv-quasi/align-tools:1.0"
echo ""
echo "Verify with:"
echo "  docker run --platform linux/amd64 hcv-quasi/align-tools:1.0 minimap2 --version"
echo "  docker run --platform linux/amd64 hcv-quasi/align-tools:1.0 samtools --version"
echo "  docker run --platform linux/amd64 hcv-quasi/align-tools:1.0 bcftools --version"
