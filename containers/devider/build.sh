#!/usr/bin/env bash
# containers/devider/build.sh
#
# Build the DEVIDER container image from the Dockerfile in this directory.
#
# Usage:
#   bash containers/devider/build.sh
#
# The resulting image is tagged hcv-quasi/devider:0.0.1 — the exact name
# referenced in modules/local/devider.nf.
#
# Build target: linux/amd64
#   On Apple Silicon (M-series) Docker uses Rosetta to emulate amd64;
#   this is acceptable for development.  For production, build on a native
#   amd64 host or use docker buildx with --platform linux/amd64.
#
# The build fetches DEVIDER v0.0.1 from GitHub via cargo install.
# Ensure network access is available during the build.
#
# Estimated build time: 5–15 minutes (Rust compile time dominates).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

docker build \
    --platform linux/amd64 \
    -t hcv-quasi/devider:0.0.1 \
    "${SCRIPT_DIR}"

echo "Built hcv-quasi/devider:0.0.1"
echo ""
echo "Verify with:"
echo "  docker run --platform linux/amd64 hcv-quasi/devider:0.0.1 --help"
