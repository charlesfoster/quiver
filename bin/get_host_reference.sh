#!/usr/bin/env bash
# get_host_reference.sh — Download and index GRCh38 no-alt reference for host depletion.
#
# Usage:
#   get_host_reference.sh <output_prefix>
#
# Arguments:
#   output_prefix — Path prefix for the output files (without extension).
#                   The script writes:
#                     <output_prefix>.fna.gz   — downloaded reference (FASTA, gzip)
#                     <output_prefix>.mmi       — minimap2 index for map-ont
#
# Behaviour:
#   - If <output_prefix>.mmi already exists, the script exits immediately (skip).
#   - If <output_prefix>.fna.gz exists but is not a valid gzip, it is deleted and
#     re-downloaded.
#   - Uses wget if available, otherwise curl.
#
# Requirements:
#   minimap2 (on PATH), wget or curl, gzip (for validation)
#
# Exit codes:
#   0 — success (index written or already existed)
#   1 — failure (download or indexing error)

set -euo pipefail

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
GRCH38_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
THREADS="${NTHREADS:-4}"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
if [ "$#" -ne 1 ]; then
    echo "Usage: $(basename "$0") <output_prefix>" >&2
    echo "Example: $(basename "$0") /data/references/grch38_no_alt" >&2
    exit 1
fi

OUTPUT_PREFIX="$1"
FASTA_GZ="${OUTPUT_PREFIX}.fna.gz"
MMI="${OUTPUT_PREFIX}.mmi"

# ---------------------------------------------------------------------------
# Helper: download a URL to a target path
# ---------------------------------------------------------------------------
download() {
    local url="$1"
    local dest="$2"
    echo "[get_host_reference] Downloading: ${url}" >&2
    if command -v wget &>/dev/null; then
        wget --no-verbose --show-progress -O "${dest}" "${url}"
    elif command -v curl &>/dev/null; then
        curl --progress-bar -L -o "${dest}" "${url}"
    else
        echo "ERROR: Neither wget nor curl is available. Cannot download reference." >&2
        exit 1
    fi
}

# ---------------------------------------------------------------------------
# Skip if index already exists
# ---------------------------------------------------------------------------
if [ -f "${MMI}" ]; then
    echo "[get_host_reference] Index already exists: ${MMI} — skipping download/indexing." >&2
    exit 0
fi

# ---------------------------------------------------------------------------
# Ensure output directory exists
# ---------------------------------------------------------------------------
OUTPUT_DIR="$(dirname "${OUTPUT_PREFIX}")"
mkdir -p "${OUTPUT_DIR}"

# ---------------------------------------------------------------------------
# Download FASTA if not present or corrupt
# ---------------------------------------------------------------------------
if [ -f "${FASTA_GZ}" ]; then
    echo "[get_host_reference] Checking existing FASTA: ${FASTA_GZ}" >&2
    if ! gzip -t "${FASTA_GZ}" &>/dev/null; then
        echo "[get_host_reference] FASTA is corrupt — removing and re-downloading." >&2
        rm -f "${FASTA_GZ}"
    else
        echo "[get_host_reference] Existing FASTA is valid." >&2
    fi
fi

if [ ! -f "${FASTA_GZ}" ]; then
    download "${GRCH38_URL}" "${FASTA_GZ}"
    # Validate download
    if ! gzip -t "${FASTA_GZ}" &>/dev/null; then
        echo "ERROR: Downloaded file is not a valid gzip: ${FASTA_GZ}" >&2
        rm -f "${FASTA_GZ}"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Build minimap2 index
# ---------------------------------------------------------------------------
echo "[get_host_reference] Building minimap2 index: ${MMI}" >&2
echo "[get_host_reference] This may take 10–20 minutes and requires ~8 GB RAM." >&2

minimap2 \
    -x map-ont \
    -t "${THREADS}" \
    -d "${MMI}" \
    "${FASTA_GZ}"

if [ ! -f "${MMI}" ]; then
    echo "ERROR: minimap2 index was not created: ${MMI}" >&2
    exit 1
fi

echo "[get_host_reference] Done. Index written to: ${MMI}" >&2
exit 0
