#!/usr/bin/env bash
# test_integration.sh — Integration test harness for hcv-quasi
#
# Usage:
#   tests/test_integration.sh              # stub-run only (fast, no tools needed)
#   tests/test_integration.sh --full       # full pipeline run (requires all tools)
#   tests/test_integration.sh --full --outdir results/test
#
# Exit codes:
#   0 — all checks passed
#   1 — one or more checks failed

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NXF_VER="${NXF_VER:-24.10.5}"

FULL_RUN=false
OUTDIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --full)
            FULL_RUN=true
            shift
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            echo "Usage: $0 [--full] [--outdir DIR]" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$OUTDIR" ]]; then
    OUTDIR="${PROJECT_DIR}/test_out_$(date +%Y%m%d_%H%M%S)"
fi

echo "=== hcv-quasi integration test ==="
echo "Project : ${PROJECT_DIR}"
echo "Nextflow: NXF_VER=${NXF_VER}"
echo "Outdir  : ${OUTDIR}"
echo "Mode    : $(${FULL_RUN} && echo 'FULL' || echo 'stub-run')"
echo ""

# ── Step 1: Check Nextflow is available ───────────────────────────────────────
echo "[1/4] Checking Nextflow availability..."
if ! command -v nextflow &>/dev/null; then
    echo "FAIL: 'nextflow' not found in PATH." >&2
    echo "      Install Nextflow >= 24.10.5 or set NXF_VER and ensure nextflow is on PATH." >&2
    exit 1
fi
NF_ACTUAL_VER=$(NXF_VER="${NXF_VER}" nextflow -version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo "unknown")
echo "    Nextflow version: ${NF_ACTUAL_VER}  (requested: ${NXF_VER})"

# ── Step 2: Stub-run (always) ─────────────────────────────────────────────────
echo ""
echo "[2/4] Running stub-run to validate workflow structure..."
STUB_LOG="${PROJECT_DIR}/nextflow_stub.log"
if NXF_VER="${NXF_VER}" nextflow run "${PROJECT_DIR}/main.nf" \
        -profile test \
        -stub-run \
        --outdir "${OUTDIR}_stub" \
        -ansi-log false \
        2>&1 | tee "${STUB_LOG}"; then
    echo "    Stub-run: PASS"
else
    echo ""
    echo "FAIL: Stub-run exited non-zero." >&2
    echo "      Last 30 lines of ${STUB_LOG}:" >&2
    tail -30 "${STUB_LOG}" >&2
    exit 1
fi

# ── Step 3: Full run (optional) ───────────────────────────────────────────────
if ${FULL_RUN}; then
    echo ""
    echo "[3/4] Running full pipeline (this may take a few minutes)..."
    FULL_LOG="${PROJECT_DIR}/nextflow_test.log"
    if NXF_VER="${NXF_VER}" nextflow run "${PROJECT_DIR}/main.nf" \
            -profile test \
            --outdir "${OUTDIR}" \
            -ansi-log false \
            2>&1 | tee "${FULL_LOG}"; then
        echo "    Full run: PASS"
    else
        echo ""
        echo "FAIL: Full pipeline run exited non-zero." >&2
        echo "      Last 50 lines of ${FULL_LOG}:" >&2
        tail -50 "${FULL_LOG}" >&2
        exit 1
    fi

    # ── Step 4: Validate outputs ───────────────────────────────────────────────
    echo ""
    echo "[4/4] Validating outputs against expected values..."
    if python3 "${SCRIPT_DIR}/check_outputs.py" \
            "${OUTDIR}" \
            "${PROJECT_DIR}/test_data/expected_outputs"; then
        echo ""
        echo "PASS: Integration test complete."
        exit 0
    else
        echo ""
        echo "FAIL: Output validation failed (see above)." >&2
        exit 1
    fi
else
    echo ""
    echo "[3/4] Skipping full pipeline run (use --full to enable)."
    echo "[4/4] Skipping output validation (requires full run)."
    echo ""
    echo "PASS: Stub-run validation complete."
    echo "      Re-run with --full to execute the complete pipeline and validate outputs."
    exit 0
fi
