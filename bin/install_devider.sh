#!/usr/bin/env bash
# bin/install_devider.sh
#
# Install DEVIDER from source using cargo.
#
# This script is the conda-profile fallback for DEVIDER — used when Docker/
# Singularity containers are not available (e.g. --profile conda_local).
#
# No bioconda package exists for DEVIDER (as of 2026-05-13); installation
# requires a working Rust toolchain.  Install Rust via https://rustup.rs/
# if cargo is not present.
#
# Usage:
#   bash bin/install_devider.sh [INSTALL_PREFIX]
#
#   INSTALL_PREFIX  Optional.  The --root prefix passed to cargo install.
#                   Defaults to $HOME/.local  (binary lands at ~/.local/bin/devider).
#                   The directory must already exist or cargo will create it.
#
# After installation, ensure $INSTALL_PREFIX/bin is in your PATH.
#
# DEVIDER version: v0.0.1 (only tagged release as of 2026-05-13).
# Source: https://github.com/bluenote-1577/devider
#
# If a newer tag exists at implementation time, update DEVIDER_VERSION below,
# then verify CLI flags against `devider --help` before updating devider.nf.

set -euo pipefail

DEVIDER_VERSION="v0.0.1"
INSTALL_PREFIX="${1:-$HOME/.local}"

# ----------------------------------------------------------------
# Guard: skip if a compatible devider binary is already on PATH.
# ----------------------------------------------------------------
if command -v devider &>/dev/null; then
    installed=$(devider --version 2>&1 | head -1 || echo "unknown version")
    echo "devider is already installed: ${installed}"
    echo "Location: $(command -v devider)"
    exit 0
fi

# ----------------------------------------------------------------
# Prerequisite: cargo must be available.
# ----------------------------------------------------------------
if ! command -v cargo &>/dev/null; then
    echo "ERROR: cargo not found." >&2
    echo "Install the Rust toolchain via https://rustup.rs/ and re-run this script." >&2
    exit 1
fi

echo "Installing DEVIDER ${DEVIDER_VERSION} into ${INSTALL_PREFIX}/bin ..."
echo "(This compiles from source — expect 5–15 minutes depending on hardware.)"

cargo install \
    --git https://github.com/bluenote-1577/devider \
    --tag "${DEVIDER_VERSION}" \
    --root "${INSTALL_PREFIX}"

# ----------------------------------------------------------------
# Verify installation.
# ----------------------------------------------------------------
DEVIDER_BIN="${INSTALL_PREFIX}/bin/devider"

if [[ ! -x "${DEVIDER_BIN}" ]]; then
    echo "ERROR: cargo install completed but ${DEVIDER_BIN} is not executable." >&2
    exit 1
fi

echo ""
echo "DEVIDER installed at ${DEVIDER_BIN}"
"${DEVIDER_BIN}" --version 2>&1 | head -1

echo ""
echo "If ${INSTALL_PREFIX}/bin is not already on your PATH, add:"
echo "  export PATH=\"${INSTALL_PREFIX}/bin:\$PATH\""
