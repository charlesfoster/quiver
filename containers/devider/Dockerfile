# containers/devider/Dockerfile
#
# Multi-stage build: compile DEVIDER from source (Rust), copy only the binary
# to a slim Debian runtime image.
#
# Build target: linux/amd64 (primary).
# M5 Max Mac runs amd64 containers via Rosetta — acceptable for development.
#
# Production note: pin base images by digest before deploying to shared
# infrastructure.  Replace the tagged FROM lines with:
#   FROM rust:1.78-slim@sha256:<digest> AS builder
#   FROM debian:bookworm-slim@sha256:<digest>
# Digest pinning is omitted here to keep the development workflow simple;
# tags are sufficient for local testing.
#
# DEVIDER version: v0.0.1 (only tagged release as of 2026-05-13).
# Source: https://github.com/bluenote-1577/devider
#
# If a newer tag exists at implementation time, update DEVIDER_TAG, rebuild,
# and verify CLI flags against `devider --help` before changing modules/local/devider.nf.

# ---------------------------------------------------------------------------
# Stage 1: build
# Compile DEVIDER from source using the official Rust toolchain image.
# ---------------------------------------------------------------------------
FROM rust:1.78-slim AS builder

# Build dependencies:
#   git        — cargo install fetches the repo via git
#   cmake      — some crate build scripts require cmake
#   pkg-config — needed to locate system libraries
#   libssl-dev — openssl-sys crate links against libssl
RUN apt-get update && apt-get install -y \
    git \
    cmake \
    pkg-config \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install DEVIDER v0.0.1 from the tagged release on GitHub.
# --root /usr/local places the binary at /usr/local/bin/devider.
RUN cargo install \
    --git https://github.com/bluenote-1577/devider \
    --tag v0.0.1 \
    --root /usr/local

# ---------------------------------------------------------------------------
# Stage 2: runtime
# Slim Debian image containing only the compiled binary and its runtime libs.
# ---------------------------------------------------------------------------
FROM debian:bookworm-slim

# Runtime dependencies:
#   libssl3      — dynamic link target for the openssl-sys crate
#   ca-certificates — required if devider makes any HTTPS calls at runtime
RUN apt-get update && apt-get install -y \
    libssl3 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/bin/devider /usr/local/bin/devider

ENTRYPOINT ["/usr/local/bin/devider"]
