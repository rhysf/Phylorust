# Use the official Rust image to build
FROM rust:1.81-slim AS builder

# Install required system deps (FastTree + R + libraries)
RUN apt-get update && apt-get install -y \
    build-essential \
    fasttree \
    r-base \
    r-cran-ggplot2 \
    r-cran-dplyr \
    r-cran-readr \
    && rm -rf /var/lib/apt/lists/*

# Create app directory
WORKDIR /app

# Copy sources
COPY . .

# Build the binary in release mode
RUN cargo build --release

# Final minimal image
FROM debian:bookworm-slim

# Install runtime dependencies (FastTree, R)
RUN apt-get update && apt-get install -y \
    fasttree \
    r-base \
    r-cran-ggplot2 \
    r-cran-dplyr \
    r-cran-readr \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/local/bin

# Copy Phylorust binary only
COPY --from=builder /app/target/release/phylorust .

# Default entrypoint
ENTRYPOINT ["phylorust"]
