# Stage 1: Build the application in a Rust container
FROM rust:1.85-slim as builder

# Install runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        libssl-dev pkg-config \
        curl 
# Create app directory
RUN mkdir -p /app
WORKDIR /app
COPY . .
# Build the application
RUN cargo build --release
# Stage 2: Create a minimal runtime image
FROM debian:bookworm-slim
# Install runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        pkg-config libssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/pathogenic_variant_finder /usr/bin/pathogenic
