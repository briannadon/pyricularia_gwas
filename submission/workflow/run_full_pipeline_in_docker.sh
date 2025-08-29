#!/bin/bash

# Script to run the full RNA-seq pipeline inside the bioinf-bench Docker container

set -e  # Exit on error

# Configuration
CONTAINER_NAME="bioinf-bench"
WORKDIR="/workspace/"
LOCAL_DIR="../"

# Function to log messages with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    log "ERROR: Docker is not running. Please start Docker and try again."
    exit 1
fi

# Check if the container exists
if ! docker inspect "$CONTAINER_NAME" > /dev/null 2>&1; then
    log "ERROR: Container '$CONTAINER_NAME' not found. Please build it first."
    log "You can build it with: docker build -t bioinf-bench ."
    exit 1
fi

log "Starting full pipeline in Docker container '$CONTAINER_NAME'..."

# Run the pipeline inside the container
docker run --rm -it \
    -v "$LOCAL_DIR:/workspace/" \
    -w "/workspace" \
    "$CONTAINER_NAME" \
    bash -c "
    # Source conda initialization
    source /opt/conda/etc/profile.d/conda.sh
    source /opt/conda/etc/profile.d/mamba.sh
    
    echo '=== Starting Full Pipeline ==='
    echo 'Current directory: $(pwd)'
    echo 'Listing contents:'
    ls -la
    
    echo 'Changing to workflow directory...'
    cd workflow
    
    echo 'Running: ./pipeline.sh'
    
    # Execute the full pipeline
    bash ./pipeline.sh
    "

log "Full pipeline completed in Docker container!"