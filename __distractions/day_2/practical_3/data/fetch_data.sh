#! /usr/bin/env bash

# Exit immediately if a command exits with a non-zero status,
# Use unbound variables as errors, and fail in a pipeline if any command fails.
set -euo pipefail

# --- Configuration ---

DATA_DIR="Zhuang-ABCA-1_Dataset"

# Base URL for the Allen Brain Cell Atlas S3 bucket
BASE_URL="https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com"

# Array of file paths relative to the base URL
FILES=(
    "expression_matrices/Zhuang-ABCA-1/20230830/Zhuang-ABCA-1-log2.h5ad"
    "expression_matrices/Zhuang-ABCA-1/20230830/Zhuang-ABCA-1-raw.h5ad"
    "metadata/Zhuang-ABCA-1/20241115/cell_metadata.csv"
    "metadata/Zhuang-ABCA-1/20241115/gene.csv"
    "metadata/Zhuang-ABCA-1/20241115/views/cell_metadata_with_cluster_annotation.csv"
    "metadata/Zhuang-ABCA-1-CCF/20230830/ccf_coordinates.csv"
)

# --- Main Download Logic ---

if [ -d "$DATA_DIR" ]; then
    echo "${DATA_DIR}: Data directory already exists. Checking for missing files."
else
    echo "Downloading Zhuang-ABCA-1 data."
    mkdir -p "$DATA_DIR"
fi

for relative_path in "${FILES[@]}"; do
    # Construct the full URL
    FULL_URL="$BASE_URL/$relative_path"

    # Determine the local file path by replacing the slashes with underscores for simplicity,
    # or use the full path structure if preferred (see note below).
    # We will use the full path structure for better organization.

    # Extract the directory path from the relative_path (e.g., 'expression_matrices/Zhuang-ABCA-1/20230830')
    FILE_DIR=$(dirname "$relative_path")
    # Extract the filename (e.g., 'Zhuang-ABCA-1-log2.h5ad')
    FILENAME=$(basename "$relative_path")

    # Create the full local path for the file, maintaining the directory structure
    LOCAL_FILE_DIR="$DATA_DIR/$FILE_DIR"
    LOCAL_FILE_PATH="$LOCAL_FILE_DIR/$FILENAME"

    # Ensure the local directory structure exists
    mkdir -p "$LOCAL_FILE_DIR"

    echo "Downloading: ${FILENAME}..."

    # -nc: no clobber, skip download if file exists
    # -O: output document to specific file (using the full path)
    wget -nc -O "$LOCAL_FILE_PATH" "$FULL_URL"
done

echo "Data download complete. Files are saved in the '$DATA_DIR' directory, maintaining their original path structure."