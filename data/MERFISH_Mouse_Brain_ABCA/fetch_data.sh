#!/usr/bin/env bash

# -------------------------------------------------------------------------------------
# DATASET: MERFISH_Mouse_Brain_ABCA (Allen Brain Cell Atlas)
# -------------------------------------------------------------------------------------
#
# ## Purpose
# This script downloads a primary dataset from the Allen Brain Cell Atlas (ABCA).
# The data uses MERFISH spatial transcriptomics to map gene expression at single-cell
# resolution and is aligned to the 3D Common Coordinate Framework (CCF) for the
# adult mouse brain.
#
# ## Sample Details
# * **Tissue Type:** Adult Mouse Brain (sectioned).
# * **Assay:** MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization).
# * **Resolution:** Single-cell spatial resolution.
# * **Alignment:** Aligned to the Allen Mouse Brain CCF (v3).
#
# ## Source & Download Files
# Source Name: Zhuang-ABCA-1_Dataset (Hosted on Allen S3)
# Base URL: https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com
#
# -------------------------------------------------------------------------------------

set -euo pipefail

# --- CONFIGURATION ---
BASE_URL="https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com"

FILES=(
    "expression_matrices/Zhuang-ABCA-1/20230830/Zhuang-ABCA-1-log2.h5ad"
    "expression_matrices/Zhuang-ABCA-1/20230830/Zhuang-ABCA-1-raw.h5ad"
    "metadata/Zhuang-ABCA-1/20241115/cell_metadata.csv"
    "metadata/Zhuang-ABCA-1/20241115/gene.csv"
    "metadata/Zhuang-ABCA-1/20241115/views/cell_metadata_with_cluster_annotation.csv"
    "metadata/Zhuang-ABCA-1-CCF/20230830/ccf_coordinates.csv"
)

# --- DOWNLOAD STAGE ---
echo "Starting download"
if [ -f "${FILES[0]}" ]; then
    echo "${FILES[0]} exists. Skipping download and extraction."
else
    for relative_path in "${FILES[@]}"; do
        filename=$(basename "$relative_path")
        filedir=$(dirname "$relative_path")
        echo "  Downloading $filename..."
        # Ensure the directory exists
        mkdir -p "$filedir"
        wget -nc "$BASE_URL/$relative_path" -O "$relative_path"
    done
    echo "Data download complete."
fi
