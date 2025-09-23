#!/usr/bin/env bash

# -------------------------------------------------------------------------------------
# DATASET: Visium_Mouse_Brain_Coronal (10x Visium v1)
# -------------------------------------------------------------------------------------
#
# ## Purpose
# This script downloads the Visium v1 Spatial Gene Expression data for a Fresh Frozen
# Adult Mouse Brain Coronal Section.
#
# ## Sample Details
# * **Tissue Type:** Fresh Frozen Adult Mouse Brain (Coronal Section).
# * **Assay:** Visium Spatial Gene Expression (v1 - 55 µm spots).
# * **Staining:** Hematoxylin & Eosin (H&E).
# * **Organism Details:** Male, C57BL/6 strain, >8 weeks old.
#
# ## Source & Download Files
# Source Name: V1_Adult_Mouse_Brain
# Direct Link: https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_web_summary.html
#
# -------------------------------------------------------------------------------------

set -euo pipefail

# --- CONFIGURATION ---
BASE_URL="https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain"

FILES=(
    "V1_Adult_Mouse_Brain_molecule_info.h5"
    "V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5"
    "V1_Adult_Mouse_Brain_raw_feature_bc_matrix.h5"
    "V1_Adult_Mouse_Brain_analysis.tar.gz"
    "V1_Adult_Mouse_Brain_spatial.tar.gz"
    "V1_Adult_Mouse_Brain_metrics_summary.csv"
    "V1_Adult_Mouse_Brain_cloupe.cloupe"
)

# --- DOWNLOAD STAGE ---
echo "Starting download"
# Check for a marker file to skip if already downloaded.
if [ -f "${FILES[0]}" ]; then
    echo "${FILES[0]} exists. Skipping download and extraction."
else
    for file in "${FILES[@]}"; do
        echo "  Downloading $file..."
        wget -nc "$BASE_URL/$file"
    done
    echo "Data download complete."

    # --- EXTRACTION STAGE ---
    echo "Extracting tarballs..."
    for tar in *.tar.gz; do
        tar xvzf "$tar" -C .
        # Clean up the large tarball after extraction
        rm "$tar"
    done
    echo "Extraction complete."
fi
