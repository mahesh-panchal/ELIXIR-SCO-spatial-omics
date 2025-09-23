#! /usr/bin/env bash

# -------------------------------------------------------------------------------------
# DATASET: SpatialHD_Mouse_Small_Intestine (Visium HD)
# -------------------------------------------------------------------------------------
#
# ## Purpose
# This script downloads the Visium High Definition (HD) Spatial Gene Expression
# data for a Mouse Small Intestine section.
#
# ## Sample Details
# * **Tissue Type:** Formalin-Fixed Paraffin-Embedded (FFPE) Mouse Small Intestine (Swiss roll).
# * **Assay:** Visium HD Spatial Gene Expression (Whole Transcriptome Probe Set v2.0).
# * **Resolution:** High Definition (HD) - Continuous 2x2 µm barcoded squares.
# * **Staining:** H&E stained.
# * **Organism:** Mouse (C57BL/6), 8 weeks old, Male.
#
# ## Source & Download Files
# The raw Space Ranger outputs contain the 2µm, 8µm, and 16µm binned data.
# Direct Link: https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_web_summary.html
#
# -------------------------------------------------------------------------------------

set -euo pipefail

# --- CONFIGURATION ---
BASE_URL="https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine"

# Files to download:
FILES=(
    "Visium_HD_Mouse_Small_Intestine_feature_slice.h5"
    "Visium_HD_Mouse_Small_Intestine_molecule_info.h5"
    "Visium_HD_Mouse_Small_Intestine_spatial.tar.gz"
    "Visium_HD_Mouse_Small_Intestine_binned_outputs.tar.gz"
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
