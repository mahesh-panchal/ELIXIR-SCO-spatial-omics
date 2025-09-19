#! /usr/bin/env bash

# -------------------------------------------------------------------------------
# DATASET: SpatialGEP_Human_Glioblastoma (10x Genomics Visium)
# -------------------------------------------------------------------------------
# 
# ## Purpose
# This script downloads the Visium CytAssist Spatial Gene and Protein Expression data
# for a Human Glioblastoma multiforme (Grade IV brain cancer) tumor section.
#
# ## Modules
# This dataset is used in:
# * 01_Data_types
#
# ## Sample Details
# * **Tissue Type:** Formalin-Fixed Paraffin-Embedded (FFPE) Human Brain.
# * **Disease:** Glioblastoma multiforme.
# * **Assay:** Visium CytAssist Spatial Gene and Protein Expression (simultaneous RNA and protein).
# * **Staining:** Immunofluorescence (IF) stained.
# * **Patient:** Female, Age 60s.
# * **Key Metrics:** 5,756 Spots detected under tissue; Median 7,629 genes per spot.
#
# ## Source & Download Files
# The raw output files (FASTQs, HDF5, Loupe files) are available from 10x Genomics.
# The primary analysis file is the Space Ranger output archive.
# Direct link: https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_web_summary.html
#
# -------------------------------------------------------------------------------------

set -euo pipefail

# --- CONFIGURATION ---
BASE_URL="https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma"

# Files to download: H5 matrix (RNA/Protein) and the spatial tarball (images/metadata)
FILES=(
    "CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_filtered_feature_bc_matrix.h5"
    "CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_spatial.tar.gz"
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
