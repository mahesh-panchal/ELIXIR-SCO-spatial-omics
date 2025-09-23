#! /usr/bin/env bash

set -euo pipefail

VISIUM_DATA="visium_2.1.0_2_io_subset"
if [ -d "$VISIUM_DATA" ]; then
    echo "${VISIUM_DATA}: Data directory already exists. Skipping extraction."
else
    echo "Downloading data for Practical 0"
    mkdir -p "$VISIUM_DATA"
    BASE_URL="https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma"

    FILES=(
        "CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_filtered_feature_bc_matrix.h5"
        "CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_spatial.tar.gz"
    )

    for file in "${FILES[@]}"; do
        wget -nc -O "$VISIUM_DATA/$file" "$BASE_URL/$file"
    done
    echo "Data download complete."

    echo "Extracting tarballs..."
    for tar in "$VISIUM_DATA"/*.tar.gz; do
        tar xvzf "$tar" -C "$VISIUM_DATA"
    done
    echo "Extraction complete."
fi

VISIUMHD_DATA="visium_hd_3.0.0_io_subset"
if [ -d "$VISIUMHD_DATA" ]; then
    echo "${VISIUMHD_DATA}: Data directory already exists. Skipping extraction."
else
    echo "Downloading data for Practical 0"
    mkdir -p "$VISIUMHD_DATA"
    BASE_URL="https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine"
    FILES=(
        "Visium_HD_Mouse_Small_Intestine_feature_slice.h5"
        "Visium_HD_Mouse_Small_Intestine_molecule_info.h5"
        "Visium_HD_Mouse_Small_Intestine_spatial.tar.gz"
        "Visium_HD_Mouse_Small_Intestine_binned_outputs.tar.gz"
    )

    for file in "${FILES[@]}"; do
        wget -nc -O "$VISIUMHD_DATA/$file" "$BASE_URL/$file"
    done
    echo "Data download complete."

    echo "Extracting tarballs..."
    for tar in "$VISIUMHD_DATA"/*.tar.gz; do
        tar xvzf "$tar" -C "$VISIUMHD_DATA"
    done
    echo "Extraction complete."
fi
