#! /usr/bin/env bash

set -euo pipefail

VISIUM_DATA="visium_2.1.0_2_io_subset"
if [ -d "$VISIUM_DATA" ]; then
    echo "${VISIUM_DATA}: Data directory already exists. Skipping extraction."
else
    echo "Downloading data for Practical 0"
    mkdir -p "$VISIUM_DATA"
    wget -nc -O "$VISIUM_DATA/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_filtered_feature_bc_matrix.h5" "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_filtered_feature_bc_matrix.h5"
    wget -nc -O "$VISIUM_DATA/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_spatial.tar.gz" "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_spatial.tar.gz"
    echo "Data download complete."

    echo "Extracting downloaded data..."
    tar -xvzf "$VISIUM_DATA/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_spatial.tar.gz" -C "$VISIUM_DATA"
    echo "Extraction complete."
fi

VISIUMHD_DATA="visium_hd_3.0.0_io_subset"
if [ -d "$VISIUMHD_DATA" ]; then
    echo "${VISIUMHD_DATA}: Data directory already exists. Skipping extraction."
else
    echo "Downloading data for Practical 0"
    mkdir -p "$VISIUMHD_DATA"
    wget -nc -O "$VISIUMHD_DATA/Visium_HD_Mouse_Small_Intestine_feature_slice.h5" "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_feature_slice.h5"
    wget -nc -O "$VISIUMHD_DATA/Visium_HD_Mouse_Small_Intestine_molecule_info.h5" "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_molecule_info.h5"
    wget -nc -O "$VISIUMHD_DATA/Visium_HD_Mouse_Small_Intestine_spatial.tar.gz" "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_spatial.tar.gz"
    wget -nc -O "$VISIUMHD_DATA/Visium_HD_Mouse_Small_Intestine_binned_outputs.tar.gz" "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Small_Intestine/Visium_HD_Mouse_Small_Intestine_binned_outputs.tar.gz"
    echo "Data download complete."

    echo "Extracting downloaded data..."
    tar -xvzf "$VISIUMHD_DATA/Visium_HD_Mouse_Small_Intestine_spatial.tar.gz" -C "$VISIUMHD_DATA"
    tar -xvzf "$VISIUMHD_DATA/Visium_HD_Mouse_Small_Intestine_binned_outputs.tar.gz" -C "$VISIUMHD_DATA"
    echo "Extraction complete."
fi
