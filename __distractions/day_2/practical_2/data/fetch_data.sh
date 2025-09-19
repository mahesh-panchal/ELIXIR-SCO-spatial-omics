#! /usr/bin/env bash

set -euo pipefail

VISIUM_DATA="V1_Adult_Mouse_Brain_data"
if [ -d "$VISIUM_DATA" ]; then
    echo "${VISIUM_DATA}: Data directory already exists. Skipping download."
else
    echo "Downloading V1_Adult_Mouse_Brain data."
    mkdir -p "$VISIUM_DATA"
    BASE_URL="https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain"

    FILES=(
        "V1_Adult_Mouse_Brain_molecule_info.h5"
        "V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5"
        "V1_Adult_Mouse_Brain_raw_feature_bc_matrix.h5"
        "V1_Adult_Mouse_Brain_analysis.tar.gz"
        "V1_Adult_Mouse_Brain_spatial.tar.gz"
        "V1_Adult_Mouse_Brain_metrics_summary.csv"
        "V1_Adult_Mouse_Brain_web_summary.html"
        "V1_Adult_Mouse_Brain_cloupe.cloupe"
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