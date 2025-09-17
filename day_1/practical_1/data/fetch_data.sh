#! /usr/bin/env bash

set -euo pipefail

XENIUM_DATA="Xenium_V1_FFPE_TgCRND8_17_9_months_outs"
if [ -d "$XENIUM_DATA" ]; then
    echo "${XENIUM_DATA}: Data directory already exists. Skipping extraction."
else
    echo "Downloading data for Practical 1"
    mkdir -p "$XENIUM_DATA"
    BASE_URL="https://cf.10xgenomics.com/samples/xenium/1.4.0/Xenium_V1_FFPE_TgCRND8_17_9_months"

    FILES=(
        "Xenium_V1_FFPE_TgCRND8_17_9_months_outs.zip"
    )

    for file in "${FILES[@]}"; do
        wget -nc -O "$XENIUM_DATA/$file" "$BASE_URL/$file"
    done
    echo "Data download complete."

    echo "Extracting downloaded data..."
    unzip "$XENIUM_DATA/Xenium_V1_FFPE_TgCRND8_17_9_months_outs.zip" -d "$XENIUM_DATA"
    echo "Extraction complete."
fi
