#!/usr/bin/env bash

# -------------------------------------------------------------------------------------
# DATASET: XeniumIS_Mouse_Brain_TgCRND8 (10x Xenium In Situ)
# -------------------------------------------------------------------------------------
#
# ## Purpose
# This script downloads the 10x Genomics Xenium In Situ Gene Expression data for a
# mouse brain section. This dataset is from a transgenic mouse model (TgCRND8) that
# serves as a model for Alzheimer's disease.
#
# ## Sample Details
# * **Tissue Type:** Formalin-Fixed Paraffin-Embedded (FFPE) Mouse Brain.
# * **Assay:** Xenium In Situ (targeted gene panel).
# * **Key Genes:** This assay measures a panel of key genes relevant to the mouse
#   brain, including those involved in neurodegeneration.
#
# ## Source & Download Files
# Source Name: Xenium_V1_FFPE_TgCRND8_17_9_months
# Direct Link: https://cf.10xgenomics.com/samples/xenium/1.4.0/Xenium_V1_FFPE_TgCRND8_17_9_months/Xenium_V1_FFPE_TgCRND8_17_9_months_web_summary.html
#
# -------------------------------------------------------------------------------------

set -euo pipefail

# --- CONFIGURATION ---
BASE_URL="https://cf.10xgenomics.com/samples/xenium/1.4.0/Xenium_V1_FFPE_TgCRND8_17_9_months"

# Files to download: A single ZIP archive
FILES=(
    "Xenium_V1_FFPE_TgCRND8_17_9_months_outs.zip"
)

# --- DOWNLOAD STAGE ---
echo "Starting download"
# Check for the existence of the output directory to skip if already present
OUTPUT_DIR="${FILES[0]%.zip}"
if [ -d "$OUTPUT_DIR" ]; then
    echo "'$OUTPUT_DIR' already exists. Skipping download and extraction."
else
    for file in "${FILES[@]}"; do
        echo "  Downloading $file..."
        wget -nc "$BASE_URL/$file"
    done
    echo "Data download complete."

    # --- EXTRACTION STAGE ---
    echo "Extracting archives..."
    for zip in *.zip; do
        unzip -q "$zip" -d .
        # Clean up the large archive after extraction
        rm "$zip"
    done
    echo "Extraction complete."
fi
