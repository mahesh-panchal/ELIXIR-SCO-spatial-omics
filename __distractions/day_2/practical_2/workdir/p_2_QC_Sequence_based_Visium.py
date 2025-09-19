import marimo

__generated_with = "0.15.2"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# **ELIXIR Spatial Transcriptomics Course**""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Author: George Gavriilidis""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Date: 2025-02-22""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Author email: ggeorav@certh.gr""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Practical 2: Quality Control and Visualisation for Sequence-based data analysis

    ### Quiz: Pre-Workshop Check


    1. What is the importance of QC metrics on spatial data?

    2. Evaluating library size - why?

    3. Thresholding across various QC metrics - how to decide?

    4. Cells per spot - why does it matter?

    5. Mitochondria - what do they represent for QC?

    6. Normalisation - why do we need that for downstream analysis?

    7. Spotting Highly Variable Genes - does it matter and how QC is important?

    8. What are some advanced visualisations that `SpatialData` can facilitate for sequence-based analysis?

    Take a moment to reflect on these questions before moving on to the hands-on exercises.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Dependencies""")
    return


@app.cell
def _():
    import spatialdata as sd
    import spatialdata_plot as sdp
    from spatialdata_io import visium, visium_hd, xenium
    import matplotlib.pyplot as plt  # for multi-panel plots later
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from scipy.sparse import issparse

    sc.settings.verbosity = 3
    return issparse, np, pd, plt, sc, sns, visium_hd


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# 2.0 Load Visium dataset""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""The dataset used here consists of a Visium slide of a coronal section of the mouse brain. The original dataset is publicly available at the 10x Genomics dataset portal <https://support.10xgenomics.com/spatial-gene-expression/datasets>_""")
    return


@app.cell
def _(sc):
    # Define the path to the organized Visium data directory
    visium_data_path = "day_2/practical_2/data/V1_Adult_Mouse_Brain_data"

    # Read the Visium data into an AnnData object
    adata = sc.read_visium(
        visium_data_path,
        count_file="V1_Adult_Mouse_Brain_raw_feature_bc_matrix.h5",
    )

    # Display basic information about the AnnData object
    print(adata)
    return (adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# 2.1 Spot-level Quality Control""")
    return


@app.cell
def _():
    # An R implementation of these analytical steps can pe found in the following URL: https://bookdown.org/sjcockell/ismb-tutorial-2023/
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Introduction to Spot-Level Quality Control (sQC)

    Quality control (QC) is essential for analyzing high-throughput molecular biology data. Removing noise and low-quality data from complex datasets enhances the reliability of downstream analyses. Spatial transcriptomics (STx) is no exception, with QC typically performed at two main levels: **spot-level** and **gene-level**. This document focuses on **spot-level QC**.

    ### Spot-Level Quality Control (sQC)

    Spot-level QC eliminates low-quality spots prior to analysis. These low-quality spots can result from issues during library preparation or experimental procedures, such as:
    - **High percentage of dead cells** due to cell damage during library preparation.
    - **Low mRNA capture efficiency** caused by ineffective reverse transcription or PCR amplification.

    Keeping low-quality spots can introduce noise and compromise the accuracy of downstream analyses.

    ### Characteristics of Low-Quality Spots

    Low-quality spots are identified using criteria similar to those used in cell-level QC for single-cell RNA sequencing (scRNA-seq) data:
    - **Library size**: Total UMI counts per spot may vary due to sequencing or the number of cells in the spot.
    - **Number of expressed genes**: The number of genes with non-zero UMI counts per spot.
    - **Proportion of mitochondrial reads**: A high proportion indicates potential cell damage.

    #### Key Indicators:
    1. **Low Library Size / Expressed Features**:
       - Reflects poor mRNA capture rates due to cell damage or low reaction efficiency.
    2. **High Mitochondrial Proportion**:
       - Indicates cell damage, e.g., partial lysis, leading to cytoplasmic mRNA leakage and a higher concentration of relatively protected mitochondrial mRNAs.
    3. **Unusually High Cell Numbers per Spot**:
       - May indicate problems during cell segmentation.

    ### Relevance of scRNA-seq QC Metrics in STx

    The use of scRNA-seq QC metrics in STx data is based on the similarity between datasets. When spatial information is ignored, each spot can be treated as a single cell. However, STx spots differ from scRNA-seq cells because they may contain zero, one, or multiple cells, leading to distinct distributions for high-quality spots compared to high-quality cells in scRNA-seq.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## 2.1.1 Plot tissue map""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Initially, we will try to get a sense of the data by plotting High-Resolution Image with Spatial Coordinates Overlay""")
    return


@app.cell
def _(adata, sc):
    # Plot spatial coordinates without annotations
    sc.pl.spatial(adata, color=None, size=1.5)
    return


@app.cell
def _(adata, sc):
    # Overlay spatial coordinates on the tissue image
    sc.pl.spatial(
        adata,
        color=None,  # No color (default spots)
        size=1.5,  # Adjust spot size
        img_key="hires",  # Use the high-resolution image
        alpha_img=0.6,  # Adjust image transparency
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""At present, the dataset contains both on- and off-tissue spots - we plotted these in the previous practical. For any future analysis though we are only interested in the on-tissue spots. Therefore, before we run any calculations we want to remove the off-tissue spots""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## 2.1.2 Calculating QC metrics""")
    return


@app.cell
def _(adata):
    print(f"Dataset dimensions before filtering: {adata.shape}")
    adata_1 = adata[adata.obs["in_tissue"] == 1]
    print(f"Dataset dimensions after filtering: {adata_1.shape}")
    return (adata_1,)


@app.cell
def _(adata_1):
    print("First 20 gene names in the dataset:")
    print(adata_1.var_names[:20].tolist())
    mito_genes_alt = adata_1.var_names.str.contains(
        "(^MT-)|(^mt-)|(^mt\\-)", regex=True
    )
    print(f"Number of mitochondrial genes detected: {mito_genes_alt.sum()}")
    mito_genes_list = adata_1.var_names[mito_genes_alt]
    print("Mitochondrial genes detected:")
    print(mito_genes_list.tolist())
    return


@app.cell
def _(adata_1, sc):
    adata_1.var["mt"] = adata_1.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
        adata_1, qc_vars=["mt"], percent_top=None, inplace=True
    )
    print(
        adata_1.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].head()
    )
    return


@app.cell
def _(adata_1, plt, sns):
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    sns.histplot(adata_1.obs["total_counts"], kde=False, bins=50, ax=axs[0])
    axs[0].set_title("Total Counts")
    axs[0].set_xlabel("Counts")
    axs[0].set_ylabel("Frequency")
    sns.histplot(adata_1.obs["n_genes_by_counts"], kde=False, bins=50, ax=axs[1])
    axs[1].set_title("Number of Genes by Counts")
    axs[1].set_xlabel("Number of Genes")
    axs[1].set_ylabel("Frequency")
    sns.histplot(adata_1.obs["pct_counts_mt"], kde=False, bins=50, ax=axs[2])
    axs[2].set_title("Mitochondrial Content")
    axs[2].set_xlabel("Mitochondrial %")
    axs[2].set_ylabel("Frequency")
    plt.tight_layout()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""This is an overview of what is happening to our dataset let's focus on each metric below.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## 2.1.3 Library size threshold plot""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We can plot a histogram of the library sizes across spots. The library size is the number of UMI counts in each spot.""")
    return


@app.cell
def _(adata_1, plt, sns):
    plt.figure(figsize=(8, 6))
    sns.histplot(
        adata_1.obs["total_counts"], kde=True, color="gray", label="Density"
    )
    plt.title("Distribution of Library Sizes")
    plt.xlabel("Library Size (Total Counts)")
    plt.ylabel("Density")
    plt.legend()
    plt.show()
    return


@app.cell
def _(adata_1, plt, sns):
    plt.figure(figsize=(8, 6))
    sns.histplot(
        adata_1.obs["total_counts"], kde=True, color="gray", label="Density"
    )
    plt.xlim(0, 20000)
    plt.title("Distribution of Library Sizes (up to 2000)")
    plt.xlabel("Library Size (Total Counts)")
    plt.ylabel("Density")
    plt.legend()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Scatter Plot of Library Size vs. Cell Count""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    <div style="border: 1px solid #ffa6a6; padding: 10px; border-radius: 5px;">
    <span style="color: #ff6666; font-size: 20px;"><b>Reflection Point:</b></span> <span style="font-size: 20px;">Library size threshold plot</span>  
    <ul>
        <li>What is density in the yy axis?</li>

    </div>
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""To get a glimpse of how library size changes spatially:""")
    return


@app.cell
def _(adata_1, sc):
    sc.pl.spatial(
        adata_1, img_key="hires", color=["total_counts", "n_genes_by_counts"]
    )
    return


@app.cell
def _(adata_1):
    adata_1.obs["qc_lib_size"] = adata_1.obs["total_counts"] < 10000
    low_library_count = adata_1.obs["qc_lib_size"].sum()
    print(f"Number of spots with low library size: {low_library_count}")
    return


@app.cell
def _(adata_1, sc):
    sc.pl.spatial(
        adata_1,
        color="qc_lib_size",
        size=1.5,
        title="Spatial Distribution of Low-Library-Size Spots",
    )
    return


@app.cell
def _(adata_1, sc):
    adata_1.obs["qc_lib_size_high"] = adata_1.obs["total_counts"] < 20000
    high_threshold_count = adata_1.obs["qc_lib_size_high"].sum()
    print(
        f"Number of spots flagged with stricter threshold (2000): {high_threshold_count}"
    )
    sc.pl.spatial(
        adata_1,
        color="qc_lib_size_high",
        size=1.5,
        title="Spatial Distribution of High-Threshold Low-Library-Size Spots",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We need to keep in mind here that the threshold is, to an extent, arbitrary. It is therefore important to look at the number of spots that are left out of the dataset by this choice of cut-off value, and also have a look at their putative spatial patterns. If we filtered out spots with biological relevance, then we should observe some patterns on the tissue map that correlate with some of the known biological structures of the tissue. If we do observe such a phenomenon, we have probably set our threshold too high (i.e. not permissive enough).""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    <div style="border: 1px solid #ffa6a6; padding: 10px; border-radius: 5px;">
    <span style="color: #ff6666; font-size: 20px;"><b>Reflection Point:</b></span> <span style="font-size: 20px;">Library size threshold plot</span>  
    <ul>
        <li>Do you see any spatial patterning while we are thresholding?</li>

    </div>
    """
    )
    return


@app.cell
def _(adata_1, sc):
    adata_1.obs["qc_lib_size_high"] = adata_1.obs["total_counts"] < 7000
    high_threshold_count_1 = adata_1.obs["qc_lib_size_high"].sum()
    print(
        f"Number of spots flagged with stricter threshold (2000): {high_threshold_count_1}"
    )
    sc.pl.spatial(
        adata_1,
        color="qc_lib_size_high",
        size=1.5,
        title="Spatial Distribution of High-Threshold Low-Library-Size Spots",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## 2.1.4 Number of expressed genes""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""As we did with the library sizes, we can plot a histogram of the number of expressed genes across spots. A gene is “expressed” in a spot if it has at least one count in it.""")
    return


@app.cell
def _(adata_1, plt, sns):
    plt.figure(figsize=(8, 6))
    sns.histplot(adata_1.obs["n_genes_by_counts"], kde=True, color="gray")
    plt.title("Distribution of Expressed Genes per Spot")
    plt.xlabel("Number of Expressed Genes")
    plt.ylabel("Density")
    plt.legend()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Finally, again as before, we apply the chosen threshold to flag spots with (in this case) fewer than 500 expressed genes.""")
    return


@app.cell
def _(adata_1):
    threshold_expressed_genes = 2000
    adata_1.obs["qc_expressed_genes"] = (
        adata_1.obs["n_genes_by_counts"] < threshold_expressed_genes
    )
    low_quality_genes_count = adata_1.obs["qc_expressed_genes"].sum()
    print(
        f"Number of spots with fewer than {threshold_expressed_genes} expressed genes: {low_quality_genes_count}"
    )
    return


@app.cell
def _(adata_1, sc):
    sc.pl.spatial(
        adata_1,
        color="qc_expressed_genes",
        size=1.5,
        title="Spatial Distribution of Spots with Low Expressed Genes",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Optional: Adjust Threshold to Test Over-Enthusiastic Filtering""")
    return


@app.cell
def _(adata_1, sc):
    strict_threshold_expressed_genes = 4000
    adata_1.obs["qc_expressed_genes_strict"] = (
        adata_1.obs["n_genes_by_counts"] < strict_threshold_expressed_genes
    )
    strict_flagged_count = adata_1.obs["qc_expressed_genes_strict"].sum()
    print(
        f"Number of spots flagged with stricter threshold ({strict_threshold_expressed_genes}): {strict_flagged_count}"
    )
    sc.pl.spatial(
        adata_1,
        color="qc_expressed_genes_strict",
        size=1.5,
        title=f"Spatial Distribution of Spots with Fewer than {strict_threshold_expressed_genes} Expressed Genes",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## 2.1.5 Percentage of mitochondrial expression""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    As we briefly touched on at the beginning, a high proportion of mitochondrial reads indicates low cell quality, probably due to cell damage.

    We calculated this data earlier on in this session, and can now investigate the percentage of mitochondrial expression across spots
    """
    )
    return


@app.cell
def _(adata_1, plt, sns):
    if "pct_counts_mt" in adata_1.obs.columns:
        plt.figure(figsize=(8, 6))
        sns.histplot(adata_1.obs["pct_counts_mt"], kde=True, color="gray")
        plt.title("Distribution of Mitochondrial Expression")
        plt.xlabel("Percentage of Mitochondrial Expression")
        plt.ylabel("Density")
        plt.axvline(x=25, color="red", linestyle="--", label="Threshold (25%)")
        plt.legend()
        plt.show()
    else:
        print(
            "The 'pct_counts_mt' column does not exist in adata.obs. Skipping the histogram plot."
        )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""In this instance, a higher percentage of mitochondrial expression is the thing to avoid, so the threshold is an upper bound, rather than the lower bounds we have observed so far. Our suggestion this time is to cut-off at 28%.""")
    return


@app.cell
def _(adata_1):
    threshold_mito_percent = 25
    if "pct_counts_mt" in adata_1.obs.columns:
        adata_1.obs["qc_mito"] = (
            adata_1.obs["pct_counts_mt"] > threshold_mito_percent
        )
        flagged_mito_count = adata_1.obs["qc_mito"].sum()
        print(
            f"Number of spots with mitochondrial expression above {threshold_mito_percent}%: {flagged_mito_count}"
        )
    else:
        print(
            "The 'pct_counts_mt' column does not exist in adata.obs. Skipping thresholding."
        )
    return


@app.cell
def _(adata_1, sc):
    if "qc_mito" in adata_1.obs.columns:
        sc.pl.spatial(
            adata_1,
            color="qc_mito",
            size=1.5,
            title="Spatial Distribution of High Mitochondrial Expression Spots",
        )
    else:
        print(
            "The 'qc_mito' column does not exist in adata.obs. Skipping spatial visualization."
        )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Again, try to illustrate what happens if we set the threshold too low (i.e., 20 0r 25%).""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## 2.1.6 Remove low-quality spots""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""All the steps so far have flagged spots with potential issues - before proceeding with the analysis; we want to remove these spots from our adata object. Since we have calculated different spot-level QC metrics and selected thresholds for each one, we can combine them to identify a set of low-quality spots, and remove them from our spe object in a single step.""")
    return


@app.cell
def _(adata_1, sc):
    metrics_to_check = ["qc_lib_size", "qc_expressed_genes", "qc_mito"]
    discarded_counts = {}
    for metric in metrics_to_check:
        if metric in adata_1.obs.columns:
            discarded_counts[metric] = adata_1.obs[metric].sum()
        else:
            print(f"Warning: '{metric}' column not found in adata.obs.")
    print("Number of discarded spots for each metric:")
    print(discarded_counts)
    if all((metric in adata_1.obs.columns for metric in metrics_to_check)):
        adata_1.obs["discard"] = (
            adata_1.obs["qc_lib_size"]
            | adata_1.obs["qc_expressed_genes"]
            | adata_1.obs["qc_mito"]
        )
        print(f"Total spots marked as discarded: {adata_1.obs['discard'].sum()}")
    else:
        print(
            "Error: One or more required QC metrics are missing. Cannot combine discarded spots."
        )
    if "discard" in adata_1.obs.columns:
        sc.pl.spatial(
            adata_1,
            color="discard",
            size=1.5,
            title="Spatial Distribution of Combined Low-Quality Spots",
        )
    else:
        print("Error: 'discard' column not found in adata.obs. Cannot visualize.")
    if "ground_truth" in adata_1.obs.columns:
        adata_1.obs["qc_NA_spots"] = adata_1.obs["ground_truth"].isna()
        adata_1.obs["discard"] = (
            adata_1.obs["discard"] | adata_1.obs["qc_NA_spots"]
        )
        print(
            f"Number of spots without annotation: {adata_1.obs['qc_NA_spots'].sum()}"
        )
        print(
            f"Total spots marked as discarded (after including NA spots): {adata_1.obs['discard'].sum()}"
        )
        sc.pl.spatial(
            adata_1,
            color="discard",
            size=1.5,
            title="Spatial Distribution of Discarded Spots (Including NA Annotations)",
        )
    else:
        print(
            "No 'ground_truth' column found in adata.obs. Skipping NA spot filtering."
        )
    return


@app.cell
def _(adata_1):
    adata_1
    return


@app.cell
def _(adata_1, sc):
    sc.pp.filter_cells(adata_1, min_counts=7000)
    sc.pp.filter_cells(adata_1, max_counts=80000)
    adata_1
    return


@app.cell
def _(adata_1):
    adata_2 = adata_1[adata_1.obs["pct_counts_mt"] < 25].copy()
    print(f"#cells after MT filter: {adata_2.n_obs}")
    return (adata_2,)


@app.cell
def _(adata_2, sc):
    sc.pp.filter_genes(adata_2, min_cells=10)
    adata_2
    return


@app.cell
def _(adata_2):
    threshold_expressed_genes_1 = 2000
    adata_2.obs["qc_expressed_genes"] = (
        adata_2.obs["n_genes_by_counts"] >= threshold_expressed_genes_1
    )
    high_quality_genes_count = adata_2.obs["qc_expressed_genes"].sum()
    print(
        f"Number of spots with at least {threshold_expressed_genes_1} expressed genes: {high_quality_genes_count}"
    )
    adata_3 = adata_2[adata_2.obs["qc_expressed_genes"], :]
    print(f"Number of spots retained after filtering: {adata_3.shape[0]}")
    adata_3
    return (adata_3,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    <div style="border: 1px solid #ffa6a6; padding: 10px; border-radius: 5px;">
    <span style="color: #ff6666; font-size: 20px;"><b>Reflection Point:</b></span> <span style="font-size: 20px;">Library size threshold plot</span>  
    <ul>
        <li> What about plotting after QC thresholding?</li>

    </div>
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# 2.2 Normalisation of counts""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # *Normalization in Spatial Transcriptomics (STx)*

    Normalization is essential in Spatial Transcriptomics (STx) to account for systematic effects in count data, such as variations in library size (counts/UMIs per spot), which are not biologically meaningful. Proper normalization ensures accurate comparisons of expression levels and reliable downstream analyses.

    ---

    ## *Library Size Scaling*
    - **Objective:** Adjust library sizes across all spots so that their mean becomes 1.
    - **Method:** 
      - Counts are scaled using library size factors.
      - This adjustment corrects for differences in sequencing depth or total counts.
    - **Outcome:** Enables more accurate comparisons of expression levels between spots.

    ---

    ## *Log-Transformation*
    - **Purpose:** Stabilizes variance and prevents dominance by highly expressed genes.
    - **Methodology:**
      - A log-transformation is applied after adding a pseudocount of 1 to avoid undefined values for `log2(0)`.
      - Helps render expression changes symmetrical for easier comparisons.
    - **Benefit:** Facilitates downstream analyses by reducing biases caused by extreme values.

    ---

    ## *STx vs. scRNA-seq Normalization*
    Normalization in STx differs from scRNA-seq due to the following challenges:
    - **Spots Containing Multiple Cell Types:**
      - Unlike single-cell data, STx spots often represent a mixture of cell types, making cell-level normalization less effective.
    - **Clustering Variability:**
      - Datasets often span multiple tissue samples, introducing variability in clustering and expression patterns.
    - **Importance in STx:**
      - Proper normalization is foundational for reliable downstream analyses, as it ensures systematic biases are minimized.

    ---
    """
    )
    return


@app.cell
def _(adata_3, pd, sc):
    adata_3.raw = adata_3.raw if adata_3.raw else adata_3.copy()
    normalization_methods = {
        "Raw Counts": adata_3.raw.to_adata(),
        "Normalize Total (default)": adata_3.raw.to_adata(),
        "Exclude Highly Expressed": adata_3.raw.to_adata(),
        "SCTransform": adata_3.raw.to_adata(),
        "Log-Normalization": adata_3.raw.to_adata(),
        "Quantile Normalization": adata_3.raw.to_adata(),
    }
    sc.pp.normalize_total(
        normalization_methods["Normalize Total (default)"], target_sum=1000000.0
    )
    sc.pp.normalize_total(
        normalization_methods["Exclude Highly Expressed"],
        target_sum=1000000.0,
        exclude_highly_expressed=True,
    )
    try:
        sc.experimental.pp.normalize_pearson_residuals(
            normalization_methods["SCTransform"]
        )
    except RuntimeWarning:
        sc.pp.filter_genes(normalization_methods["SCTransform"], min_counts=1)
        sc.experimental.pp.normalize_pearson_residuals(
            normalization_methods["SCTransform"]
        )
    sc.pp.log1p(normalization_methods["Log-Normalization"])
    adata_qnorm = normalization_methods["Quantile Normalization"]
    adata_qnorm.X = (
        pd.DataFrame(adata_qnorm.X.toarray())
        .rank(method="average", axis=0)
        .to_numpy()
    )
    return (normalization_methods,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    <div style="border: 1px solid #ffa6a6; padding: 10px; border-radius: 5px;">
    <span style="color: #ff6666; font-size: 20px;"><b>Reflection Point:</b></span> <span style="font-size: 20px;">Library size threshold plot</span>  
    <ul>
        <li>Can you check in Enrichr how relevant these genes are https://maayanlab.cloud/Enrichr/?</li>

    </div>
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# 2.3 Selecting genes""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### 2.3.1 Background""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Gene selection - or alternatively “feature selection” - is applied to identify genes that are likely to be informative for downstream analyses. The most common feature selection method is the definition of highly variable genes (HVGs). The assumption is that since we quality-controlled and normalised our dataset, the genes with high variability are the ones that contain high levels of biological variability too. Since here we have a spatial dataset we can also try to identify spatially variable genes too (SVGs).

    It is important to note that HVGs are identified solely from the gene expression data. Spatial information does not play a role in finding HVGs. STx data pose a dilemma; does the meaningful spatial information reflect only spatial distribution of major cell types or does it reflect additional important spatial features? If we believe the former, relying on HVGs can be enough. If the second also holds true though, it is important to identify SVGs as well.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### 2.3.2 Highly Variable Genes (HVGs)""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Right now, let's check how HVG are affected by the previous normalisation strategies:""")
    return


@app.cell
def _(issparse, normalization_methods, plt, sc):
    for method, norm_adata in normalization_methods.items():
        if method in ["Raw Counts", "Quantile Normalization"]:
            print(
                f"Skipping {method} due to data issues or unsuitability for HVG detection."
            )
            continue
        print(f"Processing {method}:")
        if issparse(norm_adata.X):
            gene_filter = norm_adata.X.mean(axis=0).A1 > 1e-12
        else:
            gene_filter = norm_adata.X.mean(axis=0) > 1e-12
        norm_adata = norm_adata[:, gene_filter]
        sc.pp.highly_variable_genes(
            norm_adata,
            flavor="cell_ranger",
            n_top_genes=int(0.1 * norm_adata.shape[1]),
            subset=False,
        )
        hvg_count = norm_adata.var["highly_variable"].sum()
        print(f"Number of Highly Variable Genes for {method}: {hvg_count}")
        top_20_hvgs = (
            norm_adata.var[norm_adata.var["highly_variable"]]
            .sort_values(by="dispersions_norm", ascending=False)
            .head(20)
        )
        print(f"Top 20 Highly Variable Genes for {method}:")
        print(top_20_hvgs.index.tolist())
        plt.figure(figsize=(10, 6))
        plt.scatter(
            norm_adata.var["means"],
            norm_adata.var["dispersions_norm"],
            alpha=0.6,
            s=10,
            label="All Genes",
        )
        plt.scatter(
            top_20_hvgs["means"],
            top_20_hvgs["dispersions_norm"],
            color="red",
            label="Top 20 HVGs",
        )
        for gene in top_20_hvgs.index:
            mean = top_20_hvgs.loc[gene, "means"]
            disp = top_20_hvgs.loc[gene, "dispersions_norm"]
            plt.text(mean, disp, gene, fontsize=8, alpha=0.75)
        plt.xlabel("Mean Expression")
        plt.ylabel("Normalized Dispersion")
        plt.title(f"HVGs - {method}")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.show()
        print(f"Top 20 HVGs for {method} (printed as a list):")
        for i, gene in enumerate(top_20_hvgs.index.tolist(), start=1):
            print(f"{i}. {gene}")
        print("-" * 50)
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    <div style="border: 1px solid #ffa6a6; padding: 10px; border-radius: 5px;">
    <span style="color: #ff6666; font-size: 20px;"><b>Reflection Point:</b></span> <span style="font-size: 20px;">HVG and biological on-the-fly evaluation</span>  
    <ul>
        <li>Which of these gene-sets appear more biological plausible based on prior knowledge (e.g., https://maayanlab.cloud/Enrichr/)?</li>

    </div>
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# 2.4 Dimensionality reduction and clustering""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Now we can check how the different normalisation strategies are affecting dimensionality reduction and clustering:""")
    return


@app.cell
def _(normalization_methods, plt, sc):
    umap_results = {}
    for method_1, norm_adata_1 in normalization_methods.items():
        print(f"Processing normalization method: {method_1}")
        sc.tl.pca(norm_adata_1)
        norm_adata_1.obsm["X_pca"] = norm_adata_1.obsm["X_pca"][:, :2]
        sc.pp.neighbors(norm_adata_1, n_neighbors=15)
        sc.tl.leiden(norm_adata_1, resolution=0.2)
        sc.tl.umap(norm_adata_1)
        umap_results[method_1] = {
            "pca": norm_adata_1.obsm["X_pca"],
            "umap": norm_adata_1.obsm["X_umap"],
            "clusters": norm_adata_1.obs["leiden"],
        }
        plt.figure(figsize=(6, 6))
        plt.scatter(
            norm_adata_1.obsm["X_pca"][:, 0],
            norm_adata_1.obsm["X_pca"][:, 1],
            c=norm_adata_1.obs["leiden"].astype(int),
            cmap="tab20",
            s=10,
            alpha=0.8,
        )
        plt.title(f"PCA Clustering ({method_1})")
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.colorbar(label="Cluster")
        plt.show()
        sc.pl.umap(
            norm_adata_1,
            color="leiden",
            title=f"UMAP Clustering ({method_1})",
            show=False,
            legend_loc="on data",
        )
        plt.show()
        if "spatial" in norm_adata_1.uns:
            sc.pl.spatial(
                norm_adata_1,
                img_key="hires",
                color="leiden",
                size=1.5,
                title=f"Spatial Clustering ({method_1})",
            )
        norm_adata_1.var_names_make_unique()
        sc.tl.rank_genes_groups(norm_adata_1, groupby="leiden", method="t-test")
        print(f"Plotting heatmap for cluster ({method_1})...")
        sc.pl.rank_genes_groups_heatmap(
            norm_adata_1, n_genes=10, groupby="leiden", show=False
        )
        plt.title(f"Top Marker Genes for Cluster ({method_1})")
        plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    <div style="border: 1px solid #ffa6a6; padding: 10px; border-radius: 5px;">
    <span style="color: #ff6666; font-size: 20px;"><b>Reflection Point:</b></span> <span style="font-size: 20px;">Normalisation and downstream processing</span>  
    <ul>
        <li>With which normalisation strategy do you see more plausible clusterings based on the ground-truth image?</li>

    </div>
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    <div style="border: 1px solid #ffa6a6; padding: 10px; border-radius: 5px;">
    <span style="color: #ff6666; font-size: 20px;"><b>Reflection Point:</b></span> <span style="font-size: 20px;">Normalisation and downstream processing</span>  
    <ul>
        <li>Re-run the pipeline with different QC; more conservative or more extreme, what do you see?</li>

    </div>
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# 2.5 SpatialData and Visium_HD""")
    return


@app.cell
def _(visium_hd):
    sdata_visium_hd = visium_hd(
        "day_1/practical_0/data/visium_hd_3.0.0_io_subset",
        filtered_counts_file=False,
    )
    return (sdata_visium_hd,)


@app.cell
def _():
    # sdata_visium_hd.tables["square_008um"].var_names_make_unique()  # modifies in-place
    return


@app.cell
def _():
    # for table in sdata_visium_hd.tables.values():
    #     table.var_names_make_unique()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We will use the .query.bounding_box() function to subset the data so that we can visually inspect the result easier. This function allows us to subset the data to a specific bounding box, which is defined by the top-left and bottom-right coordinates. This is useful when working with large datasets, as it allows us to only load the data that we're interested in.""")
    return


@app.cell
def _():
    # sdata_visium_hd.tables["square_008um"].to_df().apply(
    #     lambda col: col.var(), axis=0
    # ).sort_values()
    return


@app.cell
def _():
    # sdata_visium_hd_crop.pl.render_shapes(color="Igha", method="matplotlib").pl.show(
    #     "Visium_HD_Mouse_Small_Intestine_downscaled_hires", figsize=(10, 10)
    # )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We can see that these bins are continous, so when overlaying the expression over an image, we would no longer see the image. To counter this, we will modify the colormap we'll use so that true 0s in the expression matrix are shown fully transparent. The spatialdata-plot library provides a helper function for this.""")
    return


@app.cell
def _(plt):
    from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent

    cmap = set_zero_in_cmap_to_transparent(plt.cm.viridis)
    cmap
    return cmap, set_zero_in_cmap_to_transparent


@app.cell
def _(sdata_visium_hd):
    sdata_visium_hd_crop = sdata_visium_hd.query.bounding_box(
        axes=["x", "y"],
        min_coordinate=[4500, 0],
        max_coordinate=[6000, 1500],
        target_coordinate_system="Visium_HD_Mouse_Small_Intestine_downscaled_hires",
    )
    return (sdata_visium_hd_crop,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We can now use this cmap to visualise the spatial expression of the genes. We see that the gene is only expressed in certain areas of the slide.""")
    return


@app.cell
def _(cmap, sdata_visium_hd_crop):
    sdata_visium_hd_crop.pl.render_images().pl.render_shapes(
        color="Igha", method="matplotlib", cmap=cmap
    ).pl.show("Visium_HD_Mouse_Small_Intestine_downscaled_hires", figsize=(8, 6))
    return


@app.function
# Define a function to crop a specific region
def crop_region(sdata, min_coord, max_coord):
    return sdata.query.bounding_box(
        axes=["x", "y"],
        min_coordinate=min_coord,
        max_coordinate=max_coord,
        target_coordinate_system="Visium_HD_Mouse_Small_Intestine_downscaled_hires",
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Let’s now pre-compute a lazy-rasterized version of the data using `rasterize_bins()`. This operation takes just a few seconds and unlocks very fast on-demand channel-wise rasterization.""")
    return


@app.cell
def _(np, sdata_visium_hd):
    # magic command not supported in marimo; please file an issue to add support
    # %%time
    from spatialdata import rasterize_bins

    for bin_size in ["016", "008", "002"]:
        table_key = f"square_{bin_size}um"

        # 1. Get the AnnData's X matrix (it should be sparse)
        X_matrix = sdata_visium_hd.tables[table_key].X

        # 2. Convert the matrix to CSC format
        #    (This should work as it's a basic SciPy function)
        X_matrix_csc = X_matrix.tocsc()

        # 3. Force the internal index arrays to 64-bit integers manually
        #    This prevents the original integer overflow error without using asformat().
        X_matrix_csc.indptr = X_matrix_csc.indptr.astype(np.int64)
        X_matrix_csc.indices = X_matrix_csc.indices.astype(np.int64)

        # 4. Assign the modified CSC matrix back to the AnnData object
        sdata_visium_hd.tables[table_key].X = X_matrix_csc

        # 5. Proceed with rasterization
        rasterized = rasterize_bins(
            sdata_visium_hd,
            f"Visium_HD_Mouse_Small_Intestine_square_{bin_size}um",
            table_key,
            "array_col",
            "array_row",
        )
        sdata_visium_hd[f"rasterized_{bin_size}um"] = rasterized
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    This produces lazy image objects that can be accessed gene-wise exceptionally efficiently.

    Importantly, these objects should not be computed as a whole, because this would lead to the unnecessary computation of hundreds of GB of memory.
    """
    )
    return


@app.cell
def _(sdata_visium_hd):
    # magic command not supported in marimo; please file an issue to add support
    # %%time
    # this is very fast
    gene_name = "AA986860"
    sdata_visium_hd["rasterized_002um"].sel(c=gene_name).compute()
    # this must not be called
    # sdata["rasterized_002um"].compute()
    return (gene_name,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Here is an example of how to plot the rasterized data. Please pay attention to setting scale="full"; this is essential. If it is not set, spatialdata-plot will try to re-rasterize the data to fit the canvas size and try to compute the whole object.""")
    return


@app.cell
def _(gene_name, plt, sdata_visium_hd):
    plt.figure(figsize=(10, 10))
    ax = plt.gca()
    sdata_visium_hd.pl.render_images(
        "rasterized_016um", channel=gene_name, scale="full"
    ).pl.show(coordinate_systems="Visium_HD_Mouse_Small_Intestine", ax=ax)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    As you can see below, we can plot the 2µm bins very performantly, which would not be possible with the matplotlib based approach shown before.

    Note: when we plot the full data below the plot appears as uniformly violet. This happens because the data is very sparse and the bins are too small for the target figure size; since most of the non-zero bins are 1s the final interpolated colors are very close to the zero value and difficult to see.

    Let’s therefore:

    change the limits of the plot to show a portion of the data to avoid interpolation artifacts

    plot a binary mask of the full data to avoid the interpolated data to be to close to 0.
    """
    )
    return


@app.cell
def _(gene_name, plt, sdata_visium_hd):
    from matplotlib.colors import ListedColormap

    colors = ["#000000", "#ffffff"]
    cmap_1 = ListedColormap(colors)
    plt.figure(figsize=(10, 10))
    plt.imshow(
        sdata_visium_hd["rasterized_002um"].sel(c=gene_name).data.compute(),
        cmap=cmap_1,
        vmin=0,
        vmax=0.01,
    )
    plt.colorbar()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""As a final not on the on-the-fly rasterization approach. Please, by looking at the corners of the last plot, notice how the data is on a grid that is actually sligthly rotated. The advantage of using rasterize_bins() is that the produced object contains the coordinate transformations necessary to align (rotation and scale) the rasterized data together with the high-resolution images.""")
    return


@app.cell
def _(sdata_visium_hd):
    sdata_small = sdata_visium_hd.query.bounding_box(
        min_coordinate=[7000, 11000],
        max_coordinate=[10000, 14000],
        axes=("x", "y"),
        target_coordinate_system="Visium_HD_Mouse_Small_Intestine",
    )
    return (sdata_small,)


@app.cell
def _(gene_name, sdata_small):
    sdata_small.pl.render_shapes(
        "Visium_HD_Mouse_Small_Intestine_square_016um", color=gene_name
    ).pl.show(coordinate_systems="Visium_HD_Mouse_Small_Intestine")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Notice how a Moiré pattern is visible. This is due to the fact that the grid is not axis-aligned but presents a small rotation. A solution is to switch to datashader as a backend (which is enabled by default when the number of geometries is large). This will create some artifacts (bins of different sizes, some “holes” in the plot), but generally is expected to mitigate the effects as opposed the plot above.""")
    return


@app.cell
def _(sdata_small):
    gene_name_1 = "AA986860"
    sdata_small.pl.render_shapes(
        "Visium_HD_Mouse_Small_Intestine_square_016um",
        color=gene_name_1,
        method="datashader",
    ).pl.show(
        coordinate_systems="Visium_HD_Mouse_Small_Intestine_downscaled_hires"
    )
    return


@app.cell
def _(sdata_small):
    gene_name_2 = "AA986860"
    for bin_size_1 in [16, 8, 2]:
        sdata_small.pl.render_shapes(
            f"Visium_HD_Mouse_Small_Intestine_square_{bin_size_1:03}um",
            color=gene_name_2,
            method="datashader",
        ).pl.show(
            coordinate_systems="Visium_HD_Mouse_Small_Intestine_downscaled_hires",
            title=f"bin_size={bin_size_1}µm",
            figsize=(10, 10),
        )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    The data present a lot of sparsity. Let’s remake the plots above by visualizing only the non-zero entries and using the full-resolution image as a background.

    We will do this by modifying the viridis colormap so that 0 is plotted as transparent. Let’s also truncate the viridis colormap so that the highest value is colored green and not yellow since green has a better contrast against the pink of the H&E microscopy image.
    """
    )
    return


@app.cell
def _(set_zero_in_cmap_to_transparent):
    # let's display the areas where no expression is detected as transparent
    new_cmap = set_zero_in_cmap_to_transparent(cmap="viridis")
    new_cmap
    return (new_cmap,)


@app.cell
def _(sdata_small):
    sdata_small
    return


@app.cell
def _(new_cmap, sdata_small):
    gene_name_3 = "AA986860"
    for bin_size_2 in [16, 8]:
        sdata_small.pl.render_images(
            "Visium_HD_Mouse_Small_Intestine_hires_image"
        ).pl.render_shapes(
            f"Visium_HD_Mouse_Small_Intestine_square_{bin_size_2:03}um",
            color=gene_name_3,
            cmap=new_cmap,
        ).pl.show(
            coordinate_systems="Visium_HD_Mouse_Small_Intestine_downscaled_hires",
            title=f"bin_size={bin_size_2}µm",
            figsize=(10, 10),
        )
    return (gene_name_3,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Let’s make a zoomed version of the plot for the 2µm bins to better visualize them. Please notice that we pass method='matplotlib' as currently this is required in order to visualize the background as transparent by means of the modified colormap.""")
    return


@app.cell
def _():
    from spatialdata import bounding_box_query
    return (bounding_box_query,)


@app.cell
def _(bounding_box_query, gene_name_3, new_cmap, sdata_small):
    crop1 = lambda x: bounding_box_query(
        x,
        min_coordinate=[9000, 11000],
        max_coordinate=[10000, 12000],
        axes=("x", "y"),
        target_coordinate_system="Visium_HD_Mouse_Small_Intestine",
    )
    crop1(sdata_small).pl.render_images(
        "Visium_HD_Mouse_Small_Intestine_hires_image"
    ).pl.render_shapes(
        "Visium_HD_Mouse_Small_Intestine_square_002um",
        color=gene_name_3,
        cmap=new_cmap,
        method="matplotlib",
    ).pl.show(
        coordinate_systems="Visium_HD_Mouse_Small_Intestine_downscaled_hires",
        title=f"bin_size=2µm",
        figsize=(10, 10),
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""As you can see the 8µm bins are convenient for looking at gene expression distribution from a broad perspective (same for the 16µm bins, where some resolution can be sacrificed in exchange for a faster visualization). On the other hand, the 2µm bins allow to precisely locate the expressed genes in the tissue.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Let’s now color the 16µm bins by cluster identity. Let’s reuse the clusters gene_expression_graphclust computed from 10x Genomics and available with the raw data from the 10x Genomics website.""")
    return


@app.cell
def _(pd):
    import os
    from tempfile import TemporaryDirectory

    import requests

    # For convenience we rehost the single file containing the clusters we are interested in.
    # Let's download it in a temporary directory and read it in a pandas DataFrame. The file is 2 MB.
    clusters_file_url = "https://s3.embl.de/spatialdata/misc/visium_hd_mouse_intestine_16um_graphclust.csv"

    with TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "data.csv")
        response = requests.get(clusters_file_url)
        with open(path, "wb") as f:
            f.write(response.content)
        df = pd.read_csv(path)
    return (df,)


@app.cell
def _(df):
    df.head(3)
    print(df.columns.tolist())
    return


@app.cell
def _(df):
    # let's convert the Cluster dtype from int64 to categorical since later we want the plots to use a categorical colormap
    df["Cluster"] = df["Cluster"].astype("category")
    # df.set_index("Barcode", inplace=True)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Let’s merge the data.""")
    return


@app.cell
def _(df, sdata_visium_hd):
    sdata_visium_hd["square_016um"].obs["Cluster"] = df["Cluster"]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Let’s plot the clusters on one of the data crops we used before.""")
    return


@app.cell
def _(sdata_visium_hd):
    sdata_visium_hd
    return


@app.cell
def _(bounding_box_query, plt, sdata_visium_hd):
    fig_1, ax_1 = plt.subplots(1, 1, figsize=(10, 5))
    crop0 = lambda x: bounding_box_query(
        x,
        min_coordinate=[5000, 8000],
        max_coordinate=[10000, 13000],
        axes=("x", "y"),
        target_coordinate_system="Visium_HD_Mouse_Small_Intestine",
    )
    crop0(sdata_visium_hd).pl.render_images(
        "Visium_HD_Mouse_Small_Intestine_hires_image"
    ).pl.show(
        ax=ax_1,
        title="Full image",
        coordinate_systems="Visium_HD_Mouse_Small_Intestine_downscaled_hires",
    )
    return (crop0,)


@app.cell
def _(crop0, sdata_visium_hd):
    # magic command not supported in marimo; please file an issue to add support
    # %%time
    crop0(sdata_visium_hd).pl.render_images(
        "Visium_HD_Mouse_Small_Intestine_downscaled_hires_image"
    ).pl.render_shapes(
        "Visium_HD_Mouse_Small_Intestine_square_016um", color="Cluster"
    ).pl.show(
        coordinate_systems="Visium_HD_Mouse_Small_Intestine_downscaled_hires",
        title=f"bin_size=016µm",
        figsize=(10, 10),
    )
    return


@app.cell
def _(sdata_visium_hd):
    # Access the specific AnnData table
    adata_table3 = sdata_visium_hd.tables[
        "square_016um"
    ]  # Replace "table_name" with the actual table name

    # Make variable names unique
    adata_table3.var_names_make_unique()

    # Reassign the updated table back to the SpatialData object (optional if you need to reuse it)
    sdata_visium_hd.tables["square_016um"] = adata_table3
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## QC experimental""")
    return


@app.cell
def _(sdata_visium_hd):
    sdata_small_1 = sdata_visium_hd.query.bounding_box(
        min_coordinate=[7000, 11000],
        max_coordinate=[10000, 14000],
        axes=("x", "y"),
        target_coordinate_system="Visium_HD_Mouse_Small_Intestine",
    )
    return


@app.cell
def _(sdata_visium_hd_crop):
    adata_table3_1 = sdata_visium_hd_crop.tables["square_016um"]
    adata_table3_1.var_names_make_unique()
    sdata_visium_hd_crop.tables["square_016um"] = adata_table3_1
    return


@app.cell
def _(sdata_visium_hd):
    # Check initial dimensions
    print(f"Initial dataset dimensions: {sdata_visium_hd['square_016um'].shape}")

    # Filter on-tissue spots
    on_tissue_mask = sdata_visium_hd["square_016um"].obs["in_tissue"] == 1
    sdata_visium_hd["square_016um"] = sdata_visium_hd["square_016um"][
        on_tissue_mask, :
    ]

    # Verify dimensions after filtering
    print(f"Filtered dataset dimensions: {sdata_visium_hd['square_016um'].shape}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### 1. Mito-tracking""")
    return


@app.cell
def _(sdata_visium_hd):
    print("First 20 gene names in the dataset:")
    print(sdata_visium_hd["square_016um"].var.index[:20].tolist())
    mito_genes_alt_1 = sdata_visium_hd["square_016um"].var.index.str.contains(
        "(^MT-)|(^mt-)|(^mt\\\\-)", regex=True
    )
    print(f"Number of mitochondrial genes detected: {mito_genes_alt_1.sum()}")
    mito_genes_list_1 = sdata_visium_hd["square_016um"].var.index[mito_genes_alt_1]
    print("Mitochondrial genes detected:")
    print(mito_genes_list_1.tolist())
    return


@app.cell
def _(sdata_visium_hd):
    bin_sizes = ["square_002um", "square_008um", "square_016um"]
    for bin_size_3 in bin_sizes:
        print(f"Processing bin: {bin_size_3}")
        adata_4 = sdata_visium_hd[bin_size_3]
        print(f"First 20 gene names in {bin_size_3}:")
        print(adata_4.var.index[:20].tolist())
        mito_genes_alt_2 = adata_4.var.index.str.contains(
            "(^MT-)|(^mt-)|(^mt\\\\-)", regex=True
        )
        print(
            f"Number of mitochondrial genes detected in {bin_size_3}: {mito_genes_alt_2.sum()}"
        )
        mito_genes_list_2 = adata_4.var.index[mito_genes_alt_2]
        print(f"Mitochondrial genes detected in {bin_size_3}:")
        print(mito_genes_list_2.tolist())
        print("-" * 40)
    return (mito_genes_alt_2,)


@app.cell
def _(mito_genes_alt_2, sdata_visium_hd):
    sdata_visium_hd["square_016um"].var["mt"] = mito_genes_alt_2
    from scanpy.preprocessing import calculate_qc_metrics

    calculate_qc_metrics(
        sdata_visium_hd["square_016um"],
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    print("Mitochondrial QC metrics:")
    print(
        sdata_visium_hd["square_016um"]
        .obs[["pct_counts_mt", "total_counts", "n_genes_by_counts"]]
        .head()
    )
    return (calculate_qc_metrics,)


@app.cell
def _(calculate_qc_metrics, sdata_visium_hd):
    bin_sizes_1 = ["square_002um", "square_008um", "square_016um"]
    for bin_size_4 in bin_sizes_1:
        print(f"Processing bin: {bin_size_4}")
        adata_5 = sdata_visium_hd[bin_size_4]
        mito_genes_alt_3 = adata_5.var.index.str.contains(
            "(^MT-)|(^mt-)|(^mt\\\\-)", regex=True
        )
        adata_5.var["mt"] = mito_genes_alt_3
        calculate_qc_metrics(
            adata_5, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
        print(f"Mitochondrial QC metrics for {bin_size_4}:")
        print(
            adata_5.obs[
                ["pct_counts_mt", "total_counts", "n_genes_by_counts"]
            ].head()
        )
        print("-" * 40)
    return


@app.cell
def _(plt, sdata_visium_hd, sns):
    plt.figure(figsize=(8, 6))
    sns.histplot(
        sdata_visium_hd["square_016um"].obs["pct_counts_mt"],
        kde=True,
        bins=50,
        color="gray",
    )
    plt.title("Percentage of Mitochondrial Counts per Spot")
    plt.xlabel("Mitochondrial Counts (%)")
    plt.ylabel("Frequency")
    plt.legend()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### 2. Total counts""")
    return


@app.cell
def _(sdata_visium_hd):
    # List all available elements in the SpatialData object
    print(sdata_visium_hd)

    # Check if 'total_counts' exists in any table
    if "square_016um" in sdata_visium_hd.tables:
        print("Columns in 'square_016um':")
        print(sdata_visium_hd["square_016um"].obs.columns)
    return


@app.cell
def _(plt, sdata_visium_hd, sns):
    total_counts = sdata_visium_hd["square_016um"].obs["total_counts"]
    sns.histplot(total_counts, kde=True, bins=50, color="blue")
    plt.title("Total Counts per Spot")
    plt.xlabel("Total Counts")
    plt.ylabel("Density")
    plt.legend()
    plt.show()
    return


@app.cell
def _(sdata_visium_hd):
    for bin_size_5 in ["square_002um", "square_008um", "square_016um"]:
        print(f"\nProcessing bin size: {bin_size_5}")
        adata_6 = sdata_visium_hd[bin_size_5]
        if "total_counts" in adata_6.obs.columns:
            total_counts_1 = adata_6.obs["total_counts"]
            print("Basic statistics for total counts per spot:")
            print(total_counts_1.describe())
            low_library_size_threshold = 1000
            low_library_spots = total_counts_1[
                total_counts_1 < low_library_size_threshold
            ]
            print(
                f"Number of spots below the threshold of {low_library_size_threshold}: {len(low_library_spots)}"
            )
            print("Example spots with low total counts:")
            print(low_library_spots.head())
        else:
            print(
                f"'total_counts' column not found for bin size {bin_size_5}. Skipping."
            )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### 3. Genes""")
    return


@app.cell
def _(plt, sdata_visium_hd, sns):
    plt.figure(figsize=(8, 6))
    sns.histplot(
        sdata_visium_hd["square_016um"].obs["n_genes_by_counts"],
        kde=True,
        bins=50,
        color="gray",
    )
    plt.title("Number of Detected Genes per Spot")
    plt.xlabel("Number of Genes")
    plt.ylabel("Frequency")
    plt.legend()
    plt.show()
    return


@app.cell
def _(plt, sdata_visium_hd):
    plt.figure(figsize=(10, 10))
    ax_2 = plt.gca()
    sdata_visium_hd.pl.render_shapes(
        "Visium_HD_Mouse_Small_Intestine_square_016um",
        color="total_counts",
        method="datashader",
    ).pl.show(
        coordinate_systems="Visium_HD_Mouse_Small_Intestine_downscaled_hires",
        ax=ax_2,
    )
    return


@app.cell
def _(plt, sdata_visium_hd):
    plt.figure(figsize=(10, 10))
    ax_3 = plt.gca()
    sdata_visium_hd.pl.render_shapes(
        "Visium_HD_Mouse_Small_Intestine_square_016um",
        color="n_genes_by_counts",
        method="datashader",
    ).pl.show(
        coordinate_systems="Visium_HD_Mouse_Small_Intestine_downscaled_hires",
        ax=ax_3,
    )
    plt.title("Number of Genes Detected")
    plt.figure(figsize=(10, 10))
    ax_3 = plt.gca()
    sdata_visium_hd.pl.render_shapes(
        "Visium_HD_Mouse_Small_Intestine_square_016um",
        color="pct_counts_mt",
        method="datashader",
    ).pl.show(
        coordinate_systems="Visium_HD_Mouse_Small_Intestine_downscaled_hires",
        ax=ax_3,
    )
    plt.title("Percentage of Mitochondrial Reads")
    return


@app.cell
def _(pd, sdata_visium_hd):
    bin_sizes_2 = ["square_002um", "square_008um", "square_016um"]
    metrics = {}
    for bin_size_6 in bin_sizes_2:
        adata_7 = sdata_visium_hd[bin_size_6]
        metrics[bin_size_6] = pd.DataFrame(
            {
                "total_counts": adata_7.obs["total_counts"],
                "n_genes_by_counts": adata_7.obs["n_genes_by_counts"],
                "pct_counts_mt": adata_7.obs["pct_counts_mt"],
            }
        )
    return bin_sizes_2, metrics


@app.cell
def _(bin_sizes_2, metrics, plt, sns):
    for metric_1 in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
        plt.figure(figsize=(10, 6))
        for bin_size_7 in bin_sizes_2:
            sns.kdeplot(
                metrics[bin_size_7][metric_1],
                label=bin_size_7,
                fill=True,
                alpha=0.4,
            )
        plt.title(f"Distribution of {metric_1} Across Bin Sizes")
        plt.xlabel(metric_1)
        plt.ylabel("Density")
        plt.legend(title="Bin Size")
        plt.show()
    return


@app.cell
def _(bin_sizes_2, metrics, pd, plt, sns):
    combined_data = pd.concat(
        [metrics[bin_size].assign(bin_size=bin_size) for bin_size in bin_sizes_2]
    )
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=combined_data, x="bin_size", y="total_counts", palette="Set2")
    plt.title("Total Counts Across Bin Sizes")
    plt.xlabel("Bin Size")
    plt.ylabel("Total Counts")
    plt.show()
    return


@app.cell
def _(bin_sizes_2, metrics):
    for metric_2 in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
        print(f"Variance of {metric_2} across bin sizes:")
        for bin_size_8 in bin_sizes_2:
            variance = metrics[bin_size_8][metric_2].var()
            print(f"  {bin_size_8}: {variance:.2f}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    <div style="border: 1px solid #ffa6a6; padding: 10px; border-radius: 5px;">
    <span style="color: #ff6666; font-size: 20px;"><b>Reflection Point:</b></span> <span style="font-size: 20px;">QC across bins</span>  
    <ul>
        <li>Are you seeing differences across bins?</li>

    </div>
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Let's remove:""")
    return


@app.cell
def _(sdata_visium_hd):
    sdata_visium_hd
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Remove cells""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Below, some indicative filtering (please change based on more accurate observation of QCs)""")
    return


@app.cell
def _(sdata_visium_hd):
    # Define thresholds for filtering
    mito_threshold = 20  # Maximum percentage of mitochondrial genes allowed
    min_counts = 500  # Minimum total counts per spot
    max_counts = 50000  # Maximum total counts per spot
    min_genes = 200  # Minimum number of genes detected per spot

    # Apply filters to the data
    filtered_spots = (
        (sdata_visium_hd["square_016um"].obs["pct_counts_mt"] < mito_threshold)
        & (sdata_visium_hd["square_016um"].obs["total_counts"] > min_counts)
        & (sdata_visium_hd["square_016um"].obs["total_counts"] < max_counts)
        & (sdata_visium_hd["square_016um"].obs["n_genes_by_counts"] > min_genes)
    )

    # Subset the dataset to keep only high-quality spots
    sdata_visium_hd["square_016um"] = sdata_visium_hd["square_016um"][
        filtered_spots
    ]

    # Confirm filtering results
    print("After filtering:")
    print(f"Number of spots retained: {sdata_visium_hd['square_016um'].n_obs}")
    print(f"Number of genes retained: {sdata_visium_hd['square_016um'].n_vars}")
    return


@app.cell
def _(sdata_visium_hd):
    sdata_visium_hd
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    <div style="border: 1px solid #ffa6a6; padding: 10px; border-radius: 5px;">
    <span style="color: #ff6666; font-size: 20px;"><b>Reflection Point:</b></span> <span style="font-size: 20px;">Last challenge</span>  
    <ul>
        <li>Can you make the following code work for Visium HD?</li>

    </div>
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Normalisation""")
    return


@app.cell
def _(sdata_visium_hd):
    sadata = sdata_visium_hd["square_016um"]
    return (sadata,)


@app.cell
def _(sadata):
    sadata
    return


@app.cell
def _(pd, sadata, sc):
    normalization_methods_1 = {
        "Raw Counts": sadata.copy(),
        "Normalize Total (default)": sadata.copy(),
        "Exclude Highly Expressed": sadata.copy(),
        "SCTransform": sadata.copy(),
        "Log-Normalization": sadata.copy(),
        "Quantile Normalization": sadata.copy(),
    }
    sc.pp.normalize_total(
        normalization_methods_1["Normalize Total (default)"], target_sum=1000000.0
    )
    sc.pp.normalize_total(
        normalization_methods_1["Exclude Highly Expressed"],
        target_sum=1000000.0,
        exclude_highly_expressed=True,
    )
    try:
        sc.experimental.pp.normalize_pearson_residuals(
            normalization_methods_1["SCTransform"]
        )
    except RuntimeWarning:
        sc.pp.filter_genes(normalization_methods_1["SCTransform"], min_counts=1)
        sc.experimental.pp.normalize_pearson_residuals(
            normalization_methods_1["SCTransform"]
        )
    sc.pp.log1p(normalization_methods_1["Log-Normalization"])
    quantile_sadata = normalization_methods_1["Quantile Normalization"]
    quantile_sadata.X = (
        pd.DataFrame(quantile_sadata.X.toarray())
        .rank(method="average", axis=0)
        .to_numpy()
    )
    return (normalization_methods_1,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## HVG""")
    return


@app.cell
def _(issparse, normalization_methods_1, plt, sc):
    for method_2, sadata_1 in normalization_methods_1.items():
        if method_2 in ["Raw Counts", "Quantile Normalization"]:
            print(
                f"Skipping {method_2} due to data issues or unsuitability for HVG detection."
            )
            continue
        print(f"Processing {method_2}:")
        if issparse(sadata_1.X):
            gene_filter_1 = sadata_1.X.mean(axis=0).A1 > 1e-12
        else:
            gene_filter_1 = sadata_1.X.mean(axis=0) > 1e-12
        sadata_1 = sadata_1[:, gene_filter_1]
        sc.pp.highly_variable_genes(
            sadata_1,
            flavor="cell_ranger",
            n_top_genes=int(0.1 * sadata_1.shape[1]),
            subset=False,
        )
        hvg_count_1 = sadata_1.var["highly_variable"].sum()
        print(f"Number of Highly Variable Genes for {method_2}: {hvg_count_1}")
        top_20_hvgs_1 = (
            sadata_1.var[sadata_1.var["highly_variable"]]
            .sort_values(by="dispersions_norm", ascending=False)
            .head(20)
        )
        print(f"Top 20 Highly Variable Genes for {method_2}:")
        print(top_20_hvgs_1.index.tolist())
        plt.figure(figsize=(10, 6))
        plt.scatter(
            sadata_1.var["means"],
            sadata_1.var["dispersions_norm"],
            alpha=0.6,
            s=10,
            label="All Genes",
        )
        plt.scatter(
            top_20_hvgs_1["means"],
            top_20_hvgs_1["dispersions_norm"],
            color="red",
            label="Top 20 HVGs",
        )
        for gene_1 in top_20_hvgs_1.index:
            mean_1 = top_20_hvgs_1.loc[gene_1, "means"]
            disp_1 = top_20_hvgs_1.loc[gene_1, "dispersions_norm"]
            plt.text(mean_1, disp_1, gene_1, fontsize=8, alpha=0.75)
        plt.xlabel("Mean Expression")
        plt.ylabel("Normalized Dispersion")
        plt.title(f"HVGs - {method_2}")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.show()
        print(f"Top 20 HVGs for {method_2} (printed as a list):")
        for i_1, gene_1 in enumerate(top_20_hvgs_1.index.tolist(), start=1):
            print(f"{i_1}. {gene_1}")
        print("-" * 50)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Dimensionality Reduction""")
    return


@app.cell
def _(normalization_methods_1, plt, sc):
    umap_results_1 = {}
    for method_3, sadata_2 in normalization_methods_1.items():
        print(f"Processing normalization method: {method_3}")
        sc.tl.pca(sadata_2)
        sadata_2.obsm["X_pca"] = sadata_2.obsm["X_pca"][:, :2]
        sc.pp.neighbors(sadata_2, n_neighbors=15)
        sc.tl.leiden(sadata_2, resolution=0.2)
        sc.tl.umap(sadata_2)
        umap_results_1[method_3] = {
            "pca": sadata_2.obsm["X_pca"],
            "umap": sadata_2.obsm["X_umap"],
            "clusters": sadata_2.obs["leiden"],
        }
        plt.figure(figsize=(6, 6))
        plt.scatter(
            sadata_2.obsm["X_pca"][:, 0],
            sadata_2.obsm["X_pca"][:, 1],
            c=sadata_2.obs["leiden"].astype(int),
            cmap="tab20",
            s=10,
            alpha=0.8,
        )
        plt.title(f"PCA Clustering ({method_3})")
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.colorbar(label="Cluster")
        plt.show()
        sc.pl.umap(
            sadata_2,
            color="leiden",
            title=f"UMAP Clustering ({method_3})",
            show=False,
            legend_loc="on data",
        )
        plt.show()
        if "spatial" in sadata_2.uns:
            sc.pl.spatial(
                sadata_2,
                img_key="hires",
                color="leiden",
                size=1.5,
                title=f"Spatial Clustering ({method_3})",
            )
        sadata_2.var_names_make_unique()
        sc.tl.rank_genes_groups(sadata_2, groupby="leiden", method="t-test")
        cluster_to_plot = "9"
        print(f"Plotting heatmap for cluster {cluster_to_plot} ({method_3})...")
        sc.pl.rank_genes_groups_heatmap(
            sadata_2,
            groups=cluster_to_plot,
            n_genes=10,
            groupby="leiden",
            show=False,
        )
        plt.title(f"Top Marker Genes for Cluster {cluster_to_plot} ({method_3})")
        plt.show()
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
