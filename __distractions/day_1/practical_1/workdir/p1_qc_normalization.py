import marimo

__generated_with = "0.15.2"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## ELIXIR Spatial Transcriptomics Course
    ### Practical 1a: Imaging-based spatial transcriptomics data QC and normalization 
    Date: 2025-01-21

    Author(s): Rasool Saghaleyni, Åsa Björklund

    Author(s) email: <rasool.saghaleyni@scilifelab.se>, <asa.bjorklund@scilifelab.se>

    ⚠️ Note: The proper environment for this notebook is `p1_qc_normalization`. It can be activated by selecting the kernel in the Jupyter notebook.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Loading packages""")
    return


@app.cell
def _():
    import numpy as np
    import pandas as pd
    import os
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import seaborn as sns
    import scanpy as sc
    import squidpy as sq
    import scipy.sparse as sp
    from scipy.stats import gaussian_kde
    from scipy.signal import find_peaks, argrelextrema
    from scipy.sparse import issparse, csr_matrix
    from diptest import diptest
    import statsmodels.api as sm
    from shapely.geometry import Point, Polygon, MultiPoint
    from scipy.spatial import ConvexHull
    import sys
    sys.path.append('day_1/practical_1/workdir/custom')
    import tenx_method_nb_helper_functions as hf
    return hf, np, os, pd, plt, sc, sns, sp


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Dataset description and the design of the experiment
    Before starting the analysis, please make that you know about the biological background, experimental design and the data structure for the datset that we will do the analysis on it. Here you can find good information about this: https://pages.10xgenomics.com/rs/446-PBO-704/images/10x_LIT000210_App-Note_Xenium-In-Situ_Letter_Digital.pdf
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Loading data and primary inspections""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Making adata object considering the transcripts that should be used for the analysis.
    Running analysis on the transcripts that are only in the nucleus. 

    Data was downloaded from 10x at https://cf.10xgenomics.com/samples/xenium/1.4.0/Xenium_V1_FFPE_TgCRND8_17_9_months/Xenium_V1_FFPE_TgCRND8_17_9_months_outs.zip and unzipped in folder `data/`.
    """
    )
    return


@app.cell
def _(hf, os, pd):
    sample_path = "day_1/practical_1/data/Xenium_V1_FFPE_TgCRND8_17_9_months_outs"
    transcripts_csv_path = os.path.join(sample_path, "transcripts.csv.gz")
    transcripts_df = pd.read_csv(transcripts_csv_path, compression='gzip')
    nucleus_boundaries_gz_path = os.path.join(sample_path, "nucleus_boundaries.csv.gz")
    nucleus_df = pd.read_csv(nucleus_boundaries_gz_path, compression='gzip')
    if not os.path.exists(os.path.join(sample_path, "cells.csv")):
        hf.decompress_file(os.path.join(sample_path, "cells.csv.gz"))
    return (sample_path,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Making the adata object""")
    return


@app.cell
def _(hf, sample_path):
    adata = hf.create_adata(sample_path, nucleus_genes_only = False)
    adata
    return (adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""The `AnnData` object, `adata`, contains a structured dataset with cells as observations (`n_obs = 62268`) and genes as variables (`n_vars = 347`). Here's a breakdown of the main components within `adata`:""")
    return


@app.cell
def _(adata):
    adata.obs
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    `obs` (Observations): This table contains metadata about each cell, where each row corresponds to a cell, and each column holds information about a specific attribute:

    - `cell_id`: Unique identifier for each cell.

    - `x_centroid` and `y_centroid`: Coordinates of each cell’s center in the spatial layout, indicating where each cell is located within the tissue.

    - `transcript_counts`: Total transcript counts for each cell, showing the overall gene expression level.
    - `control_probe_counts` and `control_codeword_counts`: Counts related to control probes and codewords, which are often used for quality control in spatial transcriptomics.
    - `unassigned_codeword_counts` and `deprecated_codeword_counts`: Counts of unassigned or deprecated codewords, indicating low-confidence or outdated identifiers.
    - `total_counts`: Total counts across all measured attributes, representing the cell’s total signal.
    - `cell_area` and `nucleus_area`: Physical measurements of the cell and its nucleus area, in pixels or micrometers.
    """
    )
    return


@app.cell
def _(adata):
    adata.var
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    `var` (Variables): This table contains metadata about each gene, where each row is a gene and each column is an attribute:

    - `gene_ids`: Unique identifiers for each gene, often in Ensembl or another standardized format.

    - `feature_types`: Type of feature associated with each gene, such as "gene" or "transcript."

    - `genome`: Information about the genome source of each gene, like "human" or "mouse."

    `obsm` (Multi-dimensional Observations): This slot contains multi-dimensional data related to cells. Here, `spatial` stores spatial coordinates for each cell, allowing visualization and spatial analysis of cells in their tissue context.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Detection statistics

    We calculate the total transcript counts for each cell by summing gene expression values across all genes and the number of genes detected in each cell. The resulting distributions gives us an overview of the data quality and cellular diversity in terms of RNA content.

    These provide a measure of each cell's transcriptional activity or RNA content. High total counts typically indicate cells with higher transcriptional activity, while very low total counts may suggest low-quality cells or empty spots with minimal RNA.
    Examining this distribution helps us assess data quality and identify potential outliers:

    - Cells with Very Low Counts: These may represent low-quality cells or background noise, which could be filtered out in subsequent steps to improve analysis accuracy.
    - Cells with Very High Counts: High total counts may indicate cell types with naturally high transcriptional activity or potential doublets (two cells counted as one).

    This plot provides a quick check of the dataset’s quality and helps inform any initial filtering steps. A typical goal is to ensure the data has a reasonable distribution of RNA counts per cell, without excessive noise or artifacts that might skew downstream analysis.
    """
    )
    return


@app.cell
def _(adata, plt, sns):
    adata.obs['total_counts'] = adata.X.sum(axis=1)
    plt.figure(figsize=(8, 5))
    sns.histplot(adata.obs['total_counts'], kde='True', bins=50)
    plt.xlabel("Total Transcript Counts per Cell")
    plt.ylabel("Number of Cells")
    plt.title("Distribution of Total Transcript Counts per Cell")
    plt.show()
    return


@app.cell
def _(adata, plt, sns):
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
    plt.figure(figsize=(8, 5))
    sns.histplot(adata.obs['n_genes'], kde=True, bins=50)
    plt.xlabel("Number of Genes Detected per Cell")
    plt.ylabel("Number of Cells")
    plt.title("Distribution of Genes Detected per Cell")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Now lets find genes with the highest expression across the whole dataset. We use the scanpy function `pl.highest_expr_genes` to understand which genes dominate the transcriptional landscape. This helps us quickly identify genes with the highest abundance, which are often either essential housekeeping genes or specific markers that define particular cell types. Examining the top expressed genes serves both technical and biological purposes: it allows us to check for any potential technical artifacts (e.g., genes with unusually high background expression) and offers biological insight by highlighting key genes likely involved in core cellular functions or distinguishing cell types.""")
    return


@app.cell
def _(adata, sc):
    sc.pl.highest_expr_genes(adata)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Segmentation statistics

    We can also look into the distributions of cell and nucleus areas to assess segmentation quality and examine cell size diversity across the dataset. These distributions provide an overview of the range of cell and nucleus sizes, helping to identify segmentation artifacts or inconsistencies, such as unusually small areas (which may indicate partial cells or segmentation errors) or large areas (potentially indicating doublets or multiplets of cells).
    """
    )
    return


@app.cell
def _(adata, plt, sns):
    #distribution of cell and nucleus areas
    plt.figure(figsize=(8, 5))
    sns.histplot(adata.obs['cell_area'], kde=True, bins=50)
    plt.xlabel("Cell Area")
    plt.ylabel("Number of Cells")
    plt.title("Distribution of Cell Areas")
    plt.show()

    plt.figure(figsize=(8, 5))
    sns.histplot(adata.obs['nucleus_area'], kde=True, bins=50)
    plt.xlabel("Nucleus Area")
    plt.ylabel("Number of Cells")
    plt.title("Distribution of Nucleus Areas")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Cell Area vs Nucleus Area""")
    return


@app.cell
def _(adata, plt):
    plt.figure(figsize=(8, 6))
    plt.scatter(adata.obs['cell_area'], adata.obs['nucleus_area'], alpha=0.5)
    plt.xlabel("Cell Area")
    plt.ylabel("Nucleus Area")
    plt.title("Cell Area vs Nucleus Area")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Calculate cell-to-nucleus area ratio and plot it. This ratio provides a measure of the relative size of the nucleus compared to the whole cell, which can be informative for understanding cell morphology and the distribution of nuclear material within cells. A high ratio may indicate cells with large nuclei relative to their overall size, which could be relevant for cell type classification or biological interpretation. Examining this ratio helps identify potential outliers or unusual cell morphologies that may require further investigation or filtering.""")
    return


@app.cell
def _(adata, plt, sns):
    adata.obs['area_ratio'] = adata.obs['nucleus_area'] / adata.obs['cell_area']
    plt.figure(figsize=(8, 5))
    sns.histplot(adata.obs['area_ratio'], kde=True, bins=50)
    plt.xlabel("Nucleus-to-Cell Area Ratio")
    plt.ylabel("Number of Cells")
    plt.title("Distribution of Nucleus-to-Cell Area Ratios")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    We can also look at cell area vs total counts, we would normally expect more detected transcripts in the larger cell areas. But that of course depends on if the probeset covers all celltypes in the tissue well. 

    In this case we have quite a few large segmented areas with very low counts, these are most likely low quality regions of the tissue.
    """
    )
    return


@app.cell
def _(adata, plt):
    plt.figure(figsize=(8, 6))
    plt.scatter(adata.obs['cell_area'], adata.obs['total_counts'], alpha=0.5)
    plt.xlabel("Cell Area")
    plt.ylabel("Total Transcript Counts")
    plt.title("Cell Area vs Total Transcript Counts")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Any statistics you have may also be visualized on the tissue to understand more about the quality of different regions. You may expect smaller/larger cell sizes in different regions of the tissue due to different cell compositions, but if you see unexpected behaviours you should consider filtering out those regions.""")
    return


@app.cell
def _(adata, sc):
    #cell area
    sc.pl.spatial(adata, color='cell_area', spot_size=10, title="Spatial Distribution of Cell Area", cmap='viridis_r')

    #nucleus area
    sc.pl.spatial(adata, color='nucleus_area', spot_size=10, title="Spatial Distribution of Nucleus Area", cmap='viridis_r')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    In the same way, we can look at the spatial distibution of one of the genes. For example the gene "Neurod6"
    """
    )
    return


@app.cell
def _(adata, sc):
    sc.pl.spatial(adata, color=['Neurod6'], spot_size=10, title="Neurod6 Gene Expression",cmap='viridis_r')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Plot distribution of total transcript counts per cell and distribution of Nucleus area after filtering""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Control probes and decoding metrics.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We can also look into the control probe counts and their spatial distribution to assess background noise and technical quality in the dataset. First, we plot the distribution of control probe counts across cells to understand how frequently control probes are detected. This distribution provides insights into potential technical noise or background signals, as higher-than-expected control counts may suggest artifacts or contamination. Next, we create a spatial plot of control probe counts, which allows us to see if any regions in the tissue exhibit unexpectedly high control counts, possibly indicating localized technical issues.""")
    return


@app.cell
def _(adata, plt, sc, sns):
    #plot distribution of control probe counts
    plt.figure(figsize=(8, 5))
    sns.histplot(adata.obs['control_probe_counts'], kde=True, bins=50)
    plt.xlabel("Control Probe Counts")
    plt.ylabel("Number of Cells")
    plt.title("Distribution of Control Probe Counts")
    plt.show()

    #spatial plot of control probe counts
    sc.pl.spatial(adata, color='control_probe_counts', spot_size=10, title="Spatial Distribution of Control Probe Counts", cmap = 'viridis_r')
    sc.pl.spatial(adata, color='unassigned_codeword_counts', spot_size=15, title="Spatial Distribution of Unassigned codeword Counts", cmap='viridis_r')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""In this case we have very low control probe counts and they seem to be scattered randomly over the tissue, so it would indicate a quite good quality.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Filtering


    ### Explore filtering criteria

    We will explore a bit further what suitable cutoffs for filtering may be, so first we look at several of the stats in paralell on the tissue. 

    Some stats that you may want to consider for filtering are:

    * Total counts per cell
    * Detected genes per cell
    * Cell area
    * Nuclei area
    """
    )
    return


@app.cell
def _(adata, sc):
    sc.pl.spatial(adata, color=['cell_area','nucleus_area'], spot_size=15, color_map = "viridis_r")
    sc.pl.spatial(adata, color=['total_counts','n_genes'], spot_size=15,  color_map = "viridis_r")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    As you can see, the top region of the section has very low counts, but still looks normal with regards to segmentation.  There are also a lot of scattered fragments of tissue outside of the secttion. 

    In such a case it would be adviced to remove those part of the section before analysis. If the region clearly deviates in quality we can remove it by simple cutoffs on number of counts and genes, if that is not possible, you may have to manually select the region to remove. For this purpose you would need an interactive browser like Napari. 

    #### Cutoff for total counts

    Here, we try a few different cutoffs on number of total counts and visualize it on the tissue. It is evident that we have to remove cells that have very few counts as they will be hard to analyse, also the cells with extremely high counts are most likely not single cells.
    """
    )
    return


@app.cell
def _(adata, sc, sns):
    sns.histplot(adata.obs['total_counts'], kde=True, bins=100)

    adata.obs['min10'] = adata.obs['total_counts'] < 10
    print(adata.obs['min10'].value_counts())
    adata.obs['min30'] = adata.obs['total_counts'] < 30
    print(adata.obs['min30'].value_counts())
    adata.obs['max600'] = adata.obs['total_counts'] > 600
    print(adata.obs['max600'].value_counts())

    sc.pl.spatial(adata, color=['min10','min30','max600'], spot_size=15, palette = ["lightgray","blue"])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    It looks like we get rid of most of the problematic regions with a cutoff of 30 counts. The really high count cells are scattered all over, but we will still filter those as well for now.


    **OBS!** You can do similar plots for all filtering criteria, if you have time.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Filter out cells and genes

    What do you think could be a good cutoff for filtering out lowly expressed genes and cells with low transcript count?
    Lets removes low-quality cells and genes from the dataset, which helps reduce noise and computational load in downstream analyses. This filtering step ensures that the dataset is focused on cells and genes with a minimum level of expression, which are more likely to be biologically relevant.

    For now, we will keep cells with these critera:

    * Nuclei size > 7
    * Cell area 20 - 1000
    * Total counts 20 - 600
    * Number of genes 4 - 150

    It is also adviced to filter out genes that are only detected in a few cells, in the code we set the cutoff to 5 cells, but in this dataset all genes are detected above the threshold.
    """
    )
    return


@app.cell
def _(adata, sc):
    initial_cells_count = adata.n_obs
    initial_genes_count = adata.n_vars
    sc.pp.filter_cells(adata, min_counts=30, inplace=True)
    sc.pp.filter_cells(adata, max_counts=600, inplace=True)
    sc.pp.filter_cells(adata, min_genes=5, inplace=True)
    sc.pp.filter_cells(adata, max_genes=150, inplace=True)
    adata_1 = adata[adata.obs['cell_area'] > 20, :]
    adata_1 = adata_1[adata_1.obs['cell_area'] < 1000, :]
    adata_1 = adata_1[adata_1.obs['nucleus_area'] > 7, :]
    filtered_cells_count = adata_1.n_obs
    filtered_genes_count = adata_1.n_vars
    print(f'Filtered {initial_cells_count - filtered_cells_count} (out of intial {initial_cells_count} cells)')
    print(f'Filtered {initial_genes_count - filtered_genes_count} (out of intial {initial_genes_count} genes)')
    return (adata_1,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Lets plot the distributions again.""")
    return


@app.cell
def _(adata_1, plt, sns):
    _fig, _axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    sns.histplot(adata_1.obs['total_counts'], kde=True, bins=50, ax=_axs[0])
    _axs[0].set_title('Distribution of Total Transcript Counts per Cell After Filtering')
    sns.histplot(adata_1.obs['nucleus_area'], kde=True, bins=50, ax=_axs[1])
    _axs[1].set_title('Distribution of Nucleus Area After Filtering')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Lets plot the top expressed genes again. It is still the same genes, but as you can see, no gene now makes up 100% of the transcripts in any cell. So we have a dataset that we can use for clustering of the data.""")
    return


@app.cell
def _(adata_1, sc):
    sc.pl.highest_expr_genes(adata_1)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""⚠️ What differences you see in the distributions of top expressed genes before and after filtering?""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Dimensionality reduction and clustering.

    We use dimensionality reduction, clustering, and visualization techniques to analyze the structure of our dataset, group similar cells, and visualize their relationships. Here, we borrow a lot of the methods that are used in scRNAseq analysis. 

    First, we apply Principal Component Analysis (PCA) to reduce the high-dimensional gene expression data into a smaller set of principal components, retaining the main patterns of variation while reducing noise. 

    Then, we create a neighbors graph, which identifies connections between cells based on their similarity in the reduced PCA space, capturing local relationships essential for clustering. Using the Leiden clustering algorithm, we group cells into clusters based on these connections, allowing us to identify groups of similar cells that may represent distinct cell types or states. 

    We then apply UMAP (Uniform Manifold Approximation and Projection), a technique that reduces the data to two dimensions for visualization, preserving both local and global structures in the data. 

    Finally, we generate a UMAP plot with cells colored by their assigned clusters, providing an overview of the dataset's structure and making it easy to spot distinct cell populations or clusters. This visualization gives us an interpretable view of the relationships within our data, helping us understand its organization before further analysis.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Using raw counts

    To see how normalization affects diemntionality reduction and clustering, first we will run PCA, clustering and UMAP without any normalization. Then we can compare the results with the normalized data later.
    """
    )
    return


@app.cell
def _(adata_1, sc):
    sc.pp.pca(adata_1)
    sc.pp.neighbors(adata_1)
    sc.tl.leiden(adata_1)
    sc.tl.umap(adata_1)
    sc.pl.umap(adata_1, color='leiden', title='UMAP raw counts')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We save the resulting umap and clustering in new slots in `.obs` and `.obsm` as they will be overwritten when we run new analyses.""")
    return


@app.cell
def _(adata_1):
    adata_1.obsm['X_umap_counts'] = adata_1.obsm['X_umap']
    adata_1.obs['leiden_counts'] = adata_1.obs['leiden']
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Normaliation 
    Normalization adjusts for cell-specific technical differences, and log transformation makes the data easier to analyze by reducing the effects of extreme values.
    This normalization and transformation make the dataset more appropriate for dimensionality reduction, clustering, and other analyses.

    Here, we will explore dimensionality reduction and clustering using 3 different methods:

    * No normalization - using the raw counts
    * Count normalization and logtransformation
    * Pearson residuals
    * Normalization by cell area
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Lognorm""")
    return


@app.cell
def _(adata_1, sc):
    adata_1.layers['raw'] = adata_1.X.copy()
    sc.pp.normalize_total(adata_1, target_sum=10000.0)
    sc.pp.log1p(adata_1)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Comparing the effect of normalization""")
    return


@app.cell
def _(adata_1, plt, sns):
    original_counts = adata_1.layers['raw'].sum(axis=1)
    normalized_counts = adata_1.X.sum(axis=1)
    plt.figure(figsize=(10, 5))
    sns.histplot(original_counts.A1, color='blue', label='Before Normalization', kde=True)
    sns.histplot(normalized_counts.A1, color='orange', label='After Normalization', kde=True)
    plt.xlabel('Total Expression per Cell')
    plt.ylabel('Number of Cells')
    plt.title('Effect of Normalization on Expression Distribution')
    plt.legend()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Now we again run pca, clustering and umap and plot again.""")
    return


@app.cell
def _(adata_1, sc):
    sc.pp.pca(adata_1)
    sc.pp.neighbors(adata_1)
    sc.tl.leiden(adata_1)
    sc.tl.umap(adata_1)
    sc.pl.umap(adata_1, color='leiden', title='UMAP lognorm')
    return


@app.cell
def _(adata_1):
    adata_1.obsm['X_umap_lognorm'] = adata_1.obsm['X_umap']
    adata_1.obs['leiden_lognorm'] = adata_1.obs['leiden']
    adata_1.layers['lognorm'] = adata_1.X.copy()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Pearson residuals normalization

    As suggested in the Lause et al 2021 https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02451-7
    """
    )
    return


@app.cell
def _(adata_1, sc):
    adata_1.X = adata_1.layers['raw'].copy()
    sc.experimental.pp.normalize_pearson_residuals(adata_1)
    return


@app.cell
def _(adata_1, sc):
    sc.pp.pca(adata_1)
    sc.pp.neighbors(adata_1)
    sc.tl.leiden(adata_1)
    sc.tl.umap(adata_1)
    sc.pl.umap(adata_1, color='leiden', title='UMAP pearson R')
    return


@app.cell
def _(adata_1):
    adata_1.obsm['X_umap_pearson'] = adata_1.obsm['X_umap']
    adata_1.obs['leiden_pearson'] = adata_1.obs['leiden']
    adata_1.layers['pearson'] = adata_1.X.copy()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Cell size normalization""")
    return


@app.cell
def _(adata_1, np):
    X = adata_1.layers['raw'].copy()
    Xnorm = np.divide(X.todense().T, adata_1.obs['cell_area'].values).T * 1000
    Xnorm = np.log1p(Xnorm)
    from scipy import sparse
    adata_1.X = sparse.csr_matrix(Xnorm.copy())
    return


@app.cell
def _(adata_1, sc):
    sc.pp.pca(adata_1)
    sc.pp.neighbors(adata_1)
    sc.tl.leiden(adata_1)
    sc.tl.umap(adata_1)
    sc.pl.umap(adata_1, color='leiden', title='UMAP size norm')
    return


@app.cell
def _(adata_1):
    adata_1.obsm['X_umap_size'] = adata_1.obsm['X_umap']
    adata_1.obs['leiden_size'] = adata_1.obs['leiden']
    adata_1.layers['size'] = adata_1.X.copy()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Now we can plot all of the different umaps together and compare them. 

    First lets plot them with their individual clusterings:
    """
    )
    return


@app.cell
def _(adata_1, plt, sc):
    _fig, _axs = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    sc.pl.embedding(adata_1, color='leiden_counts', basis='umap_counts', title='Umap counts', ax=_axs[0, 0], show=False, legend_loc='on data')
    sc.pl.embedding(adata_1, color='leiden_lognorm', basis='umap_lognorm', title='Umap lognorm', ax=_axs[0, 1], show=False, legend_loc='on data')
    sc.pl.embedding(adata_1, color='leiden_pearson', basis='umap_pearson', title='Umap pearson', ax=_axs[1, 0], show=False, legend_loc='on data')
    sc.pl.embedding(adata_1, color='leiden_size', basis='umap_size', title='Umap size norm', ax=_axs[1, 1], show=False, legend_loc='on data')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""But to see how the cells are shifting place in the umaps, we can also plot them all with the same clustering. In this case we use the clustering with lognorm.""")
    return


@app.cell
def _(adata_1, plt, sc):
    _fig, _axs = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    sc.pl.embedding(adata_1, color='leiden_lognorm', basis='umap_counts', title='Umap counts', ax=_axs[0, 0], show=False, legend_loc='on data')
    sc.pl.embedding(adata_1, color='leiden_lognorm', basis='umap_lognorm', title='Umap lognorm', ax=_axs[0, 1], show=False, legend_loc='on data')
    sc.pl.embedding(adata_1, color='leiden_lognorm', basis='umap_pearson', title='Umap pearson', ax=_axs[1, 0], show=False, legend_loc='on data')
    sc.pl.embedding(adata_1, color='leiden_lognorm', basis='umap_size', title='Umap size norm', ax=_axs[1, 1], show=False, legend_loc='on data')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Compare UMAP plots with different normalization methods. What differences do you see in the UMAP plots? What do you think is the best normalization method for this dataset?""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    As you can tell, we have more distinct clusters after normalization compared to the raw counts. And to a large degree they agree between the methods but there are some regions where the clustering may differ. 

    We can also visualize the clustering onto the sections:
    """
    )
    return


@app.cell
def _(adata_1, sc):
    sc.pl.spatial(adata_1, color='leiden_counts', spot_size=15, title='Scaled', palette=adata_1.uns['leiden_pearson_colors'])
    sc.pl.spatial(adata_1, color='leiden_lognorm', spot_size=15, title='Lognorm', palette=adata_1.uns['leiden_pearson_colors'])
    sc.pl.spatial(adata_1, color='leiden_pearson', spot_size=15, title='Pearson')
    sc.pl.spatial(adata_1, color='leiden_size', spot_size=15, title='Size', palette=adata_1.uns['leiden_pearson_colors'])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    It looks like we have better separation of cortical layers with Size or Pearson normalization.



    Using vioilin plots to see how the clusters are distributed in terms of gene expression (number of genes and the total counts) and cell size (cell and nucleus area). This helps us understand the characteristics of each cluster in terms of gene expression levels, cell size, and other features, providing insights into the biological properties of each cluster.

    It may also indicate a technical bias in the clustering if cells of similar quality are grouped together.
    """
    )
    return


@app.cell
def _(adata_1, plt, sc):
    _fig, _axs = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    sc.pl.violin(adata_1, 'total_counts', groupby='leiden_counts', ax=_axs[0, 0], show=False)
    sc.pl.violin(adata_1, 'total_counts', groupby='leiden_lognorm', ax=_axs[0, 1], show=False)
    sc.pl.violin(adata_1, 'total_counts', groupby='leiden_pearson', ax=_axs[1, 0], show=False)
    sc.pl.violin(adata_1, 'total_counts', groupby='leiden_size', ax=_axs[1, 1], show=False)
    return


@app.cell
def _(adata_1, plt, sc):
    _fig, _axs = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    sc.pl.violin(adata_1, 'cell_area', groupby='leiden_counts', ax=_axs[0, 0], show=False)
    sc.pl.violin(adata_1, 'cell_area', groupby='leiden_lognorm', ax=_axs[0, 1], show=False)
    sc.pl.violin(adata_1, 'cell_area', groupby='leiden_pearson', ax=_axs[1, 0], show=False)
    sc.pl.violin(adata_1, 'cell_area', groupby='leiden_size', ax=_axs[1, 1], show=False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Systematic Filtering of Suspected False Positives
    In Xenium data, some genes may appear to be "expressed" across large areas or in many cells due to background noise or technical artifacts, rather than true biological expression. To address this, we are implementing a gene-specific filtering approach to identify and remove suspected false positives systematically. By filtering out these spurious signals, we can focus our analysis on more reliable gene expression patterns, improving the quality of downstream analyses.

    We now have expression of genes with different normalizations, in this case we will select the log-normalized data to work with, but it is suggested to use whatever normalization that you are using for the analysis.
    """
    )
    return


@app.cell
def _(adata_1):
    adata_1.X = adata_1.layers['lognorm'].copy()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    #### 1- Calculate the Mean Expression per Gene per Cluster
    In the first step, we calculate the mean expression of each gene within each cluster. Clustering organizes cells into groups that likely share biological characteristics, and taking the average expression of each gene within clusters provides a baseline for typical expression levels. This allows us to identify clusters where a gene has unusually low expression, which might indicate noise rather than true expression.
    """
    )
    return


@app.cell
def _(adata_1):
    mean_expression = adata_1.to_df().groupby(adata_1.obs['leiden_lognorm']).mean()
    return (mean_expression,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    #### 2- Determine the Maximum Expression Level of Each Gene Across Clusters
    Then, we find the maximum expression level for each gene across all clusters. This maximum value serves as a reference point for each gene’s typical expression level in the dataset. The highest expression level of each gene is assumed to represent meaningful expression, while lower values may be more likely to represent noise or background.
    """
    )
    return


@app.cell
def _(mean_expression):
    max_expression_levels = mean_expression.max(axis=0)
    return (max_expression_levels,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    #### 3- Calculate Spurious Expression Threshold for Each Gene
    We define a threshold, here set at 5%, which will be used to identify potential false positives. This threshold means that if a gene's expression in a given cluster is below 5% of its highest expression level across all clusters, we will consider it to be a likely false positive. This step allows us to systematically identify clusters where the gene's expression is likely due to background noise rather than true biological signal. Then, we calculate a `spurious expression threshold` for each gene, which is 5% of its maximum expression level. This threshold provides a cut-off below which we consider expression values to be suspected false positives. By applying this threshold, we can filter out clusters where a gene’s expression level is unlikely to be biologically meaningful.
    """
    )
    return


@app.cell
def _(max_expression_levels):
    _threshold = 0.05
    spurious_expression_threshold = max_expression_levels * _threshold
    return (spurious_expression_threshold,)


@app.cell
def _(adata_1, mean_expression, plt, sc, sns, spurious_expression_threshold):
    for _gene in adata_1.var_names[:5]:
        cut = spurious_expression_threshold[_gene]
        _fig, _axs = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)
        sc.pl.violin(adata_1, _gene, groupby='leiden_lognorm', use_raw=False, ax=_axs[0], show=False, ylabel=str(cut))
        sns.barplot(mean_expression[_gene], ax=_axs[1])
        plt.axhline(y=cut, color='black')
        print(_gene, ':', cut)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    #### 4- Identify Suspected False Positives for Each Cluster and Gene
    Next we creat a dictionary, `suspected_false_positives`, to store suspected false positive genes for each cluster.
    For each cluster, it checks each gene’s mean expression level within that cluster. If the gene's mean expression is below the spurious expression threshold (5% of its maximum expression), the gene is flagged as a suspected false positive in that cluster.
    """
    )
    return


@app.cell
def _(adata_1, mean_expression, spurious_expression_threshold):
    suspected_false_positives = {}
    for _cluster in mean_expression.index:
        suspected_genes = []
        for _gene in adata_1.var_names:
            if mean_expression.at[_cluster, _gene] < spurious_expression_threshold[_gene]:
                suspected_genes.append(_gene)
        suspected_false_positives[_cluster] = suspected_genes
    return (suspected_false_positives,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Now we can see a quick summary of the suspected false positives across clusters""")
    return


@app.cell
def _(plt, sns, suspected_false_positives):
    for _cluster, _genes in suspected_false_positives.items():
        print(f'Cluster {_cluster} has {len(_genes)} suspected false positive genes: {_genes[:10]}')
    cluster_counts = {_cluster: len(_genes) for _cluster, _genes in suspected_false_positives.items()}
    _clusters = list(cluster_counts.keys())
    counts = list(cluster_counts.values())
    plt.figure(figsize=(12, 6))
    sns.barplot(x=_clusters, y=counts, palette='viridis')
    plt.title('Count of Suspected False Positive Genes per Cluster')
    plt.xlabel('Cluster')
    plt.ylabel('Count of Suspected False Positive Genes')
    plt.xticks(rotation=45)
    plt.show()
    return


@app.cell
def _(pd, plt, sns, suspected_false_positives):
    all_suspected_fp_genes = set()
    for _genes in suspected_false_positives.values():
        all_suspected_fp_genes.update(_genes)
    all_suspected_fp_genes_list = sorted(all_suspected_fp_genes)
    _clusters = list(suspected_false_positives.keys())
    heatmap_data = pd.DataFrame(0, index=list(all_suspected_fp_genes_list), columns=_clusters)
    for _cluster, _genes in suspected_false_positives.items():
        heatmap_data.loc[_genes, _cluster] = 1
    plt.figure(figsize=(12, 10))
    sns.heatmap(heatmap_data, cmap='viridis', cbar_kws={'label': 'Presence of Suspected False Positives'})
    plt.title('Heatmap of Suspected False Positive Genes Across Clusters')
    plt.xlabel('Cluster')
    plt.ylabel('Gene')
    plt.show()
    return


@app.cell
def _(adata_1, mean_expression, spurious_expression_threshold):
    background_noise_counts = {_gene: 0 for _gene in adata_1.var_names}
    for _gene in adata_1.var_names:
        for _cluster in mean_expression.index:
            if mean_expression.at[_cluster, _gene] < spurious_expression_threshold[_gene]:
                background_noise_counts[_gene] = background_noise_counts[_gene] + 1
    highest_noise_gene = max(background_noise_counts, key=background_noise_counts.get)
    highest_noise_count = background_noise_counts[highest_noise_gene]
    print(f"The gene with the highest background noise is '{highest_noise_gene}', detected below the threshold in {highest_noise_count} clusters.")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    #### 5- Adjusting gene expression values across the dataset (background noise reduction)
    Now, we aim to reduce background noise by adjusting gene expression levels across all cells based on suspected false positives. Specifically, if a gene has low expression (below the spurious expression threshold) in any cluster, we take the highest of those low expression levels and subtract it from all cells for that gene. This method removes unspecific, potentially spurious expression, which can help improve cluster separation and make biologically relevant expression patterns clearer.
    """
    )
    return


@app.cell
def _(adata_1, mean_expression, np, sc, spurious_expression_threshold):
    adata_app_1 = adata_1.copy()
    max_filtered_expression = {}
    adata_df = adata_app_1.to_df()
    for _gene in adata_app_1.var_names:
        below_threshold_values = []
        for _cluster in mean_expression.index:
            expression_level = mean_expression.at[_cluster, _gene]
            if expression_level < spurious_expression_threshold[_gene]:
                below_threshold_values.append(expression_level)
        if below_threshold_values:
            max_filtered_expression[_gene] = max(below_threshold_values)
        else:
            max_filtered_expression[_gene] = 0
    for _gene in adata_app_1.var_names:
        max_expr = max_filtered_expression[_gene]
        if max_expr > 0:
            adata_df[_gene] = adata_df[_gene] - max_expr
            adata_df[_gene] = np.maximum(adata_df[_gene], 0)
    adata_app_1 = sc.AnnData(X=adata_df.values, obs=adata_app_1.obs, var=adata_app_1.var, obsm=adata_app_1.obsm)
    sc.pp.pca(adata_app_1)
    sc.pp.neighbors(adata_app_1)
    sc.tl.leiden(adata_app_1)
    sc.tl.umap(adata_app_1)
    return (adata_app_1,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We can now compare umap before/after correction:""")
    return


@app.cell
def _(adata_1, adata_app_1, plt, sc):
    adata_app_1.obs['old_clust'] = adata_1.obs['leiden_lognorm']
    adata_1.obs['new_clust'] = adata_app_1.obs['leiden']
    _fig, _axs = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    sc.pl.umap(adata_app_1, color='leiden', title='UMAP after filtering, clusters after filtering', ax=_axs[0, 0], show=False, legend_loc='on data')
    sc.pl.umap(adata_app_1, color='old_clust', title='UMAP after filtering, clusters before filtering', ax=_axs[0, 1], show=False, legend_loc='on data')
    sc.pl.umap(adata_1, color='new_clust', title='UMAP before filtering, clusters after filtering', ax=_axs[1, 0], show=False, legend_loc='on data')
    sc.pl.umap(adata_1, color='leiden_lognorm', title='UMAP before filtering, clusters before filtering', ax=_axs[1, 1], show=False, legend_loc='on data')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Compare expression before/after:""")
    return


@app.cell
def _(adata_1, adata_app_1, np, plt):
    normalized_counts1 = adata_1.X.sum(axis=1)
    normalized_counts2 = adata_app_1.X.sum(axis=1)
    plt.figure(figsize=(8, 6))
    plt.scatter(np.array(normalized_counts1), normalized_counts2, alpha=0.5)
    plt.xlabel('Sum expression before filtering')
    plt.ylabel('Sum expression after filtering')
    plt.title('Expression sum after filtering')
    plt.show()
    return


@app.cell
def _(adata_1, adata_app_1, sc):
    gene_of_interest = 'Calb1'
    sc.pl.spatial(adata_1, color=[gene_of_interest], spot_size=15, title=f'{gene_of_interest} Expression Before Adjustment', color_map='viridis_r')
    sc.pl.spatial(adata_app_1, color=[gene_of_interest], spot_size=15, title=f'{gene_of_interest} Expression After Adjustment', color_map='viridis_r')
    return (gene_of_interest,)


@app.cell
def _(adata_1, adata_app_1, gene_of_interest, plt, sns):
    _original_expression = adata_1.to_df()[gene_of_interest]
    _adjusted_expression = adata_app_1.to_df()[gene_of_interest]
    plt.figure(figsize=(10, 5))
    sns.kdeplot(_original_expression, label='Original', color='blue')
    sns.kdeplot(_adjusted_expression, label='Adjusted', color='orange')
    plt.xlabel(f'Expression of {gene_of_interest}')
    plt.ylabel('Density')
    plt.title(f'Expression Distribution of {gene_of_interest} Before and After Adjustment')
    plt.legend()
    plt.show()
    return


@app.cell
def _(adata_1, adata_app_1, plt):
    mean_expression_original = adata_1.to_df().mean(axis=0)
    mean_expression_adjusted = adata_app_1.to_df().mean(axis=0)
    plt.figure(figsize=(8, 8))
    plt.scatter(mean_expression_original, mean_expression_adjusted, alpha=0.5)
    plt.plot([0, max(mean_expression_original)], [0, max(mean_expression_adjusted)], color='red', linestyle='--')
    plt.xlabel('Mean Expression (Original)')
    plt.ylabel('Mean Expression (Adjusted)')
    plt.title('Mean Gene Expression Before and After Adjustment')
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Second approach
    Check the histogram of the gene expression for every gene and if it has a multimodal expression distribution to set based on that histogram a thereshold below which we consider the expression to be background. And we could remove from all cells the expression below that therahold.
    First making helper functions performing dip test and plotting the histogram of the gene expression.
    I've used dip test to check if the distribution is multimodal or not. If the p-value is smaller than 0.05 we consider the distribution to be multimodal.
    A question here would be that should we only consider non-zero expression values for this analysis or should we consider all expression values? If consider all values we might get a multimodal distribution all for genes, however if we only consider non-zero values we might miss some genes that are suspected to have multimodal distribution. What do you think?
    """
    )
    return


@app.cell
def _(adata_1, hf):
    genes_to_analyze = adata_1.var_names[0:len(adata_1.var_names)]
    thresholds = {}
    for _gene in genes_to_analyze:
        _threshold = hf.analyze_gene_expressions(adata_1, _gene, bandwidth=0.05, plot=False, filter_zeros=False)
        if _threshold is not None:
            thresholds[_gene] = _threshold
    return (thresholds,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Sort genes by threshold values in descending order and display the top genes with the highest thresholds""")
    return


@app.cell
def _(thresholds):
    sorted_thresholds = sorted(thresholds.items(), key=lambda x: x[1], reverse=True)
    print('Top genes with the highest thresholds:')
    for _gene, threshold_value in sorted_thresholds[:10]:
        print(f'{_gene}: {threshold_value}')
    return


@app.cell
def _(adata_1, hf):
    gene_of_interest_1 = 'Calb1'
    hf.analyze_gene_expressions(adata_1, gene_of_interest_1, bandwidth=0.05, plot=True, filter_zeros=False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Now we can do the subtraction of genes with a background expression from their original values in each cell.""")
    return


@app.cell
def _(adata_1, np, sc, sp, thresholds):
    adata_app_2 = adata_1.copy()
    if sp.issparse(adata_app_2.X):
        adata_app_2.X = sp.csr_matrix(np.expm1(adata_app_2.X.toarray()))
    else:
        adata_app_2.X = np.expm1(adata_app_2.X)
    thresholds_1 = {_gene: np.expm1(thresh) for _gene, thresh in thresholds.items()}
    for _gene, _threshold in thresholds_1.items():
        if _gene in adata_app_2.var_names:
            gene_index = adata_app_2.var_names.get_loc(_gene)
            gene_data = adata_app_2[:, _gene].X
            if sp.issparse(gene_data):
                gene_data = gene_data.toarray().flatten()
            updated_gene_data = np.maximum(gene_data - _threshold, 0)
            if sp.issparse(adata_app_2.X):
                adata_app_2[:, gene_index].X = sp.csr_matrix(updated_gene_data[:, np.newaxis])
            else:
                adata_app_2[:, gene_index].X = updated_gene_data[:, np.newaxis]
    sc.pp.log1p(adata_app_2)
    return (adata_app_2,)


@app.cell
def _(adata_1, adata_app_2, sc):
    gene_of_interest_2 = 'Calb1'
    sc.pl.spatial(adata_1, color=[gene_of_interest_2], spot_size=10, title=f'{gene_of_interest_2} Expression Before Adjustment', color_map='viridis_r')
    sc.pl.spatial(adata_app_2, color=[gene_of_interest_2], spot_size=10, title=f'{gene_of_interest_2} Expression After Adjustment', color_map='viridis_r')
    return (gene_of_interest_2,)


@app.cell
def _(adata_1, adata_app_2, gene_of_interest_2, plt, sns):
    _original_expression = adata_1.to_df()[gene_of_interest_2]
    _adjusted_expression = adata_app_2.to_df()[gene_of_interest_2]
    plt.figure(figsize=(10, 5))
    sns.kdeplot(_original_expression, label='Original', color='blue')
    sns.kdeplot(_adjusted_expression, label='Adjusted', color='orange')
    plt.xlabel(f'Expression of {gene_of_interest_2}')
    plt.ylabel('Density')
    plt.title(f'Expression Distribution of {gene_of_interest_2} Before and After Adjustment')
    plt.legend()
    plt.show()
    return


@app.cell
def _(adata_app_2, sc):
    # Principal component analysis for dimension reduction
    sc.pp.pca(adata_app_2)
    sc.pp.neighbors(adata_app_2)
    sc.tl.leiden(adata_app_2)
    sc.tl.umap(adata_app_2)
    sc.pl.umap(adata_app_2, color='leiden', title='UMAP after filtering suspected false positives (approach 2)')
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
