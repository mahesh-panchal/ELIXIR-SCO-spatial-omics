import marimo

__generated_with = "0.15.2"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## ELIXIR Spatial Transcriptomics Course
    ### Practical 1b: Segmentation free cell identification using `sainsc`
    Date: 2025-01-22

    Author(s): Niklas Müller-Bötticher, Rasool Saghaleyni

    Author(s) email: niklas.mueller-boetticher@bih-charite.de, rasool.saghaleyni@scilifelab.

    ⚠️ Note: The proper environment for this notebook is `p1_segmentation_sainsc`. It can be activated by selecting the kernel in the Jupyter notebook.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Imports""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    First we are going to load all necessary packages for the analysis. 

    We will use `sainsc` for the main analysis and `scanpy` to cluster the cell-types in 
    our unsupervised analysis.
    """
    )
    return


@app.cell
def _():
    # Don't try this at home!
    # Usually you do want to notice if warnings come up!
    import warnings

    from tqdm import TqdmWarning

    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=TqdmWarning)
    return


@app.cell
def _():
    from pathlib import Path

    import pandas as pd
    import scanpy as sc
    from sainsc.io import read_Xenium
    from sainsc.utils import celltype_signatures
    return Path, pd, read_Xenium, sc


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    First we define the paths to our directory where we keep the Xenium sample that we want 
    to analyze.
    """
    )
    return


@app.cell
def _(Path):
    # TODO: adjust to the correct path
    data_path = Path("day_1/practical_1/data")
    sample_path = data_path / "Xenium_V1_FFPE_TgCRND8_17_9_months_outs"
    return (sample_path,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Before we start it is good to get a brief overview of the two main Classes in `sainsc`.

    1. [`GridCounts`](https://sainsc.readthedocs.io/page/autoapi/sainsc/GridCounts.html):
    This class holds the data as a dictionary of sparse matrices of the same 
    shape. You rarely will need to interact with it directly unless you want to, filter the
    genes or crop/mask the sample. It mostly behaves like a Python dictionary but is implemented
    in Rust. Therefore, iterating over the count matrices of each gene might be slow as the 
    data needs to be transformed every time.
    2. [`LazyKDE`](https://sainsc.readthedocs.io/page/autoapi/sainsc/LazyKDE.html): This is the class that you mostly will interact with. It contains a `GridCounts`
    instance in its `counts` attribute and otherwise exposes almost all methods necessary to
    perform the analysis.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    We will use the *transcripts.csv.gz* (or *transcripts.parquet*) to load the locations of
    all identified transcripts. The control probes from the Xenium study will be automatically 
    filtered out.

    We can furthermore specify the size of the bins we will asign the transcripts into (by 
    default this is set to 0.5 um) and the number of threads we want to use to process the data.

    There are options to directly load data from common file formats/technologies
    such as Stereo-seq, Xenium, and Vizgen. If none of the options fit for your use case 
    you can have a look at [`LazyKDE.from_dataframe`](https://sainsc.rtd.io/api/) 
    or [`GridCounts.from_dataframe`](https://sainsc.rtd.io/api) methods.
    """
    )
    return


@app.cell
def _(read_Xenium, sample_path):
    brain = read_Xenium(sample_path / "transcripts.csv.gz", n_threads=8)

    brain
    return (brain,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    The `LazyKDE` object will give us some useful information when we print it; including
    the number of genes, the size of the sample in pixels and the resolution.

    Next we can get a quick overview of our sample by calculating the total mRNA and plotting
    it. If we can squint our eyes, we can notice some technical artifacts; the mRNA seems to
    be lower at certain locations that seem to form a grid, likely along the stitching borders.
    """
    )
    return


@app.cell
def _(brain):
    brain.calculate_total_mRNA()
    _ = brain.plot_genecount(im_kwargs={"vmax": 2})
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    We can crop our sample to remove some "dead" space to further speed up processing or to
    "zoom" into a smaller region of interest (ROI).

    Alternatively, we could also use the 
    [`GridCounts.filter_mask`](https://sainsc.readthedocs.io/page/autoapi/sainsc/GridCounts.filter_mask.html)
    method to use an arbitrary binary mask to determine the ROI. All the transcripts outside 
    our ROI will then be dropped. This allows us to filter the ROI to any shape desired.
    """
    )
    return


@app.cell
def _(brain):
    brain.counts.crop((500, None), (None, 10_000))

    brain.calculate_total_mRNA()
    _ = brain.plot_genecount(im_kwargs={"vmax": 2})
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""It is always a good idea to also check the distribution of transcripts detected per gene."""
    )
    return


@app.cell
def _(brain):
    _ = brain.plot_genecount_histogram()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    The kernel defines on how we will smooth the gene expression. The go-to choice is to use
    a gaussian kernel, however any square `numpy.ndarray` can be used.

    The size of the kernel can either be defined in pixels or in µm (if the resolution is set). 

    The required kernel size may depend on the technology. Here, we will use 2.5 µm.
    """
    )
    return


@app.cell
def _(brain):
    brain.gaussian_kernel(2.5, unit="um")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Now we can first smooth the total mRNA and visualize it.""")
    return


@app.cell
def _(brain):
    brain.calculate_total_mRNA_KDE()
    _ = brain.plot_KDE()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    The distribtuion of the smoothed gene expression can be used to determine a threshold to
    use for filtering out background noise. Here, a value of ~ 0.02 seems to be a good 
    first choice.
    """
    )
    return


@app.cell
def _(brain):
    _ = brain.plot_KDE_histogram(bins=200)
    return


@app.cell
def _(brain):
    _ = brain.plot_KDE_histogram(bins=100, range=(0, 0.1))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Even though in the images above it looked like the background was empty, adjusting the 
    color scale will quickly prove us that this is not the case.

    Here, masking the ROI could be used to completely remove those counts.
    Filtering the background, on the other side, will only affect the visualization 
    but not the processing.
    """
    )
    return


@app.cell
def _(brain):
    _ = brain.plot_KDE(im_kwargs={"vmax": 0.02})
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Now we can filter the background. Note, that later we can further refine this and define
    background filter on the total mRNA KDE per cell type.
    """
    )
    return


@app.cell
def _(brain):
    brain.filter_background(0.02)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Unsupervised analysis""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    To generate the cell-type map we will need a set of gene expression signatures.
    These can either be derived from previous studies e.g. scRNAseq or we can identify them 
    *de novo* from the sample we are analysing.

    The *de novo* approach works by finding the local maxima of the gene expression and 
    treating these as proxies for cells. We can then use standard single-cell/spatial workflows 
    to process and cluster the cells. The cell-types indentified in the clustering can then 
    be used to calculate the gene expression signatures.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    The first step is to identify the local maxiam, we set a minimum distance to avoid 
    sampling too many close-by spots.
    """
    )
    return


@app.cell
def _(brain):
    brain.find_local_maxima(5)
    brain
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Find cell-type signatures""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Next we can load the local maxima into and `AnnData` object and then proceed to identify
    clusters.
    """
    )
    return


@app.cell
def _(brain):
    local_max = brain.load_local_maxima()
    local_max
    return (local_max,)


@app.cell
def _():
    # for reproducibility
    random_state = 42
    return (random_state,)


@app.cell
def _(local_max, sc):
    sc.pp.normalize_total(local_max)
    return


@app.cell
def _(local_max):
    local_max.layers["counts"] = local_max.X.copy()
    return


@app.cell
def _(local_max, random_state, sc):
    sc.pp.log1p(local_max)
    sc.pp.pca(local_max, random_state=random_state)
    sc.pp.neighbors(local_max, random_state=random_state)
    sc.tl.umap(local_max, random_state=random_state)
    return


@app.cell
def _(local_max, random_state, sc):
    sc.tl.leiden(
        local_max, resolution=2, flavor="igraph", random_state=random_state
    )

    sc.pl.umap(local_max, color="leiden")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    After we have identified our cell-types/clusters we can easily calculate the gene expression
    signatures.

    Note, the gene expression signatures should be strictly positive i.e. they should not be 
    calculated from data that has been standardized or similar.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""The following code cell can be ignored (you still need to run it but do not pay to much attention), it is just a currently not yet released improvement to reduce the memory usage when calculating the cell-type signatures."""
    )
    return


@app.cell
def _(pd):
    from collections.abc import Hashable
    import anndata as ad
    import numpy as np
    from numpy.typing import DTypeLike


    def celltype_signatures_1(
        adata: ad.AnnData,
        *,
        celltype_col: str = "leiden",
        layer: str | None = None,
        dtype: DTypeLike = np.float32,
    ) -> pd.DataFrame:
        """
        Calculate gene expression signatures per 'cell type'.

        Parameters
        ----------
        adata : anndata.AnnData
        celltype_col : str, optional
            Name of column in :py:attr:`anndata.AnnData.obs` containing cell-type
            information.
        layer : str, optional
            Which :py:attr:`anndata.AnnData.layers` to use for aggregation. If `None`,
            :py:attr:`anndata.AnnData.X` is used.
        dytpe : numpy.typing.DTypeLike
            Data type to use for the signatures.

        Returns
        -------
        pandas.DataFrame
            :py:class:`pandas.DataFrame` of gene expression aggregated per 'cell type'.
        """
        X = adata.X if layer is None else adata.layers[layer]
        grouping = adata.obs.groupby(
            celltype_col, observed=True, sort=False
        ).indices
        signatures: dict[Hashable, np.ndarray] = {}
        for name, indices in grouping.items():
            mean_X_group = X[indices].mean(axis=0, dtype=dtype)
            signatures[name] = (
                mean_X_group.A1
                if isinstance(mean_X_group, np.matrix)
                else mean_X_group
            )
        return pd.DataFrame(signatures, index=adata.var_names)
    return (celltype_signatures_1,)


@app.cell
def _(celltype_signatures_1, local_max):
    signatures = celltype_signatures_1(local_max, celltype_col="leiden")
    return (signatures,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Generate cell-type map""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    To generate the cell-type map we just need to pass the signature DataFrame to the `assign_celltype` method.

    If the gene expression varies across multiple orders of magnitude across genes it might be useful to
    use log-transformation after calculating the KDE. In this case the gene expression 
    signatures should be calculated from log-transformed data, as well.
    """
    )
    return


@app.cell
def _(brain, signatures):
    brain.assign_celltype(signatures, log=True)
    return


@app.cell
def _(local_max):
    # maintain the same coloring as in UMAP
    cmap = {
        cluster: color
        for cluster, color in zip(
            local_max.obs["leiden"].cat.categories, local_max.uns["leiden_colors"]
        )
    }
    return (cmap,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We can now visualize our cell-type map.""")
    return


@app.cell
def _(brain, cmap):
    _ = brain.plot_celltype_map(cmap=cmap)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    The assignment score can be helpful to identify regions with low confidence in the cell-type
    assignment. This is especially useful when using pre-existing cell-type signatures as it 
    might highlight regions where we couldn't map any cell-type with high confidence and therefore
    might indicate that cell-types are missing in the reference.
    """
    )
    return


@app.cell
def _(brain):
    _ = brain.plot_assignment_score(remove_background=True)
    return


@app.cell
def _(brain, pd):
    kde_per_celltype = pd.DataFrame(
        {
            "kde": brain.total_mRNA_KDE.flatten(),
            "celltype": pd.Categorical.from_codes(
                brain.celltype_map.flatten(), categories=brain.celltypes
            ),
        }
    ).dropna()

    celltype_threshold = (
        kde_per_celltype.groupby("celltype", observed=True).quantile(0.5)["kde"]
        / 2
    ).to_dict()

    min_t = 0.02

    celltype_threshold = {
        ct: (t if t > min_t else min_t) for ct, t in celltype_threshold.items()
    }

    brain.filter_background(celltype_threshold)
    return


@app.cell
def _(brain, cmap):
    _ = brain.plot_celltype_map(cmap=cmap)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We can zoom-in by defining the ROI that we want to plot.""")
    return


@app.cell
def _(brain, cmap):
    roi = ((1_000, 4_500), (4_000, 6_000))

    _ = brain.plot_celltype_map(
        cmap=cmap, crop=roi, scalebar_kwargs={"box_alpha": 0.7}
    )
    return (roi,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""We can also highlight only one/few cell-types by removing the rest from the colormap."""
    )
    return


@app.cell
def _(brain, roi):
    cmap2 = {"0": "yellow"}

    _ = brain.plot_celltype_map(
        cmap=cmap2, crop=roi, scalebar_kwargs={"box_alpha": 0.7}
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""Visualizing the gene expression can help rationalizing the assigned cell-types."""
    )
    return


@app.cell
def _(brain, roi, signatures):
    _ = brain.plot_KDE(gene=signatures["0"].idxmax(), crop=roi)
    return


@app.cell
def _(brain, roi, signatures):
    _ = brain.plot_KDE(gene=signatures["18"].idxmax(), crop=roi)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Bonus task: Supervised analysis""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""Try using `sainsc` for a supervised analysis leveraging the cell-type signatures obtained from your previous segmentation-based analysis workflow."""
    )
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
