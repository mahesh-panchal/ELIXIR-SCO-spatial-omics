import marimo

__generated_with = "0.15.2"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## ELIXIR Spatial Transcriptomics Course
    ### Practical 1c: Segmentation of spatial transcriptomics data using `cellpose`
    Date: 2025-01-22

    Author(s): Rasool Saghaleyni

    Author(s) email: rasool.saghaleyni@scilifelab.se

    ⚠️ Note: The proper environment for this notebook is `p1_segmentation_cellpose`. It can be activated by selecting the kernel in the Jupyter notebook.
    """
    )
    return


@app.cell
def _():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from cellpose import models, io, plot
    from tifffile import imread
    import os
    import tifffile
    from cellpose import plot
    import shapely.geometry as geometry
    from shapely.geometry import Polygon
    from shapely.affinity import translate, scale
    from shapely.errors import TopologicalError
    from rasterio import features
    from sklearn.metrics import jaccard_score
    from skimage.measure import regionprops_table


    # Set up plotting aesthetics
    sns.set(style='whitegrid')
    # '%matplotlib inline' command supported automatically in marimo
    return (
        Polygon,
        TopologicalError,
        features,
        jaccard_score,
        models,
        np,
        pd,
        plot,
        plt,
        regionprops_table,
        scale,
        sns,
        tifffile,
        translate,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Loading and inspecting the morphology image is a foundational step in the analysis. Understanding the image's dimensions helps guide downstream processing, including selecting channels or slices for segmentation and defining an ROI if needed.""")
    return


@app.cell
def _(tifffile):
    image_data_path = 'day_1/practical_1/data/Xenium_V1_FFPE_TgCRND8_17_9_months_outs/morphology.ome.tif'
    image = tifffile.imread(image_data_path)
    print(f"Image shape: {image.shape}")
    return image, image_data_path


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    we retrieve and inspect the metadata embedded within the morphology image file. Opening the file with `tifffile.TiffFile` allows us to access not only the pixel data but also any associated metadata stored in OME (Open Microscopy Environment) format. By using a context manager (`with` statement), we ensure that the file is properly handled, meaning it will close automatically once we’re done, helping to avoid potential file-handling errors.

    The OME metadata contains crucial information about the image acquisition settings. Extracting it with `ome_metadata = tif.ome_metadata` provides us with details about the microscope settings, pixel size, and other experimental parameters. This metadata appears in XML format, which is printed for review. Examining this data is essential to understand the spatial resolution of the image, enabling us to relate image coordinates to real-world units, such as micrometers. Knowing these specifics is key for aligning segmentation results accurately with the spatial features observed in the morphology image.
    """
    )
    return


@app.cell
def _(image_data_path, tifffile):
    with tifffile.TiffFile(image_data_path) as tif:
        ome_metadata = tif.ome_metadata
        print(ome_metadata)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Lets preview the transcripts data again:""")
    return


@app.cell
def _(pd):
    transcriptomics_data_path = 'day_1/practical_1/data/Xenium_V1_FFPE_TgCRND8_17_9_months_outs/transcripts.csv.gz'
    data = pd.read_csv(transcriptomics_data_path, compression='gzip')
    print(data.head())
    return (data,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Here, we are checking and handling the dimensionality of the morphology image to extract a usable 2D channel for segmentation and analysis. To start, `print(f"Image dimensions: {image.ndim}")` reveals the number of dimensions in the image array. Images acquired from microscopy can have multiple dimensions, often representing different z-slices, time points, or channels (e.g., specific fluorescent stains). Knowing the exact number of dimensions is essential for understanding the structure of the data and selecting the specific layer or channel needed for downstream tasks.

    We then use conditional statements to select the appropriate 2D plane. If the image has five dimensions—typically representing time, Z (depth), channels, height, and width—we select the first time point, Z-slice, and channel to reduce it to 2D. Similarly, for four-dimensional images (likely Z, channels, height, and width), we choose the first Z-slice and channel. In the case of three-dimensional images, we assume they represent channels, height, and width, and extract the first channel. Finally, if the image is already 2D, we simply assign it to image_channel without further modification.

    This step ensures that we have a consistent, interpretable 2D array (image_channel) for the following analysis. By isolating a single plane or channel, we simplify the data, making it easier to overlay segmentations or spatial features without the added complexity of multiple dimensions. This also ensures that our chosen channel represents the tissue morphology effectively
    """
    )
    return


@app.cell
def _(image):
    print(f"Image dimensions: {image.ndim}")
    if image.ndim == 5:
        # Example shape: (Time, Z, Channels, Height, Width)
        # Select the first time point, z-slice, and channel
        image_channel = image[0, 0, 0, :, :]
    elif image.ndim == 4:
        # Example shape: (Z, Channels, Height, Width)
        image_channel = image[0, 0, :, :]
    elif image.ndim == 3:
        # Example shape: (Channels, Height, Width)
        image_channel = image[0, :, :]
    else:
        # Already a 2D image
        image_channel = image
    return


@app.cell
def _(image, np):
    image_channel_1 = np.max(image, axis=0)
    image_channel_1 = image_channel_1.astype(np.uint16)
    return (image_channel_1,)


@app.cell
def _(image_channel_1, plt):
    plt.figure(figsize=(8, 8))
    plt.imshow(image_channel_1)
    plt.title('Selected Image for Segmentation')
    plt.axis('off')
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Since the original image is too big here we define a region of interest (ROI) within the larger image, focusing on a smaller area for more efficient and targeted analysis. This is particularly useful for high-resolution images where analyzing the entire field of view might be computationally intensive.

    We begin by setting a scale_factor, which controls the size of the ROI as a fraction of the full image dimensions. Here, `scale_factor = 0.05` means that the ROI will cover 5% of the original image's width and height. Adjusting this factor allows flexibility in focusing on larger or smaller portions of the image, depending on the needs of the analysis.

    Using `image_channel.shape`, we extract the height and width of the full image. Then, by multiplying these dimensions by `scale_factor`, we calculate the width and height of the ROI (`roi_width` and `roi_height`). Converting these values to integers ensures that they’re compatible with image indexing.

    Finally, we make an optional adjustment to ensure that the ROI dimensions are even numbers, which can simplify image processing tasks. We do this by reducing `roi_width` and `roi_height` by 1 if they are odd, using modulo operations. This adjustment helps avoid issues when working with certain algorithms that may require even-numbered dimensions, ensuring that the ROI dimensions are compatible with a range of image processing techniques.
    """
    )
    return


@app.cell
def _(image_channel_1):
    _scale_factor = 0.05
    _image_height, _image_width = image_channel_1.shape
    _roi_width = int(_image_width * _scale_factor)
    _roi_height = int(_image_height * _scale_factor)
    _roi_width = _roi_width - _roi_width % 2
    _roi_height = _roi_height - _roi_height % 2
    _x_center = _image_width // 2
    _y_center = _image_height // 2
    x_start = _x_center - _roi_width // 2
    x_end = _x_center + _roi_width // 2
    y_start = _y_center - _roi_height // 2
    y_end = _y_center + _roi_height // 2
    return x_end, x_start, y_end, y_start


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Extract the ROI from the image""")
    return


@app.cell
def _(image_channel_1, x_end, x_start, y_end, y_start):
    roi_image = image_channel_1[y_start:y_end, x_start:x_end]
    print(f'ROI image shape: {roi_image.shape}')
    return (roi_image,)


@app.cell
def _(plt, roi_image):
    plt.figure(figsize=(8, 8))
    plt.imshow(roi_image)
    plt.title('ROI Image for Segmentation')
    plt.axis('off')
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Next we map the ROI pixel coordinates back to real-world units (micrometers) and filter the spatial transcriptomics data to include only the transcripts within the ROI. This enables precise alignment of the transcript data with the selected region in the image.

    We start by defining scaling factors for converting image pixels to micrometers. Here, `x_scale` and `y_scale` represent the conversion rate based on the pixel size provided in the image metadata: each micrometer contains approximately 4.7 pixels (1 / 0.2125). This conversion allows us to translate pixel coordinates into micrometer units, which are required for comparing and aligning data across different scales.

    Using these scaling factors, we calculate the boundaries of the ROI in micrometers. For each dimension, `x_start`, `x_end`, `y_start`, and `y_end` (which are pixel coordinates from the original image), we divide by the scaling factor to obtain the corresponding boundaries in micrometers: `x_start_um`, `x_end_um`, `y_start_um`, and `y_end_um`. This step ensures that our ROI is defined consistently in both pixel and physical units.

    Next, we filter the transcriptomics data to include only the transcripts located within the ROI. We use conditional filtering on the `x_location` and `y_location` columns of the data DataFrame, retaining only the transcripts whose coordinates fall within the calculated micrometer boundaries. The result is stored in roi_data, which represents the subset of transcripts that reside within our chosen ROI.

    Finally, by printing `len(roi_data)`, we get a quick count of the transcripts within the ROI.
    """
    )
    return


@app.cell
def _(data, x_end, x_start, y_end, y_start):
    #using the scaling factors from before
    x_scale = 1 / 0.2125  # pixels per µm
    y_scale = 1 / 0.2125  # pixels per µm

    #roi boundaries in micrometers
    x_start_um = x_start / x_scale
    x_end_um = x_end / x_scale
    y_start_um = y_start / y_scale
    y_end_um = y_end / y_scale

    #filter transcripts within the ROI boundaries
    roi_data = data[
        (data['x_location'] >= x_start_um) &
        (data['x_location'] < x_end_um) &
        (data['y_location'] >= y_start_um) &
        (data['y_location'] < y_end_um)
    ].copy()
    print(f"Number of transcripts in ROI: {len(roi_data)}")
    return roi_data, x_scale, x_start_um, y_scale, y_start_um


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Now we should subtract `x_start_um` from each transcript’s `x_location` and `y_start_um` from each `y_location` in `roi_data`. By doing so, we create new columns, `x_location_roi` and `y_location_roi`, that represent each transcript’s position relative to the top-left corner of the ROI rather than the full image.""")
    return


@app.cell
def _(roi_data, x_start_um, y_start_um):
    roi_data['x_location_roi'] = roi_data['x_location'] - x_start_um
    roi_data['y_location_roi'] = roi_data['y_location'] - y_start_um
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Now we finally set up the Cellpose segmentation model to identify cells within the ROI. Cellpose is a versatile deep learning-based tool commonly used for cell segmentation, especially on fluorescence and cytoplasmic images. Here, we are preparing the model for use in the analysis.

    First, we import the models module from the cellpose package, which provides access to pre-trained Cellpose models. Next, we initialize a model instance using `models.Cellpose()`. By setting `gpu=False`, we specify that the model will run on the CPU. This is useful if GPU resources are unavailable, though using a GPU can speed up the segmentation process if it is an option.

    We also set `model_type='cyto'`, indicating that the model should use Cellpose’s pre-trained “cyto” (cytoplasm) model, which is optimized for identifying cell boundaries in images with visible cell structures. This choice is typically well-suited for images showing cell cytoplasm, though Cellpose offers other model types, like “nuclei,” if our focus were solely on nuclear segmentation.
    """
    )
    return


@app.cell
def _(models):
    model = models.CellposeModel(gpu=False)
    return (model,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Then, we estimate the average cell diameter in pixels and then use Cellpose to perform cell segmentation on the ROI.

    We start by setting `cell_diameter_um` to an estimated cell diameter in micrometers, which is based on biological knowledge of cell sizes in the specific tissue or sample type. Here, we use 10 micrometers as the mouse cell daimeter for mouse brains cells is estrimated 7-10 micrometers, but this value can be adjusted based on the specific dataset.

    To convert this estimate into pixel units, we multiply `cell_diameter_um` by the scaling factor `x_scale` (pixels per micrometer), calculated previously. This results in `cell_diameter_pixels`, an approximation of the cell diameter in the pixel space of the image. By converting the diameter to pixels, we ensure that Cellpose can interpret the size parameter relative to the image’s resolution.

    Next, we run the Cellpose model on `roi_image`, the 2D image extracted from the ROI. The model’s eval() function applies the segmentation model to the image, using `diameter=cell_diameter_pixels` to guide the segmentation scale. The channels=[0, 0] parameter specifies that the image is grayscale; both entries as 0 indicate that there is a single channel for both input and detection purposes.

    The eval() function returns several outputs:

    - `masks`: a labeled mask array where each detected cell has a unique identifier,

    - `flows`: which provides information about cell boundary flows,

    - `styles`: representing style vectors for detected objects, and

    - `diams`: the diameter used in the model (helpful if it has been automatically adjusted).
    """
    )
    return


@app.cell
def _(model, roi_image, x_scale):
    cell_diameter_um = 10  # µm
    cell_diameter_pixels = cell_diameter_um * x_scale
    print(f"Estimated cell diameter in pixels: {cell_diameter_pixels}")
    #run segmentation
    masks, flows, styles = model.eval(
        roi_image,
        diameter=cell_diameter_pixels
    )
    return flows, masks


@app.cell
def _(masks):
    print(f"Number of cells detected in ROI: {masks.max()}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Now we convert the transcript coordinates within the ROI from micrometers to pixel indices, preparing them for alignment with the segmentation mask.

    First, we extract the x- and y-coordinates in micrometers from `roi_data`, which represents transcript locations relative to the top-left corner of the ROI. These coordinates are stored in `x_coords_um` and `y_coords_um`, making it easy to work directly with arrays of positions.

    To map these positions into the pixel space of `roi_image`, we multiply each coordinate by the scaling factor (`x_scale` and `y_scale`) previously defined. This scaling factor converts micrometers into pixel units, allowing us to obtain `x_indices` and `y_indices` the pixel indices that match the resolution of the segmentation mask.

    By converting coordinates to pixel indices, we can precisely locate each transcript in the context of the segmented cells within the ROI.
    """
    )
    return


@app.cell
def _(roi_data, x_scale, y_scale):
    x_coords_um = roi_data['x_location_roi'].values
    y_coords_um = roi_data['y_location_roi'].values

    #pixel indices
    x_indices = (x_coords_um * x_scale).astype(int)
    y_indices = (y_coords_um * y_scale).astype(int)
    return x_indices, y_indices


@app.cell
def _(np, roi_image, x_indices, y_indices):
    _roi_height, _roi_width = roi_image.shape
    x_indices_1 = np.clip(x_indices, 0, _roi_width - 1)
    y_indices_1 = np.clip(y_indices, 0, _roi_height - 1)
    return x_indices_1, y_indices_1


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Next we assign each transcript to a segmented cell based on its pixel coordinates, linking gene expression data to specific cells within the ROI.""")
    return


@app.cell
def _(masks, roi_data, x_indices_1, y_indices_1):
    _cell_labels = masks[y_indices_1, x_indices_1]
    roi_data['cellpose_cell_id'] = _cell_labels
    print(roi_data.head())
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Keep only transcripts assigned to a cell""")
    return


@app.cell
def _(roi_data):
    assigned_data = roi_data[roi_data['cellpose_cell_id'] > 0].copy()
    print(f"Total transcripts in ROI: {len(roi_data)}")
    print(f"Assigned transcripts: {len(assigned_data)}")
    return (assigned_data,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""we visualize the results of the Cellpose segmentation overlayed on the ROI image, allowing us to inspect how well the model identified individual cells in the selected region.""")
    return


@app.cell
def _(flows, masks, plot, plt, roi_image):
    fig = plt.figure(figsize=(8, 8))
    plot.show_segmentation(fig, roi_image, masks, flows[0])
    plt.title('Cellpose Segmentation on ROI')
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Group by cell and gene to get expression counts""")
    return


@app.cell
def _(assigned_data):
    expression_per_cell = assigned_data.groupby(['cellpose_cell_id', 'feature_name']).size().reset_index(name='count')
    expression_matrix = expression_per_cell.pivot(index='cellpose_cell_id', columns='feature_name', values='count').fillna(0)
    print(expression_matrix.head())
    return (expression_matrix,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""we can plot the locations of transcripts overlaid on the ROI image, specifically highlighting those that have been assigned to segmented cells. This helps us see how transcript data aligns with the detected cell boundaries within the region.""")
    return


@app.cell
def _(plt, roi_data, roi_image, x_indices_1, y_indices_1):
    plt.figure(figsize=(8, 8))
    plt.imshow(roi_image, cmap='gray')
    plt.scatter(x_indices_1[roi_data['cellpose_cell_id'] > 0], y_indices_1[roi_data['cellpose_cell_id'] > 0], c='red', s=5, label='Transcripts')
    plt.title('Transcripts Mapped to Segmented Cells')
    plt.axis('off')
    plt.legend()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Finally, we can visualize the expression of a specific gene across the detected cells, providing insights into the spatial distribution of gene expression within the ROI. This visualization can reveal patterns of gene expression, such as high expression in specific cell types or regions, helping to interpret the biological significance of the data.""")
    return


@app.cell
def _(expression_matrix, masks, np, plt, roi_image):
    gene_of_interest = 'Cst3'
    if gene_of_interest in expression_matrix.columns:
        cell_ids = expression_matrix.index.values
        expression_values = expression_matrix[gene_of_interest].values
        from skimage.measure import regionprops
        properties = regionprops(masks)
        centroids = np.array([prop.centroid for prop in properties])
        _cell_labels = np.array([prop.label for prop in properties])
        centroid_dict = {label: centroid for label, centroid in zip(_cell_labels, centroids)}
        cell_centroids = np.array([centroid_dict.get(cell_id, (np.nan, np.nan)) for cell_id in cell_ids])
        plt.figure(figsize=(8, 8))
        plt.imshow(roi_image, cmap='gray')
        plt.scatter(cell_centroids[:, 1], cell_centroids[:, 0], c=expression_values, cmap='viridis', s=50, edgecolors='k', label=f'Expression of {gene_of_interest}')
        plt.title(f'Expression of {gene_of_interest}')
        plt.axis('off')
        plt.colorbar(label='Expression Level')
        plt.show()
    else:
        print(f'{gene_of_interest} not found in expression matrix.')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""now we want to compare the cellpose segmentation with 10x segmentation, we can use the Jaccard index to quantify the similarity between the two segmentation masks. The Jaccard index, also known as the intersection-over-union (IoU), measures the overlap between two sets by dividing the size of their intersection by the size of their union. In the context of segmentation masks, the Jaccard index provides a measure of how well two masks align, with values closer to 1 indicating greater similarity.""")
    return


@app.cell
def _(pd):
    cells_data_path = 'day_1/practical_1/data/Xenium_V1_FFPE_TgCRND8_17_9_months_outs/cells.csv'
    cells_data = pd.read_csv(cells_data_path)
    print(cells_data.head())
    return (cells_data,)


@app.cell
def _(pd):
    nucleus_boundaries_path = 'day_1/practical_1/data/Xenium_V1_FFPE_TgCRND8_17_9_months_outs/nucleus_boundaries.csv.gz'  
    nucleus_boundaries = pd.read_csv(nucleus_boundaries_path)
    print(nucleus_boundaries.head())
    return (nucleus_boundaries,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""First we need to extract the region of intrest from the 10x segmentation mask, then we need to resize the 10x segmentation mask to the same size as the cellpose segmentation mask, then we can calculate the Jaccard index between the two masks.""")
    return


@app.cell
def _(image_channel_1):
    _scale_factor = 0.05
    _image_height, _image_width = image_channel_1.shape
    _roi_width = int(_image_width * _scale_factor)
    _roi_height = int(_image_height * _scale_factor)
    _roi_width = _roi_width - _roi_width % 2
    _roi_height = _roi_height - _roi_height % 2
    _x_center = _image_width // 2
    _y_center = _image_height // 2
    x_start_1 = _x_center - _roi_width // 2
    x_end_1 = _x_center + _roi_width // 2
    y_start_1 = _y_center - _roi_height // 2
    y_end_1 = _y_center + _roi_height // 2
    x_scale_1 = 1 / 0.2125
    y_scale_1 = 1 / 0.2125
    x_start_um_1 = x_start_1 / x_scale_1
    x_end_um_1 = x_end_1 / x_scale_1
    y_start_um_1 = y_start_1 / y_scale_1
    y_end_um_1 = y_end_1 / y_scale_1
    return (
        x_end_um_1,
        x_scale_1,
        x_start_um_1,
        y_end_um_1,
        y_scale_1,
        y_start_um_1,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Filter cells whose centroids are within the ROI""")
    return


@app.cell
def _(cells_data, x_end_um_1, x_start_um_1, y_end_um_1, y_start_um_1):
    cells_in_roi = cells_data[(cells_data['x_centroid'] >= x_start_um_1) & (cells_data['x_centroid'] < x_end_um_1) & (cells_data['y_centroid'] >= y_start_um_1) & (cells_data['y_centroid'] < y_end_um_1)].copy()
    print(f'Number of cells in ROI from original segmentation: {len(cells_in_roi)}')
    return (cells_in_roi,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Filter cell boundaries for cells in ROI""")
    return


@app.cell
def _(Polygon, TopologicalError, cells_in_roi, nucleus_boundaries):
    nucleus_boundaries_in_roi = nucleus_boundaries[nucleus_boundaries['cell_id'].isin(cells_in_roi['cell_id'])].copy()
    cell_polygons = {}
    for cell_id, group in nucleus_boundaries_in_roi.groupby('cell_id'):
        x_coords = group['vertex_x'].values
        y_coords = group['vertex_y'].values
        coords = list(zip(x_coords, y_coords))
        try:
            polygon = Polygon(coords)
            if not polygon.is_valid:
                # Attempt to fix invalid polygons
                polygon = polygon.buffer(0)
            cell_polygons[cell_id] = polygon
        except TopologicalError as e:
            print(f"Could not create polygon for cell {cell_id}: {e}")
    return (cell_polygons,)


@app.cell
def _(
    cell_polygons,
    scale,
    translate,
    x_scale_1,
    x_start_um_1,
    y_scale_1,
    y_start_um_1,
):
    def geometry_to_pixel_coords(geometry):
        geometry_shifted = translate(geometry, xoff=-x_start_um_1, yoff=-y_start_um_1)
        geometry_scaled = scale(geometry_shifted, xfact=x_scale_1, yfact=y_scale_1, origin=(0, 0))
        return geometry_scaled
    cell_polygons_px = {cell_id: geometry_to_pixel_coords(geom) for cell_id, geom in cell_polygons.items()}
    return (cell_polygons_px,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Map cell_id strings to integer labels""")
    return


@app.cell
def _(cell_polygons_px):
    cell_id_to_label = {cell_id: idx+1 for idx, cell_id in enumerate(cell_polygons_px.keys())}
    label_to_cell_id = {idx+1: cell_id for idx, cell_id in enumerate(cell_polygons_px.keys())}
    return (cell_id_to_label,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Prepare shapes for rasterization""")
    return


@app.cell
def _(cell_id_to_label, cell_polygons_px):
    shapes = [
        (geom, cell_id_to_label[cell_id])
        for cell_id, geom in cell_polygons_px.items()
    ]
    return (shapes,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Create an empty mask and rasterize the shapes""")
    return


@app.cell
def _(features, np, roi_image, shapes):
    original_masks = np.zeros_like(roi_image, dtype=np.uint16)
    original_masks = features.rasterize(
        shapes,
        out_shape=original_masks.shape,
        fill=0,
        all_touched=True,
        dtype=np.uint16
    )
    return (original_masks,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Lets compare the cellpose segmentation with 10x segmentation side by side""")
    return


@app.cell
def _(masks, original_masks, plt, roi_image):
    plt.figure(figsize=(16, 8))

    #Cellpose 
    plt.subplot(1, 2, 1)
    plt.imshow(roi_image, cmap='gray')
    plt.imshow(masks, alpha=0.5, cmap='jet')
    plt.title('Cellpose Segmentation')
    plt.axis('off')

    #Original 
    plt.subplot(1, 2, 2)
    plt.imshow(roi_image, cmap='gray')
    plt.imshow(original_masks, alpha=0.5, cmap='jet')
    plt.title('Original 10x Segmentation')
    plt.axis('off')

    plt.tight_layout()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""And we can overlay the cellpose segmentation on the 10x segmentation to see how well they align""")
    return


@app.cell
def _(masks, original_masks, plt, roi_image):
    # now we overlay both masks
    plt.figure(figsize=(8, 8))
    plt.imshow(roi_image, cmap='gray')
    plt.imshow((original_masks > 0).astype(int), cmap='Blues', alpha=0.5, label='Original')
    plt.imshow((masks > 0).astype(int), cmap='Reds', alpha=0.5, label='Cellpose')
    plt.title('Overlay of Segmentations')
    plt.axis('off')
    plt.show()
    return


@app.cell
def _(masks, original_masks):
    #convert masks to binary masks (cells vs background)
    cellpose_mask_binary = (masks > 0).astype(int)
    original_mask_binary = (original_masks > 0).astype(int)

    #flatten the masks for metric computation
    cellpose_mask_flat = cellpose_mask_binary.flatten()
    original_mask_flat = original_mask_binary.flatten()
    return cellpose_mask_flat, original_mask_flat


@app.cell
def _(cellpose_mask_flat, jaccard_score, original_mask_flat):
    jaccard = jaccard_score(original_mask_flat, cellpose_mask_flat)
    print(f'Jaccard Index: {jaccard:.4f}')
    return


@app.cell
def _(cellpose_mask_flat, np, original_mask_flat):
    def dice_coefficient(y_true, y_pred):
        intersection = np.sum(y_true * y_pred)
        sum_union = np.sum(y_true) + np.sum(y_pred)
        dice = 2 * intersection / sum_union
        return dice

    dice = dice_coefficient(original_mask_flat, cellpose_mask_flat)
    print(f'Dice Coefficient: {dice:.4f}')
    return


@app.cell
def _(masks, original_masks):
    cellpose_cell_count = masks.max()
    original_cell_count = original_masks.max()

    print(f"Number of cells detected by Cellpose: {cellpose_cell_count}")
    print(f"Number of cells in original segmentation: {original_cell_count}")
    return


@app.cell
def _(masks, original_masks, pd, regionprops_table):
    #cellpose cell areas
    cellpose_props = regionprops_table(masks, properties=['area'])
    cellpose_areas = pd.DataFrame(cellpose_props)
    cellpose_areas['method'] = 'Cellpose'

    #original cell areas
    original_props = regionprops_table(original_masks, properties=['area'])
    original_areas = pd.DataFrame(original_props)
    original_areas['method'] = 'Original'

    #combine data
    areas_df = pd.concat([cellpose_areas, original_areas], ignore_index=True)
    return (areas_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""compare cell sizes between the two segmentation methods""")
    return


@app.cell
def _(areas_df, plt, sns):
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=areas_df, x='area', hue='method', common_norm=False)
    plt.xlabel('Cell Area (pixels)')
    plt.ylabel('Density')
    plt.title('Cell Size Distribution Comparison')
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""and save the results to a new file""")
    return


@app.cell
def _():
    # tifffile.imwrite('/data/spatial_workshop/day1/Xenium_V1_FFPE_TgCRND8_17_9_months_outs/cellpose_masks_roi.tif', masks.astype(np.uint16))
    # tifffile.imwrite('/data/spatial_workshop/day1/Xenium_V1_FFPE_TgCRND8_17_9_months_outs/original_masks_roi.tif', original_masks.astype(np.uint16))
    return


@app.cell
def _():
    # metrics = pd.DataFrame({
    #     'Metric': ['Jaccard Index', 'Dice Coefficient'],
    #     'Value': [jaccard, dice]
    # })
    # metrics.to_csv('/data/spatial_workshop/day1/Xenium_V1_FFPE_TgCRND8_17_9_months_outs/segmentation_comparison_metrics.csv', index=False)
    return


@app.cell
def _():
    # areas_df.to_csv('/data/spatial_workshop/day1/Xenium_V1_FFPE_TgCRND8_17_9_months_outs/cell_size_comparison.csv', index=False)
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
