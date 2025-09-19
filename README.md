# Elixir Single Cell and Spatial Omics workshop materials

This repository contains modules and code to learn how to analyze spatially
resolved transcriptomics data.

## Compute and Software Requirements

This workshop uses [Pixi](https://pixi.sh/latest/) for reproducible dependency management.
Please install Pixi following the instructions on their website.
See the `pixi.toml` in the root of this repository for the list of dependencies used in this workshop.

Python exercises are run in interactive Marimo notebooks, and R exercises are run in Quarto notebooks.

The environment can also be run within a docker container available under packages on this repository.

### Supported Platforms

- Linux
- MacOS
- Windows Subsystem for Linux (WSL2)

### Hardware Requirements

- Minimum 64GB RAM
- Minimum 8 CPU cores
- Minimum 50GB free disk space

## Learning Objectives

1. **Handling spatial data:**
2. **Quality control and filtering of spatial data:**
3. **Normalization and integration of spatial data:**
4. **Dimensionality reduction and clustering of spatial data:**
5. **Visualization of spatial data:**
6. **Cell type annotation in spatial data:**
7. **Spatially variable gene identification:**
8. **Cell-cell interaction analysis in spatial data:**

## Repository Structure

- `data`: Contains scripts to fetch data used in the workshop.
- `modules`: Contains the workshop modules in Jupyter and Quarto formats.
