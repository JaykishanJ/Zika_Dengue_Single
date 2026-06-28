# Data Repository

This directory contains the essential lightweight metadata, summary statistics, and interaction prediction tables required to rapidly reproduce the downstream steps of the pipeline (such as Cross-Modal Signatures, miRNA target networks, and PPI mapping).

### Important Note on Single-Cell Data (`.h5ad`)
Due to GitHub's strict file size limits (100 MB per file), the primary annotated single-cell RNA-seq objects (`adata_annotated.h5ad` and `adata_trajectory.h5ad`), which collectively exceed 1.2 GB, cannot be hosted directly in this repository.

To fully reproduce Steps 03 through 05 (Single-Cell State, Pseudotime, and CellChat analyses), you must download the single-cell objects from our external hosting repository (e.g., Zenodo / OSF - Link pending publication) and place them in the root of your local analysis directory.
