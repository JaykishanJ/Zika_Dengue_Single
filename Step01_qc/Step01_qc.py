#!/usr/bin/env python3
"""
STEP 01 - Quality control.

This script performs quality control filtering on the raw single-cell AnnData.
It drops spike-ins, calculates QC metrics, generates pre-filter violin plots,
and filters cells based on gene counts and mitochondrial content thresholds.

Inputs:
    - adata_raw.h5ad (from Step00_build_data)

Outputs:
    - adata_qc.h5ad: Filtered AnnData object, with a raw counts layer
    - qc_violin_prefilter.png: Violin plots of QC metrics before filtering
    - qc_summary.csv: Summary of cell counts before and after filtering
    - cells_condition_x_state.csv: Crosstab of condition vs infection_state
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Perform single-cell quality control.")
    
    # Default paths
    default_here = Path(__file__).resolve().parent
    default_root = default_here.parent
    default_in = default_root / "Step00_build_data" / "adata_raw.h5ad"
    default_out_h5ad = default_here / "adata_qc.h5ad"
    default_out_plot = default_here / "qc_violin_prefilter.png"
    default_out_summary = default_here / "qc_summary.csv"
    default_out_crosstab = default_here / "cells_condition_x_state.csv"

    # IO arguments
    parser.add_argument("--in-file", type=Path, default=default_in,
                        help="Path to the input raw AnnData object (.h5ad).")
    parser.add_argument("--out-h5ad", type=Path, default=default_out_h5ad,
                        help="Path to save the output QC'd AnnData object.")
    parser.add_argument("--out-plot", type=Path, default=default_out_plot,
                        help="Path to save the QC violin plot.")
    parser.add_argument("--out-summary", type=Path, default=default_out_summary,
                        help="Path to save the QC summary CSV.")
    parser.add_argument("--out-crosstab", type=Path, default=default_out_crosstab,
                        help="Path to save the condition vs state crosstab CSV.")

    # Threshold arguments
    parser.add_argument("--min-genes", type=int, default=2000,
                        help="Minimum number of genes expressed required for a cell to pass filtering.")
    parser.add_argument("--max-genes", type=int, default=8000,
                        help="Maximum number of genes expressed required for a cell to pass filtering.")
    parser.add_argument("--max-mito", type=float, default=15.0,
                        help="Maximum percentage of mitochondrial counts for a cell to pass filtering.")
    
    # Logging
    parser.add_argument("--log-level", type=str, default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Set the logging level.")
    return parser.parse_args()

def setup_logger(level: str) -> logging.Logger:
    """Set up the logger for the script."""
    logger = logging.getLogger(__name__)
    logger.setLevel(getattr(logging, level.upper()))
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    if not logger.handlers:
        logger.addHandler(ch)
    return logger

def plot_qc_metrics(adata: sc.AnnData, out_plot: Path, min_genes: int, max_genes: int, max_mito: float, logger: logging.Logger) -> None:
    """Generate and save QC violin plots."""
    logger.info("Generating QC violin plots...")
    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    
    metrics = ["n_genes", "total_counts", "pct_mito"]
    labels = ["Genes/cell", "Total counts", "% mito"]
    
    for a, c, lab in zip(ax, metrics, labels):
        sc.pl.violin(adata, c, jitter=0.4, ax=a, show=False)
        a.set_title(lab)
        
    ax[0].axhline(min_genes, color="r", ls="--")
    ax[0].axhline(max_genes, color="r", ls="--")
    ax[2].axhline(max_mito, color="r", ls="--")
    
    fig.tight_layout()
    try:
        fig.savefig(out_plot, dpi=150)
        logger.info(f"Saved QC plots to {out_plot}")
    except Exception as e:
        logger.error(f"Failed to save QC plot: {e}")
    finally:
        plt.close(fig)

def filter_cells(adata: sc.AnnData, min_genes: int, max_genes: int, max_mito: float, logger: logging.Logger) -> sc.AnnData:
    """Filter cells based on established thresholds."""
    logger.info(f"Filtering cells: {min_genes} <= n_genes <= {max_genes}, pct_mito < {max_mito}")
    mask = (adata.obs.n_genes >= min_genes) & (adata.obs.n_genes <= max_genes) & (adata.obs.pct_mito < max_mito)
    
    adata_filtered = adata[mask.to_numpy()].copy()
    adata_filtered.layers["counts"] = adata_filtered.X.copy()
    
    return adata_filtered

def main():
    """Main execution entry point."""
    args = parse_args()
    logger = setup_logger(args.log_level)
    
    logger.info("Starting Step 01: Quality Control")
    
    if not args.in_file.exists():
        logger.error(f"Input file not found: {args.in_file}")
        sys.exit(1)
        
    logger.info(f"Loading input AnnData from {args.in_file}...")
    try:
        adata = sc.read_h5ad(args.in_file)
    except Exception as e:
        logger.error(f"Failed to read input AnnData: {e}")
        sys.exit(1)
        
    logger.info("Dropping spike-ins and calculating QC metrics...")
    adata = adata[:, ~adata.var["is_spikein"].to_numpy()].copy()
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], percent_top=None, log1p=False, inplace=True)
    
    # Store metrics under cleaner names
    adata.obs["pct_mito"] = adata.obs["pct_counts_mito"]
    adata.obs["n_genes"] = adata.obs["n_genes_by_counts"]
    
    n_before = adata.n_obs
    logger.info(f"Number of cells before filtering: {n_before}")
    
    # Generate Pre-filter plots
    plot_qc_metrics(adata, args.out_plot, args.min_genes, args.max_genes, args.max_mito, logger)
    
    # Filter cells
    adata = filter_cells(adata, args.min_genes, args.max_genes, args.max_mito, logger)
    n_after = adata.n_obs
    logger.info(f"Number of cells after filtering: {n_after} ({(n_after/n_before)*100:.1f}%)")
    
    # Save outputs
    logger.info("Saving summaries and filtered AnnData...")
    try:
        # Save summary
        summary_df = pd.DataFrame({"metric": ["cells_before", "cells_after"], "value": [n_before, n_after]})
        summary_df.to_csv(args.out_summary, index=False)
        
        # Save crosstab
        crosstab_df = pd.crosstab(adata.obs["condition"], adata.obs["infection_state"])
        crosstab_df.to_csv(args.out_crosstab)
        
        # Save h5ad
        args.out_h5ad.parent.mkdir(parents=True, exist_ok=True)
        adata.write(args.out_h5ad)
        
        logger.info(f"QC completed successfully. Outputs saved to {args.out_h5ad.parent}")
    except Exception as e:
        logger.error(f"Failed to save one or more outputs: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
