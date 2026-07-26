#!/usr/bin/env python3
"""
STEP 02 - Normalize, HVG, PCA, UMAP, Leiden.

This script normalizes the filtered AnnData, identifies highly variable genes (HVGs),
and performs dimensionality reduction (PCA, UMAP) and clustering (Leiden).
It generates UMAP plots colored by condition, infection state, clustering, and viral load.

Inputs:
    - adata_qc.h5ad (from Step01_qc)

Outputs:
    - adata_processed.h5ad: Processed AnnData object
    - umap_*.png: UMAP projection plots
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
import scanpy as sc

warnings.filterwarnings("ignore")

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Perform normalization, dimensionality reduction, and clustering.")
    
    # Default paths
    import sys
    _root_dir = Path(__file__).resolve().parent.parent.parent
    if str(_root_dir) not in sys.path: sys.path.append(str(_root_dir))
    from src import config
    default_here = config.PROCESSED_DATA_DIR / Path(__file__).resolve().parent.name
    default_root = config.PROCESSED_DATA_DIR
    default_in = default_root / "Step01_qc" / "adata_qc.h5ad"
    default_out_h5ad = default_here / "adata_processed.h5ad"

    # IO arguments
    parser.add_argument("--in-file", type=Path, default=default_in,
                        help="Path to the input QC'd AnnData object (.h5ad).")
    parser.add_argument("--out-h5ad", type=Path, default=default_out_h5ad,
                        help="Path to save the output processed AnnData object.")
    
    # Algorithm parameters
    parser.add_argument("--target-sum", type=float, default=1e4,
                        help="Target sum for total count normalization.")
    parser.add_argument("--n-top-genes", type=int, default=3000,
                        help="Number of highly variable genes to select.")
    parser.add_argument("--n-pcs", type=int, default=50,
                        help="Number of principal components to compute and use.")
    parser.add_argument("--leiden-resolution", type=float, default=1.0,
                        help="Resolution parameter for Leiden clustering.")
    parser.add_argument("--random-state", type=int, default=0,
                        help="Random seed for reproducible results.")

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

def process_data(adata: sc.AnnData, args: argparse.Namespace, logger: logging.Logger) -> sc.AnnData:
    """Normalize, compute HVGs, PCA, UMAP, and Leiden clusters."""
    logger.info("Setting .raw to the current state...")
    adata.raw = adata
    
    logger.info(f"Normalizing total counts (target_sum={args.target_sum}) and applying log1p...")
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    
    logger.info(f"Identifying {args.n_top_genes} highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, flavor="seurat")
    
    logger.info("Subsetting to HVGs for dimensionality reduction...")
    ah = adata[:, adata.var.highly_variable].copy()
    
    logger.info("Regressing out percent mitochondrial reads and scaling data...")
    sc.pp.regress_out(ah, ["pct_mito"])
    sc.pp.scale(ah, max_value=10)
    
    logger.info(f"Computing PCA (n_comps={args.n_pcs})...")
    sc.tl.pca(ah, n_comps=args.n_pcs, svd_solver="arpack", random_state=args.random_state)
    
    logger.info("Computing neighborhood graph and UMAP...")
    sc.pp.neighbors(ah, n_pcs=args.n_pcs, random_state=args.random_state)
    sc.tl.umap(ah, random_state=args.random_state)
    
    logger.info(f"Performing Leiden clustering (resolution={args.leiden_resolution})...")
    sc.tl.leiden(ah, resolution=args.leiden_resolution, random_state=args.random_state, flavor="igraph", n_iterations=2, directed=False)
    
    logger.info("Transferring UMAP coordinates and Leiden clusters back to full object...")
    adata.obsm["X_umap"] = ah.obsm["X_umap"]
    adata.obs["leiden"] = ah.obs["leiden"].values
    
    logger.info("Computing log_viral_load...")
    adata.obs["log_viral_load"] = np.log1p(adata.obs["viral_load"].astype(float))
    
    return adata

def generate_umap_plots(adata: sc.AnnData, out_dir: Path, logger: logging.Logger) -> None:
    """Generate and save UMAP plots."""
    logger.info("Generating UMAP plots...")
    
    plot_configs = [
        ("condition", "umap_condition.png", None),
        ("infection_state", "umap_infection_state.png", None),
        ("leiden", "umap_leiden.png", None),
        ("log_viral_load", "umap_viral_load.png", "viridis")
    ]
    
    for color, fn, cmap in plot_configs:
        fig, ax = plt.subplots(figsize=(6, 5))
        sc.pl.umap(adata, color=color, ax=ax, show=False, cmap=cmap, frameon=False)
        fig.tight_layout()
        
        out_path = out_dir / fn
        try:
            fig.savefig(out_path, dpi=150)
            logger.debug(f"Saved {fn}")
        except Exception as e:
            logger.error(f"Failed to save plot {fn}: {e}")
        finally:
            plt.close(fig)
            
    logger.info(f"All UMAP plots saved to {out_dir}")

def main():
    """Main execution entry point."""
    args = parse_args()
    logger = setup_logger(args.log_level)
    
    logger.info("Starting Step 02: Clustering and Dimensionality Reduction")
    
    if not args.in_file.exists():
        logger.error(f"Input file not found: {args.in_file}")
        sys.exit(1)
        
    logger.info(f"Loading QC'd AnnData from {args.in_file}...")
    try:
        adata = sc.read_h5ad(args.in_file)
    except Exception as e:
        logger.error(f"Failed to read input AnnData: {e}")
        sys.exit(1)
        
    adata = process_data(adata, args, logger)
    
    n_clusters = adata.obs['leiden'].nunique()
    logger.info(f"Identified {n_clusters} Leiden clusters.")
    
    # Generate Plots
    out_dir = args.out_h5ad.parent
    generate_umap_plots(adata, out_dir, logger)
    
    # Save output
    logger.info(f"Saving processed AnnData to {args.out_h5ad}...")
    try:
        args.out_h5ad.parent.mkdir(parents=True, exist_ok=True)
        adata.write(args.out_h5ad)
        logger.info("Step 02 completed successfully.")
    except Exception as e:
        logger.error(f"Failed to save processed AnnData: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
