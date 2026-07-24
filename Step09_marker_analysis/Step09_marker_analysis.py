#!/usr/bin/env python3
"""
STEP 09 - Cluster Marker Analysis and Cell-Specific Visualization

This script completes the single-cell analysis by:
1. Running a formal One-vs-Rest marker analysis on the 13 Leiden clusters.
2. Generating UMAP Feature Plots and Dot Plots specifically for the top
   Tier 4 (Huh7-specific) viral progression markers (e.g., SELENOK, SERP1)
   to visually demonstrate their cell-specific response.

Inputs:
    - adata_processed.h5ad (from Step 02)
    - crossmodal_validated_progression_signature.csv (from Step 06)

Outputs:
    - leiden_cluster_markers.csv
    - UMAP_Tier4_Huh7_Specific_Markers.png
    - DotPlot_Tier4_Huh7_Specific_Markers.png
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

# Suppress overly verbose warnings
import warnings
warnings.filterwarnings("ignore")

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Run cluster marker analysis and Tier 4 visualizations.")
    
    # Default paths
    default_here = Path(__file__).resolve().parent
    default_root = default_here.parent
    default_adata = default_root / "Step02_clustering" / "adata_processed.h5ad"
    default_sig = default_root / "Step06_crossmodal_validation" / "crossmodal_validated_progression_signature.csv"
    
    parser.add_argument("--adata", type=Path, default=default_adata, help="Processed AnnData.")
    parser.add_argument("--sig", type=Path, default=default_sig, help="Validated progression signature.")
    parser.add_argument("--out-dir", type=Path, default=default_here, help="Output directory.")
    
    parser.add_argument("--log-level", type=str, default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return parser.parse_args()

def setup_logger(level: str) -> logging.Logger:
    """Set up the logger."""
    logger = logging.getLogger(__name__)
    logger.setLevel(getattr(logging, level.upper()))
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    if not logger.handlers:
        logger.addHandler(ch)
    return logger

def get_top_tier4_markers(sig_path: Path, top_n: int = 5, logger: logging.Logger = None) -> List[str]:
    """Load Tier 4 genes and return the top markers sorted by mean viral load correlation."""
    if not sig_path.exists():
        if logger: logger.error(f"Signature file not found: {sig_path}")
        return []
        
    df = pd.read_csv(sig_path)
    # Tier can be string or int depending on earlier steps, handle both
    df["tier_val"] = df["tier"].astype(str).str.replace("Tier", "").str.strip()
    tier4 = df[df["tier_val"] == "4"].copy()
    
    if tier4.empty:
        if logger: logger.warning("No Tier 4 genes found.")
        return []
        
    # Calculate mean Spearman rho across both viruses
    tier4["mean_rho"] = (tier4["rho_DENV_sc"] + tier4["rho_ZIKV_sc"]) / 2
    
    # Sort by correlation and get top N symbols
    top_markers = tier4.sort_values("mean_rho", ascending=False).head(top_n)["symbol"].tolist()
    if logger: logger.info(f"Top Tier 4 markers identified: {top_markers}")
    return top_markers

def main():
    args = parse_args()
    logger = setup_logger(args.log_level)
    
    logger.info("Starting Step 09: Marker Analysis")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Load AnnData
    logger.info(f"Loading AnnData from {args.adata}...")
    try:
        adata = sc.read_h5ad(args.adata)
    except Exception as e:
        logger.error(f"Failed to read AnnData: {e}")
        sys.exit(1)
        
    # Set proper .raw
    if adata.raw is not None:
        adata = adata.raw.to_adata()
        
    # Use gene symbols instead of Ensembl IDs
    adata.var_names = adata.var["symbol"].astype(str)
    adata.var_names_make_unique()
        
    # 2. Cluster Marker Analysis (One-vs-Rest)
    logger.info("Running Wilcoxon rank-sum test to identify Leiden cluster markers (One-vs-Rest)...")
    try:
        sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
        # Extract the results to a DataFrame
        result = adata.uns["rank_genes_groups"]
        groups = result["names"].dtype.names
        df = pd.DataFrame(
            {
                group + "_" + key: result[key][group]
                for group in groups
                for key in ["names", "logfoldchanges", "pvals_adj"]
            }
        )
        csv_path = args.out_dir / "leiden_cluster_markers.csv"
        df.head(50).to_csv(csv_path, index=False)
        logger.info(f"Saved top 50 markers per cluster to {csv_path}")
    except Exception as e:
        logger.error(f"Failed cluster marker analysis: {e}")
        
    # 3. Visualizations for Tier 4 (Huh7-specific) markers
    logger.info("Retrieving top Tier 4 (Huh7-specific) markers...")
    top_tier4 = get_top_tier4_markers(args.sig, top_n=6, logger=logger)
    
    if top_tier4:
        # Ensure genes are in adata
        valid_genes = [g for g in top_tier4 if g in adata.var_names]
        
        if valid_genes:
            logger.info("Generating UMAP Feature Plots...")
            sc.pl.umap(
                adata,
                color=["log_viral_load"] + valid_genes,
                cmap="viridis",
                ncols=3,
                show=False,
                title=["Log Viral Load"] + valid_genes
            )
            plt.savefig(args.out_dir / "UMAP_Tier4_Huh7_Specific_Markers.png", dpi=300, bbox_inches="tight")
            plt.close()
            
            logger.info("Generating Dot Plot for Infection States...")
            # Reorder infection_state logically if it's a category
            states = ["Mock", "Bystander", "Low", "High"]
            # Convert string series to categorical if not already
            if not pd.api.types.is_categorical_dtype(adata.obs["infection_state"]):
                 adata.obs["infection_state"] = adata.obs["infection_state"].astype("category")
                 
            adata.obs["infection_state"] = adata.obs["infection_state"].cat.reorder_categories(
                [s for s in states if s in adata.obs["infection_state"].cat.categories]
            )
                
            sc.pl.dotplot(
                adata,
                var_names=valid_genes,
                groupby="infection_state",
                standard_scale="var", # Standardize between 0 and 1 per gene for better visual comparison
                cmap="Reds",
                show=False
            )
            plt.savefig(args.out_dir / "DotPlot_Tier4_Huh7_Specific_Markers.png", dpi=300, bbox_inches="tight")
            plt.close()
            logger.info("Marker visualizations saved successfully.")
        else:
            logger.warning("None of the top Tier 4 markers were found in the AnnData var_names.")
            
    # 4. Cell State Annotation
    logger.info("Annotating Leiden clusters with biological cell states...")
    
    # Dictionary mapping clusters to cell states based on top markers
    cluster_annotations = {
        "0": "Cell Cycle (G2/M Phase)",
        "1": "Epithelial / Adhesion",
        "2": "Core Hepatocyte / Metal Transport",
        "3": "Cell Cycle (S Phase)",
        "4": "Mesenchymal-like",
        "5": "Glycolytic / High Translation",
        "6": "Mitosis (M Phase)",
        "7": "Spindle Assembly",
        "8": "Oxidative Stress Response",
        "9": "Proliferating + Early ER Stress",
        "10": "High Ribosomal Biogenesis",
        "11": "Severe ER Stress (UPR)",
        "12": "Amino Acid Starvation (ATF4)"
    }
    
    # Map the dictionary to a new obs column
    adata.obs["cell_state"] = adata.obs["leiden"].map(cluster_annotations).astype("category")
    
    logger.info("Generating UMAP colored by cell state...")
    sc.pl.umap(
        adata,
        color="cell_state",
        title="Annotated Huh7 Cell States",
        show=False,
        palette="tab20"  # Ensures enough colors for 13 states
    )
    plt.savefig(args.out_dir / "UMAP_cell_state.png", dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info("Generating Matrix Plot for top canonical markers...")
    
    # Select canonical markers that strongly define these states to plot
    canonical_markers = {
        "Cell Cycle (G2/M Phase)": ["MKI67", "CDK1", "UBE2C"],
        "Epithelial / Adhesion": ["S100A10", "LGALS3"],
        "Core Hepatocyte / Metal Transport": ["MT2A", "APOA2"],
        "Cell Cycle (S Phase)": ["PCNA", "MCM4", "MCM6"],
        "Mesenchymal-like": ["VIM", "TUBA1A"],
        "Glycolytic / High Translation": ["LDHA", "SRP9"],
        "Mitosis (M Phase)": ["CCNB1", "CDC20"],
        "Spindle Assembly": ["PLK1", "TUBB4B"],
        "Oxidative Stress Response": ["TXN", "PRDX1"],
        "Proliferating + Early ER Stress": ["SDF2L1", "SELENOK"],
        "High Ribosomal Biogenesis": ["RPS15", "RPL31"],
        "Severe ER Stress (UPR)": ["DDIT3", "HERPUD1"],
        "Amino Acid Starvation (ATF4)": ["PSAT1", "PLIN2"]
    }
    
    # Flatten list of markers, ensuring they are in adata
    markers_to_plot = {}
    for state, genes in canonical_markers.items():
        valid_genes = [g for g in genes if g in adata.var_names]
        if valid_genes:
            markers_to_plot[state] = valid_genes
            
    sc.pl.matrixplot(
        adata,
        var_names=markers_to_plot,
        groupby="cell_state",
        dendrogram=True, # Groups similar states together on the plot
        standard_scale="var", # Standardize for clear visualization
        cmap="Blues",
        show=False
    )
    plt.savefig(args.out_dir / "MatrixPlot_Cell_States.png", dpi=300, bbox_inches="tight")
    plt.close()
            
    logger.info("Step 09 completed successfully.")

if __name__ == "__main__":
    main()
