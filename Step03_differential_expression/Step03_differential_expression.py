#!/usr/bin/env python3
"""
STEP 03 - Differential expression (pseudobulk DESeq2, batch-matched High vs Mock).

This script performs differential expression analysis using a pseudobulk approach
with DESeq2. It compares 'High' infection state cells against matched 'Mock' controls
for both DENV and ZIKV, producing DE statistics and volcano plots.

Inputs:
    - adata_processed.h5ad (from Step02_clustering)

Outputs:
    - DE_High_vs_Mock_DENV.csv: Differential expression results for DENV
    - DE_High_vs_Mock_ZIKV.csv: Differential expression results for ZIKV
    - volcano_DENV.png: Volcano plot for DENV
    - volcano_ZIKV.png: Volcano plot for ZIKV
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Run pseudobulk differential expression analysis.")
    
    # Default paths
    default_here = Path(__file__).resolve().parent
    default_root = default_here.parent
    default_in = default_root / "Step02_clustering" / "adata_processed.h5ad"

    # IO arguments
    parser.add_argument("--in-file", type=Path, default=default_in,
                        help="Path to the processed AnnData object (.h5ad).")
    parser.add_argument("--out-dir", type=Path, default=default_here,
                        help="Directory to save the differential expression outputs.")

    # Algorithm parameters
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="False Discovery Rate (padj) threshold.")
    parser.add_argument("--lfc", type=float, default=0.58,
                        help="Log2 Fold Change threshold.")
    parser.add_argument("--min-cells", type=int, default=5,
                        help="Minimum number of cells required per pseudobulk sample.")
    
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

def run_deseq2_pseudobulk(adata: sc.AnnData, virus: str, min_cells: int, logger: logging.Logger) -> pd.DataFrame:
    """Prepare pseudobulk counts and run DESeq2 for a given virus."""
    logger.info(f"Preparing pseudobulk data for {virus}...")
    
    # Pre-defined mapping of virus to its corresponding Mock experiment batch
    virus_experiment = {"DENV": "10017006", "ZIKV": "10017008"}
    
    if virus not in virus_experiment:
        raise ValueError(f"Unknown virus: {virus}. Must be 'DENV' or 'ZIKV'.")
        
    exp = virus_experiment[virus]
    
    is_high = (adata.obs["condition"] == virus) & (adata.obs["infection_state"] == "High")
    is_mock = (adata.obs["infection_state"] == "Mock") & (adata.obs["experiment"] == exp)
    
    sub = adata[(is_high | is_mock).to_numpy()].copy()
    sub.obs["grp"] = np.where(sub.obs["infection_state"] == "Mock", "Mock", "High")
    
    rows: List[np.ndarray] = []
    meta: List[Dict[str, Any]] = []
    
    for (t, g), n in sub.obs.groupby(["time_h", "grp"], observed=True).size().items():
        if n < min_cells:
            continue
        m = ((sub.obs.time_h == t) & (sub.obs.grp == g)).to_numpy()
        rows.append(np.asarray(sub.layers["counts"][m].sum(0)).ravel())
        meta.append({
            "sample": f"{g}_{t}h", 
            "time_h": str(t), 
            "grp": g, 
            "n_cells": int(n)
        })
        
    if not rows:
        raise ValueError(f"No valid pseudobulk samples constructed for {virus}.")
        
    info = pd.DataFrame(meta).set_index("sample")
    cnt = pd.DataFrame(np.vstack(rows), columns=adata.var_names, index=info.index)
    
    # Require paired timepoints
    paired = [t for t in info.time_h.unique() if info[info.time_h == t].grp.nunique() == 2]
    info = info[info.time_h.isin(paired)]
    cnt = cnt.loc[info.index]
    
    if info.empty:
        raise ValueError(f"No paired timepoints left for {virus} after filtering.")
        
    # Filter genes
    logger.info(f"Running DESeq2 on {len(info)} samples for {virus}...")
    keep = (cnt.sum(0) >= 10) & ((cnt > 0).sum(0) >= 2)
    cnt = cnt.loc[:, keep].round().astype(int)
    info["grp"] = pd.Categorical(info.grp, categories=["Mock", "High"])
    
    dds = DeseqDataSet(counts=cnt, metadata=info, design="~time_h + grp", quiet=True)
    dds.deseq2()
    
    st = DeseqStats(dds, contrast=["grp", "High", "Mock"], quiet=True)
    st.summary()
    
    res = st.results_df.copy()
    symbols = adata.var["symbol"].astype(str)
    res["symbol"] = symbols.reindex(res.index).values
    
    return res.sort_values("padj")

def plot_volcano(res: pd.DataFrame, virus: str, fdr: float, lfc: float, out_path: Path, logger: logging.Logger) -> None:
    """Generate and save a volcano plot for the DE results."""
    logger.info(f"Generating volcano plot for {virus}...")
    
    r = res.dropna(subset=["padj"]).copy()
    r["nlp"] = -np.log10(r.padj.clip(lower=1e-300))
    
    up = (r.padj < fdr) & (r.log2FoldChange > lfc)
    dn = (r.padj < fdr) & (r.log2FoldChange < -lfc)
    
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(r.log2FoldChange, r.nlp, s=5, c="#d5d8dc")
    ax.scatter(r.log2FoldChange[up], r.nlp[up], s=8, c="#c0392b", label=f"up ({int(up.sum())})")
    ax.scatter(r.log2FoldChange[dn], r.nlp[dn], s=8, c="#2471a3", label=f"down ({int(dn.sum())})")
    
    ax.axvline(lfc, ls="--", c="k", lw=0.5)
    ax.axvline(-lfc, ls="--", c="k", lw=0.5)
    
    ax.set(xlabel="log2FC (High vs Mock)", ylabel="-log10 padj", title=virus)
    ax.legend(fontsize=8)
    
    fig.tight_layout()
    try:
        fig.savefig(out_path, dpi=150)
    except Exception as e:
        logger.error(f"Failed to save volcano plot for {virus}: {e}")
    finally:
        plt.close(fig)

def main():
    """Main execution entry point."""
    args = parse_args()
    logger = setup_logger(args.log_level)
    
    logger.info("Starting Step 03: Differential Expression")
    
    if not args.in_file.exists():
        logger.error(f"Input file not found: {args.in_file}")
        sys.exit(1)
        
    logger.info(f"Loading processed AnnData from {args.in_file}...")
    try:
        adata = sc.read_h5ad(args.in_file)
    except Exception as e:
        logger.error(f"Failed to read input AnnData: {e}")
        sys.exit(1)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    for virus in ["DENV", "ZIKV"]:
        try:
            res = run_deseq2_pseudobulk(adata, virus, args.min_cells, logger)
            
            csv_path = args.out_dir / f"DE_High_vs_Mock_{virus}.csv"
            res.to_csv(csv_path)
            logger.info(f"Saved DE results to {csv_path}")
            
            # Print summary
            sig = res[(res.padj < args.fdr) & (res.log2FoldChange.abs() > args.lfc)]
            n_up = int((sig.log2FoldChange > 0).sum())
            n_down = int((sig.log2FoldChange < 0).sum())
            logger.info(f"{virus} Summary: {len(sig)} DEGs ({n_up} up, {n_down} down)")
            
            # Volcano Plot
            plot_path = args.out_dir / f"volcano_{virus}.png"
            plot_volcano(res, virus, args.fdr, args.lfc, plot_path, logger)
            
        except Exception as e:
            logger.error(f"Failed processing virus {virus}: {e}")
            continue

    logger.info("Step 03 completed successfully.")

if __name__ == "__main__":
    main()
