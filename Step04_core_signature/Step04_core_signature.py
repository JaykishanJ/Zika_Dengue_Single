#!/usr/bin/env python3
"""
STEP 04 - Convergent core signature.

This script identifies the convergent core signature across viruses by finding
genes that are significantly differentially expressed and change in the same
direction (concordant) in both DENV and ZIKV infections.

Inputs:
    - DE_High_vs_Mock_DENV.csv (from Step03_differential_expression)
    - DE_High_vs_Mock_ZIKV.csv (from Step03_differential_expression)

Outputs:
    - convergent_core_signature.csv: Dataframe of the core signature genes
    - core_heatmap.png: Heatmap of log2FC for concordant genes
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

warnings.filterwarnings("ignore")

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Identify convergent core signature across viruses.")
    
    # Default paths
    default_here = Path(__file__).resolve().parent
    default_root = default_here.parent
    default_de_dir = default_root / "Step03_differential_expression"
    default_out_csv = default_here / "convergent_core_signature.csv"
    default_out_plot = default_here / "core_heatmap.png"

    # IO arguments
    parser.add_argument("--de-dir", type=Path, default=default_de_dir,
                        help="Directory containing the differential expression CSVs.")
    parser.add_argument("--out-csv", type=Path, default=default_out_csv,
                        help="Path to save the core signature CSV.")
    parser.add_argument("--out-plot", type=Path, default=default_out_plot,
                        help="Path to save the core signature heatmap.")

    # Threshold arguments
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="False Discovery Rate (padj) threshold.")
    parser.add_argument("--lfc", type=float, default=0.58,
                        help="Log2 Fold Change threshold.")
    
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

def get_significant_genes(virus: str, de_dir: Path, fdr: float, lfc: float, logger: logging.Logger) -> pd.DataFrame:
    """Load DE results and filter for significant genes."""
    file_path = de_dir / f"DE_High_vs_Mock_{virus}.csv"
    if not file_path.exists():
        logger.error(f"DE results not found for {virus}: {file_path}")
        raise FileNotFoundError(f"Missing file: {file_path}")
        
    r = pd.read_csv(file_path)
    
    mask = (r.padj < fdr) & (r.log2FoldChange.abs() > lfc)
    s = r[mask].dropna(subset=["symbol"]).copy()
    
    s["dir"] = np.where(s.log2FoldChange > 0, "up", "down")
    s = s.drop_duplicates("symbol").set_index("symbol")
    
    logger.info(f"Identified {len(s)} significant DEGs for {virus}.")
    return s

def build_core_signature(d: pd.DataFrame, z: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """Find intersecting genes and determine concordance."""
    genes = sorted(set(d.index) & set(z.index))
    
    rows: List[Dict[str, Any]] = []
    for g in genes:
        rows.append({
            "symbol": g, 
            "DENV_log2FC": d.loc[g, "log2FoldChange"], 
            "DENV_padj": d.loc[g, "padj"], 
            "DENV_dir": d.loc[g, "dir"],
            "ZIKV_log2FC": z.loc[g, "log2FoldChange"], 
            "ZIKV_padj": z.loc[g, "padj"], 
            "ZIKV_dir": z.loc[g, "dir"]
        })
        
    core = pd.DataFrame(rows)
    
    if not core.empty:
        core["concordant"] = core["DENV_dir"] == core["ZIKV_dir"]
        core["direction"] = np.where(core["concordant"], core["DENV_dir"], "discordant")
    else:
        # Create empty dataframe with expected columns
        core = pd.DataFrame(columns=["symbol", "DENV_log2FC", "DENV_padj", "DENV_dir", "ZIKV_log2FC", "ZIKV_padj", "ZIKV_dir", "concordant", "direction"])
        
    return core

def plot_heatmap(conc: pd.DataFrame, out_plot: Path, logger: logging.Logger) -> None:
    """Generate and save a heatmap for the concordant genes."""
    if conc.empty:
        logger.warning("No concordant genes found; skipping heatmap generation.")
        return
        
    hm = conc.set_index("symbol")[["DENV_log2FC", "ZIKV_log2FC"]].sort_values("DENV_log2FC")
    
    height = max(3, 0.3 * len(hm))
    fig, ax = plt.subplots(figsize=(3.4, height))
    
    im = ax.imshow(hm.values, cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
    
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["DENV", "ZIKV"])
    ax.set_yticks(range(len(hm)))
    ax.set_yticklabels(hm.index, fontsize=8)
    
    fig.colorbar(im, ax=ax, label="log2FC", shrink=0.5)
    ax.set_title("Convergent core")
    
    fig.tight_layout()
    try:
        fig.savefig(out_plot, dpi=150)
        logger.info(f"Saved core heatmap to {out_plot}")
    except Exception as e:
        logger.error(f"Failed to save core heatmap: {e}")
    finally:
        plt.close(fig)

def main():
    """Main execution entry point."""
    args = parse_args()
    logger = setup_logger(args.log_level)
    
    logger.info("Starting Step 04: Convergent Core Signature")
    
    try:
        d = get_significant_genes("DENV", args.de_dir, args.fdr, args.lfc, logger)
        z = get_significant_genes("ZIKV", args.de_dir, args.fdr, args.lfc, logger)
    except FileNotFoundError:
        sys.exit(1)
        
    core = build_core_signature(d, z, logger)
    
    # Save CSV
    try:
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        core.to_csv(args.out_csv, index=False)
        logger.info(f"Saved core signature to {args.out_csv}")
    except Exception as e:
        logger.error(f"Failed to save core signature CSV: {e}")
        sys.exit(1)
        
    # Logging summaries
    conc = core[core["concordant"]] if not core.empty else pd.DataFrame()
    n_up = (conc["direction"] == "up").sum() if not conc.empty else 0
    n_down = (conc["direction"] == "down").sum() if not conc.empty else 0
    
    logger.info(f"Total significant in both: {len(core)}")
    logger.info(f"Concordant genes: {len(conc)} ({n_up} up, {n_down} down)")
    
    # Plotting
    plot_heatmap(conc, args.out_plot, logger)
    
    logger.info("Step 04 completed successfully.")

if __name__ == "__main__":
    main()
