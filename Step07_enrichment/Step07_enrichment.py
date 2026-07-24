#!/usr/bin/env python3
"""
STEP 07 - Pathway enrichment of the cross-modal convergent progression signature.

This script performs pathway enrichment analysis on Tier 1 and Tier 2 convergent
progression genes using Enrichr (via gseapy). It queries multiple gene set libraries
and outputs the significant enriched terms along with a summary barplot.

Inputs:
    - crossmodal_validated_progression_signature.csv (from Step06)

Outputs:
    - enrichment_up.csv: Enriched pathways for upregulated genes
    - enrichment_down.csv: Enriched pathways for downregulated genes
    - enrichment_up_barplot.png: Visualization of top upregulated pathways

Note: Requires an active internet connection to query the Enrichr API.
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import List

import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")

# Default gene set libraries for enrichment
DEFAULT_LIBS = [
    "MSigDB_Hallmark_2020",
    "GO_Biological_Process_2021",
    "KEGG_2021_Human",
    "Reactome_2022"
]

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Perform pathway enrichment analysis via Enrichr.")
    
    # Default paths
    default_here = Path(__file__).resolve().parent
    default_root = default_here.parent
    default_sig = default_root / "Step06_crossmodal_validation" / "crossmodal_validated_progression_signature.csv"

    # Input Files
    parser.add_argument("--sig-file", type=Path, default=default_sig,
                        help="Path to validated progression signature CSV.")
    
    # Output Dir
    parser.add_argument("--out-dir", type=Path, default=default_here,
                        help="Directory to save enrichment outputs.")

    # Thresholds
    parser.add_argument("--max-tier", type=int, default=2,
                        help="Maximum tier to include (default: 2, includes Tier 1 and 2).")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="Adjusted p-value threshold for enrichment results.")

    # Logging
    parser.add_argument("--log-level", type=str, default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Set the logging level.")
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

def run_enrichr(genes: List[str], tag: str, fdr: float, out_dir: Path, logger: logging.Logger) -> pd.DataFrame:
    """Query Enrichr API and return filtered significant results."""
    if not genes:
        logger.warning(f"No genes provided for '{tag}' enrichment.")
        return pd.DataFrame()
        
    logger.info(f"Querying Enrichr for {len(genes)} '{tag}' genes...")
    try:
        enr = gp.enrichr(
            gene_list=genes,
            gene_sets=DEFAULT_LIBS,
            organism="human",
            outdir=None
        )
        r = enr.results
        
        # Filter and sort
        r = r[r["Adjusted P-value"] < fdr].sort_values("Adjusted P-value")
        
        out_csv = out_dir / f"enrichment_{tag}.csv"
        r.to_csv(out_csv, index=False)
        logger.info(f"Found {len(r)} significant terms for '{tag}'. Saved to {out_csv.name}")
        
        return r
    except Exception as e:
        logger.error(f"Enrichment for '{tag}' failed: {e}")
        return pd.DataFrame()

def plot_enrichment(ru: pd.DataFrame, out_path: Path, logger: logging.Logger) -> None:
    """Plot the top 15 enriched terms."""
    logger.info("Generating enrichment barplot...")
    ru = ru.copy()
    
    ru["nlp"] = -np.log10(ru["Adjusted P-value"].clip(lower=1e-300))
    ru["short"] = ru["Term"].str.replace(r"\(GO:\d+\)", "", regex=True) \
                            .str.replace(r"R-HSA-\d+", "", regex=True) \
                            .str.strip()
                            
    pal = {
        "MSigDB_Hallmark_2020": "#c0392b",
        "GO_Biological_Process_2021": "#2471a3",
        "KEGG_2021_Human": "#27ae60",
        "Reactome_2022": "#8e44ad"
    }
    
    top = ru.drop_duplicates("short").head(15).sort_values("nlp")
    
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = [pal.get(g, "#888888") for g in top["Gene_set"]]
    
    ax.barh(range(len(top)), top["nlp"], color=colors)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["short"], fontsize=8)
    ax.set_xlabel("-log10 adj P")
    ax.set_title("Convergent progression (up) — enriched pathways")
    
    handles = [Patch(color=c, label=l.split('_')[0]) for l, c in pal.items()]
    ax.legend(handles=handles, fontsize=7, loc="lower right")
    
    fig.tight_layout()
    try:
        fig.savefig(out_path, dpi=150)
        logger.info(f"Saved barplot to {out_path}")
    except Exception as e:
        logger.error(f"Failed to save barplot: {e}")
    finally:
        plt.close(fig)

def main():
    """Main execution entry point."""
    args = parse_args()
    logger = setup_logger(args.log_level)
    
    logger.info("Starting Step 07: Pathway Enrichment")
    
    if not args.sig_file.exists():
        logger.error(f"Signature file not found: {args.sig_file}")
        sys.exit(1)
        
    try:
        sig = pd.read_csv(args.sig_file)
    except Exception as e:
        logger.error(f"Failed to read signature file: {e}")
        sys.exit(1)
        
    sig = sig[sig["tier"] <= args.max_tier]
    
    up = sorted(sig[sig["sc_direction"] == "up"]["symbol"].dropna().unique())
    down = sorted(sig[sig["sc_direction"] == "down"]["symbol"].dropna().unique())
    
    logger.info(f"Tier 1-{args.max_tier} signature -> up={len(up)}, down={len(down)}")
    
    args.out_dir.mkdir(parents=True, exist_ok=True)
    
    ru = run_enrichr(up, "up", args.fdr, args.out_dir, logger)
    run_enrichr(down, "down", args.fdr, args.out_dir, logger)
    
    if not ru.empty:
        plot_path = args.out_dir / "enrichment_up_barplot.png"
        plot_enrichment(ru, plot_path, logger)
    else:
        logger.warning("No significant up-regulated terms found, skipping plot.")
        
    logger.info("Step 07 completed successfully.")

if __name__ == "__main__":
    main()
