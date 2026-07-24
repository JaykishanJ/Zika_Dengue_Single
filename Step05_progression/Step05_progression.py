#!/usr/bin/env python3
"""
STEP 05 - Convergent progression signature.

This script identifies genes whose expression levels continuously track viral load
across single cells, in the same direction, for both DENV and ZIKV (convergent progression).
It calculates Spearman rank correlation between gene expression and log viral load.

Inputs:
    - adata_processed.h5ad (from Step02_clustering)
    - convergent_core_signature.csv (from Step04_core_signature)

Outputs:
    - convergent_progression_signature.csv: Core progression signature genes
    - progression_trends.png: Gene expression vs. viral load trends for top genes
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import Dict, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import rankdata, t as tdist

warnings.filterwarnings("ignore")

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Identify convergent progression genes.")
    
    # Default paths
    default_here = Path(__file__).resolve().parent
    default_root = default_here.parent
    default_in = default_root / "Step02_clustering" / "adata_processed.h5ad"
    default_core = default_root / "Step04_core_signature" / "convergent_core_signature.csv"
    default_out_csv = default_here / "convergent_progression_signature.csv"
    default_out_plot = default_here / "progression_trends.png"

    # Input Files
    parser.add_argument("--in-file", type=Path, default=default_in, help="Path to processed AnnData (.h5ad).")
    parser.add_argument("--core-file", type=Path, default=default_core, help="Path to convergent core signature CSV.")

    # Outputs
    parser.add_argument("--out-csv", type=Path, default=default_out_csv, help="Path to save progression signature CSV.")
    parser.add_argument("--out-plot", type=Path, default=default_out_plot, help="Path to save progression trends plot.")

    # Thresholds
    parser.add_argument("--rho-min", type=float, default=0.20, help="Minimum absolute Spearman correlation.")
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR threshold for significance.")
    
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

def spearman_vec(E: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Vectorized calculation of Spearman rank correlation and p-values."""
    ry = rankdata(y).astype(float)
    ry -= ry.mean()
    
    RE = rankdata(E, axis=1).astype(float)
    RE -= RE.mean(1, keepdims=True)
    
    den = np.sqrt((RE**2).sum(1) * (ry**2).sum())
    rho = np.where(den > 0, RE @ ry / den, 0.0)
    
    n = len(y)
    with np.errstate(divide="ignore", invalid="ignore"):
        tt = rho * np.sqrt((n - 2) / np.clip(1 - rho**2, 1e-12, None))
    
    p = 2 * tdist.sf(np.abs(tt), n - 2)
    return rho, p

def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction."""
    p = np.asarray(p)
    n = len(p)
    o = np.argsort(p)
    adj = np.empty(n)
    
    adj[o] = p[o] * n / (np.arange(n) + 1)
    run = 1.0
    for i in o[::-1]:
        run = min(run, adj[i])
        adj[i] = run
    return np.clip(adj, 0, 1)

def run_correlation(adata: sc.AnnData, virus: str, fdr: float, rho_min: float, logger: logging.Logger) -> pd.DataFrame:
    """Calculate correlation between gene expression and viral load for a virus."""
    sub = adata[(adata.obs.condition == virus).to_numpy()]
    Xc = sub.layers["counts"]
    
    # Filter genes expressed in at least 10% of cells for this condition
    det = np.asarray((Xc > 0).sum(0)).ravel() >= 0.10 * sub.n_obs
    idx = np.where(det)[0]
    
    E = sub.X[:, idx]
    E = np.asarray(E.todense()).T if hasattr(E, "todense") else np.asarray(E).T
    
    rho, p = spearman_vec(E, sub.obs.log_viral_load.to_numpy())
    padj = benjamini_hochberg(p)
    
    symbols = adata.var["symbol"].astype(str)
    res = pd.DataFrame({"symbol": symbols.values[idx], "rho": rho, "padj": padj}, index=adata.var_names[idx])
    
    sig = res[(res.padj < fdr) & (res.rho.abs() > rho_min)]
    logger.info(f"{virus}: {sub.n_obs} cells -> {len(sig)} progression-associated genes.")
    return res

def build_progression_signature(d: pd.DataFrame, z: pd.DataFrame, core_file: Path, fdr: float, rho_min: float, logger: logging.Logger) -> pd.DataFrame:
    """Find concordant genes across both viruses that meet thresholds."""
    common = d.index.intersection(z.index)
    
    m = pd.DataFrame({
        "symbol": d.loc[common, "symbol"],
        "rho_DENV": d.loc[common, "rho"],
        "padj_DENV": d.loc[common, "padj"],
        "rho_ZIKV": z.loc[common, "rho"],
        "padj_ZIKV": z.loc[common, "padj"]
    })
    
    # Significant in both
    sig_both = (m.padj_DENV < fdr) & (m.rho_DENV.abs() > rho_min) & (m.padj_ZIKV < fdr) & (m.rho_ZIKV.abs() > rho_min)
    
    # Concordant directions
    conv = m[sig_both & (np.sign(m.rho_DENV) == np.sign(m.rho_ZIKV))].copy()
    
    conv["direction"] = np.where(conv.rho_DENV > 0, "up_with_progression", "down_with_progression")
    conv["mean_abs_rho"] = (conv.rho_DENV.abs() + conv.rho_ZIKV.abs()) / 2
    conv = conv.sort_values("mean_abs_rho", ascending=False)
    
    # Map to core DE signature if available
    if core_file.exists():
        de_core = set(pd.read_csv(core_file)["symbol"].str.upper())
        conv["in_DE_core"] = conv.symbol.str.upper().isin(de_core)
    else:
        logger.warning(f"Core signature file missing: {core_file}. Skipping 'in_DE_core' annotation.")
        conv["in_DE_core"] = False
        
    n_up = (conv.direction == "up_with_progression").sum()
    n_down = (conv.direction == "down_with_progression").sum()
    logger.info(f"Identified {len(conv)} convergent progression genes ({n_up} up, {n_down} down).")
    
    return conv

def plot_trends(conv: pd.DataFrame, adata: sc.AnnData, out_plot: Path, logger: logging.Logger) -> None:
    """Plot trends for the top 9 genes."""
    if conv.empty:
        logger.warning("No convergent genes to plot.")
        return
        
    top = conv.head(9)
    symbols = adata.var["symbol"].astype(str)
    s2i = {s.upper(): i for i, s in enumerate(symbols)}
    
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))
    axes = axes.ravel()
    
    for ax_, (_, r) in zip(axes, top.iterrows()):
        gid = s2i.get(str(r.symbol).upper())
        if gid is None:
            continue
            
        for v, c in [("DENV", "#c0392b"), ("ZIKV", "#2471a3")]:
            sub = adata[(adata.obs.condition == v).to_numpy()]
            x = sub.obs.log_viral_load.to_numpy()
            
            e = sub[:, gid].X
            e = np.asarray(e.todense()).ravel() if hasattr(e, "todense") else np.asarray(e).ravel()
            
            o = np.argsort(x)
            df = pd.DataFrame({"x": x[o], "e": e[o]})
            w = max(20, len(df) // 15)
            
            ax_.scatter(df.x, df.e, s=3, alpha=0.12, c=c)
            ax_.plot(df.x, df.e.rolling(w, min_periods=1, center=True).mean(), c=c, lw=2, label=v)
            
        ax_.set_title(f"{r.symbol} (ρD={r.rho_DENV:.2f}, ρZ={r.rho_ZIKV:.2f})", fontsize=9)
        ax_.set_xlabel("log viral load", fontsize=8)
        ax_.set_ylabel("Expression", fontsize=8)
        ax_.legend(fontsize=6)
        
    fig.tight_layout()
    try:
        fig.savefig(out_plot, dpi=150)
        logger.info(f"Saved trends plot to {out_plot}")
    except Exception as e:
        logger.error(f"Failed to save plot: {e}")
    finally:
        plt.close(fig)

def main():
    """Main execution entry point."""
    args = parse_args()
    logger = setup_logger(args.log_level)
    
    logger.info("Starting Step 05: Progression Signature")
    
    if not args.in_file.exists():
        logger.error(f"Input file not found: {args.in_file}")
        sys.exit(1)
        
    logger.info(f"Loading processed AnnData from {args.in_file}...")
    try:
        adata = sc.read_h5ad(args.in_file)
    except Exception as e:
        logger.error(f"Failed to read input AnnData: {e}")
        sys.exit(1)

    if "log_viral_load" not in adata.obs:
        adata.obs["log_viral_load"] = np.log1p(adata.obs["viral_load"].astype(float))

    try:
        d = run_correlation(adata, "DENV", args.fdr, args.rho_min, logger)
        z = run_correlation(adata, "ZIKV", args.fdr, args.rho_min, logger)
    except Exception as e:
        logger.error(f"Correlation calculation failed: {e}")
        sys.exit(1)

    conv = build_progression_signature(d, z, args.core_file, args.fdr, args.rho_min, logger)
    
    try:
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        conv.to_csv(args.out_csv, index=False)
    except Exception as e:
        logger.error(f"Failed to save signature CSV: {e}")
        sys.exit(1)

    plot_trends(conv, adata, args.out_plot, logger)
    
    logger.info("Step 05 completed successfully.")

if __name__ == "__main__":
    main()
