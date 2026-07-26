import sys
from pathlib import Path
_root_dir = Path(__file__).resolve().parent.parent.parent.parent
if str(_root_dir) not in sys.path:
    sys.path.append(str(_root_dir))
from src import config

import sys
from pathlib import Path

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import sys
from pathlib import Path

# Setup paths
BASE_DIR = config.ROOT_DIR
DATA_PATH = Path("E:/Zika_Enrichment/cell_state_analysis/adata_trajectory.h5ad")
FIG_DIR = BASE_DIR / "publication_figures"

# Add pipeline root to path so we can from src.utils import pub_style
from src.utils import pub_style

def main():
    if not DATA_PATH.exists():
        print(f"Error: Could not find data file at {DATA_PATH}")
        return
        
    # Ensure output directory exists
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load dataset
    print(f"Loading data from {DATA_PATH}...")
    adata = sc.read_h5ad(DATA_PATH)
    
    # Apply publication styling
    print("Applying Nature/OSCA publication styles...")
    pub_style.apply_style()
    
    # Common plotting parameters for Scanpy to match OSCA minimalism
    # frameon=False removes the box/grid
    plot_kwargs = {
        'show': False,
        'frameon': False,
        'legend_loc': 'right margin'
    }
    
    # ---------------------------------------------------------
    # Panel A: Pseudotime
    # ---------------------------------------------------------
    if 'pseudotime' in adata.obs:
        print("Plotting Pseudotime UMAP...")
        fig, ax = plt.subplots(figsize=(4.72, 3.93)) # 120x100mm
        sc.pl.umap(adata, color='pseudotime', ax=ax, title="Pseudotime Trajectory", color_map="viridis", **plot_kwargs)
        pub_style.save_pub_fig(fig, FIG_DIR / "Python_Pseudotime.pdf")
        plt.close()
    
    # ---------------------------------------------------------
    # Panel B: Condition (ZIKV vs DENV vs Mock/Control)
    # ---------------------------------------------------------
    if 'condition' in adata.obs:
        print("Plotting Condition UMAP...")
        fig, ax = plt.subplots(figsize=(4.72, 3.93))
        sc.pl.umap(adata, color='condition', ax=ax, title="Condition Trajectory", **plot_kwargs)
        pub_style.save_pub_fig(fig, FIG_DIR / "Python_Condition.pdf")
        plt.close()
        
    # ---------------------------------------------------------
    # Panel C: Infection State (Bystander vs Infected)
    # ---------------------------------------------------------
    if 'infection_state' in adata.obs:
        print("Plotting Infection State UMAP...")
        fig, ax = plt.subplots(figsize=(4.72, 3.93))
        sc.pl.umap(adata, color='infection_state', ax=ax, title="Infection States", **plot_kwargs)
        pub_style.save_pub_fig(fig, FIG_DIR / "Python_InfectionState.pdf")
        plt.close()

    print("Successfully exported all trajectory figures to publication_figures!")

if __name__ == "__main__":
    main()
