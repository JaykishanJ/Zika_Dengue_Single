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
import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path
from src.utils import pub_style

pub_style.apply_style()

BASE_SC = config.ROOT_DIR
ADATA_PATH = BASE_SC / "01_processed_data/anndata_objects/adata_processed.h5ad"
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_03_Fig_3_scRNA_CellStates")

def main():
    if not ADATA_PATH.exists():
        pass  # Execution logging removed for final release
        return
        
    pass  # Execution logging removed for final release
    adata = sc.read_h5ad(ADATA_PATH)
    
    if 'infection_state' in adata.obs and 'condition' in adata.obs:
        pass  # Execution logging removed for final release
        
        # Calculate percentages
        counts = adata.obs.groupby(['condition', 'infection_state']).size().unstack(fill_value=0)
        percentages = counts.div(counts.sum(axis=1), axis=0) * 100
        
        fig, ax = plt.subplots()
        
        # We know states are like High_ZIKV, Bystander_ZIKV, etc. 
        # A simple stacked bar plot works best.
        percentages.plot(kind='bar', stacked=True, ax=ax, colormap='Set2', edgecolor='black', linewidth=0.5)
        
        ax.set_title("Relative Abundance of Annotated Cell States", pad=10, fontweight='bold')
        ax.set_ylabel("% Cells")
        ax.set_xlabel("")
        plt.xticks(rotation=0)
        
        # Move legend outside
        ax.legend(title='Infection State', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
        
        # Adjust layout for legend
        plt.subplots_adjust(right=0.7)
        
        pub_style.save_pub_fig(fig, OUT_DIR / "CellState_Composition_Barplot.tiff", width_mm=89, height_mm=75)
        plt.close()

if __name__ == "__main__":
    main()
