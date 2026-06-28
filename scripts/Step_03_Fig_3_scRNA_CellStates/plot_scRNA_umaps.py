import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Nature SR Aesthetics
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# Paths
DATA_PATH = Path("E:/Zika_Enrichment/cell_state_analysis/adata_annotated.h5ad")
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_03_Fig_3_scRNA_CellStates")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def add_umap_arrows(ax):
    """Add Nature SR style UMAP corner arrows"""
    # Define position and length of arrows in axes fraction coordinates
    x0, y0 = 0.05, 0.05
    length = 0.15
    
    # Create an invisible axes that spans the same area to easily draw in relative coordinates
    ax.annotate("", xy=(x0 + length, y0), xytext=(x0, y0),
                xycoords='axes fraction', arrowprops=dict(arrowstyle="->", lw=1.5, color='black'))
    ax.annotate("", xy=(x0, y0 + length), xytext=(x0, y0),
                xycoords='axes fraction', arrowprops=dict(arrowstyle="->", lw=1.5, color='black'))
                
    ax.text(x0 + length/2, y0 - 0.02, "UMAP_1", ha='center', va='top', transform=ax.transAxes, fontsize=10, fontweight='bold')
    ax.text(x0 - 0.02, y0 + length/2, "UMAP_2", ha='right', va='center', rotation=90, transform=ax.transAxes, fontsize=10, fontweight='bold')

def main():
    pass  # Execution logging removed for final release
    adata = sc.read_h5ad(DATA_PATH)
    
    columns_to_plot = []
    if 'infection_state' in adata.obs:
        columns_to_plot.append('infection_state')
    if 'phase' in adata.obs:
        columns_to_plot.append('phase')
        
    for col in columns_to_plot:
        pass  # Execution logging removed for final release
        
        # 180mm x 160mm -> 7.08 x 6.30 inches
        fig, ax = plt.subplots(figsize=(7.08, 6.30))
        
        # Plot with frameon=False to remove spines
        sc.pl.umap(adata, color=col, ax=ax, show=False, legend_loc='right margin', 
                   frameon=False, size=20, alpha=0.8, title=f"UMAP: {col.replace('_', ' ').title()}")
        
        # Add the Nature style UMAP arrows
        add_umap_arrows(ax)
        
        out_file = OUT_DIR / f"UMAP_{col}.tiff"
        plt.tight_layout()
        plt.savefig(out_file, dpi=300, bbox_inches='tight', format='tiff')
        plt.close()
        pass  # Execution logging removed for final release

if __name__ == "__main__":
    main()
