import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Nature SR Aesthetics
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

DATA_PATH = Path("E:/Zika_Enrichment/cell_state_analysis/adata_annotated.h5ad")
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_03_Fig_3_scRNA_CellStates")

def main():
    pass  # Execution logging removed for final release
    adata = sc.read_h5ad(DATA_PATH)
    
    # We want a heatmap (or dotplot) of the top marker genes for infection_state
    # Check if 'infection_state' exists
    if 'infection_state' not in adata.obs:
        pass  # Execution logging removed for final release
        return
        
    pass  # Execution logging removed for final release
    # Compute rank_genes_groups if not already present or to be sure
    sc.tl.rank_genes_groups(adata, 'infection_state', method='wilcoxon', key_added='rank_genes_infection')
    
    # 180mm x 140mm -> 7.08 x 5.51 inches
    fig, ax = plt.subplots(figsize=(7.08, 5.51))
    
    # sc.pl.rank_genes_groups_matrixplot or dotplot is very standard for single cell
    # Let's use a matrixplot with coolwarm
    pass  # Execution logging removed for final release
    sc.pl.rank_genes_groups_matrixplot(
        adata, 
        n_genes=8, 
        key='rank_genes_infection', 
        cmap='RdBu_r', 
        standard_scale='var', 
        colorbar_title='Scaled expression',
        show=False,
        ax=ax
    )
    
    # Remove some spines from the matrix plot axes if accessible
    for a in fig.axes:
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
    
    out_file = OUT_DIR / "scRNA_Heatmap_Markers.tiff"
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches='tight', format='tiff')
    plt.close()
    pass  # Execution logging removed for final release

if __name__ == "__main__":
    main()
