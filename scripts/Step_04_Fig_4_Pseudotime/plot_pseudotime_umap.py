import scanpy as sc
import matplotlib.pyplot as plt
import sys
from pathlib import Path

sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

DATA_PATH = Path("E:/Zika_Enrichment/cell_state_analysis/adata_trajectory.h5ad")
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_04_Fig_4_Pseudotime")

def main():
    if not DATA_PATH.exists():
        pass  # Execution logging removed for final release
        return
        
    pass  # Execution logging removed for final release
    adata = sc.read_h5ad(DATA_PATH)
    
    if 'pseudotime' in adata.obs:
        pass  # Execution logging removed for final release
        fig, ax = plt.subplots(figsize=(4.72, 3.93)) # 120x100mm
        sc.pl.umap(adata, color='pseudotime', ax=ax, show=False, frameon=False, 
                   title="Pseudotime Ordering", color_map="viridis")
        
        # In scanpy the colorbar is added to the figure. We can let scanpy handle it, 
        # or we can tweak it. scanpy's default is usually okay.
        
        pub_style.save_pub_fig(fig, OUT_DIR / "UMAP_Pseudotime.tiff", width_mm=120, height_mm=100)
        plt.close()
    else:
        pass  # Execution logging removed for final release

if __name__ == "__main__":
    main()
