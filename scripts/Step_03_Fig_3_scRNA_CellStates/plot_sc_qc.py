import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from pathlib import Path

sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

BASE_SC = Path("E:/Zika_antigravity/SC_RNA_Zika_DENV")
ADATA_PATH = BASE_SC / "01_processed_data/anndata_objects/adata_processed.h5ad"
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_03_Fig_3_scRNA_CellStates")

def main():
    if not ADATA_PATH.exists():
        pass  # Execution logging removed for final release
        return
        
    pass  # Execution logging removed for final release
    adata = sc.read_h5ad(ADATA_PATH)
    
    # 1. Violin plots (n_genes_by_counts, total_counts, pct_counts_mt)
    if 'n_genes_by_counts' in adata.obs and 'pct_counts_mt' in adata.obs:
        pass  # Execution logging removed for final release
        fig, axes = plt.subplots(1, 3, figsize=(7.08, 3.54)) # 180x90mm
        
        sns.violinplot(y=adata.obs['n_genes_by_counts'], ax=axes[0], color='#A8DADC', inner='quartile')
        axes[0].set_ylabel('Detected genes')
        
        sns.violinplot(y=adata.obs['total_counts'], ax=axes[1], color='#457B9D', inner='quartile')
        axes[1].set_ylabel('Total UMIs')
        
        sns.violinplot(y=adata.obs['pct_counts_mt'], ax=axes[2], color='#E63946', inner='quartile')
        axes[2].set_ylabel('Mitochondrial %')
        
        for ax in axes:
            ax.set_xticks([])
            ax.set_xlabel('')
            
        fig.suptitle('Single-cell quality control metrics', fontweight='bold', y=1.05)
        
        pub_style.save_pub_fig(fig, OUT_DIR / "scRNA_QC_Violin.tiff", width_mm=180, height_mm=90)
        plt.close()
        
    # 2. PCA Elbow
    if 'pca' in adata.uns and 'variance_ratio' in adata.uns['pca']:
        pass  # Execution logging removed for final release
        var_ratio = adata.uns['pca']['variance_ratio']
        fig, ax = plt.subplots()
        ax.plot(range(1, len(var_ratio)+1), var_ratio, marker='o', markersize=4, linestyle='-', color='black')
        ax.set_title("PCA Variance Explained", pad=10, fontweight='bold')
        ax.set_xlabel("Principal Component")
        ax.set_ylabel("Fraction of variance explained")
        
        pub_style.save_pub_fig(fig, OUT_DIR / "scRNA_PCA_Elbow.tiff", width_mm=89, height_mm=75)
        plt.close()
        
    # 3. UMAP by dataset/batch
    if 'condition' in adata.obs and 'X_umap' in adata.obsm:
        pass  # Execution logging removed for final release
        fig, ax = plt.subplots(figsize=(4.72, 3.93)) # 120x100mm
        sc.pl.umap(adata, color='condition', ax=ax, show=False, frameon=False, title="UMAP by Infection Condition")
        
        pub_style.save_pub_fig(fig, OUT_DIR / "UMAP_Batch_Dataset.tiff", width_mm=120, height_mm=100)
        plt.close()

if __name__ == "__main__":
    main()
