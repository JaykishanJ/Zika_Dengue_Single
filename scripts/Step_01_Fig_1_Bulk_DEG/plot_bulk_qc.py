import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from pathlib import Path

# Import our new publication style module
sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

BASE_SC = Path("E:/Zika_antigravity/SC_RNA_Zika_DENV")
BULK_DIR = BASE_SC / "Bulk_RNA_seq"
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_01_Fig_1_Bulk_DEG")

def plot_ma(df, title, out_path):
    df = df.dropna(subset=['baseMean', 'log2FoldChange', 'padj']).copy()
    
    # Classify
    df['Sig'] = 'Not Significant'
    df.loc[(df['padj'] < 0.05) & (df['log2FoldChange'] > 1), 'Sig'] = 'Up'
    df.loc[(df['padj'] < 0.05) & (df['log2FoldChange'] < -1), 'Sig'] = 'Down'

    # Colors
    palette = {'Up': '#E63946', 'Down': '#457B9D', 'Not Significant': '#E0E0E0'}
    
    fig, ax = plt.subplots()
    
    # Plot non-sig first
    sns.scatterplot(data=df[df['Sig']=='Not Significant'], x='baseMean', y='log2FoldChange', color=palette['Not Significant'], s=10, alpha=0.3, edgecolor=None, ax=ax)
    sns.scatterplot(data=df[df['Sig']=='Up'], x='baseMean', y='log2FoldChange', color=palette['Up'], s=15, alpha=0.8, edgecolor='black', linewidth=0.2, ax=ax)
    sns.scatterplot(data=df[df['Sig']=='Down'], x='baseMean', y='log2FoldChange', color=palette['Down'], s=15, alpha=0.8, edgecolor='black', linewidth=0.2, ax=ax)
    
    # Log scale x
    ax.set_xscale('log')
    
    # Lines
    ax.axhline(0, color='black', linewidth=1.2)
    ax.axhline(1, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axhline(-1, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
        
    ax.set_title(title, pad=10, fontweight='bold')
    ax.set_xlabel('Mean Normalized Expression (baseMean)')
    ax.set_ylabel(r'$\log_{2}$ Fold Change')
    
    pub_style.save_pub_fig(fig, out_path, width_mm=89, height_mm=75)
    plt.close()

def main():
    pass  # Execution logging removed for final release
    zikv_path = BULK_DIR / "GSE118305/Results/Tables/ZIKV_positive_vs_Mock_All_Genes.csv"
    if zikv_path.exists():
        df_zikv = pd.read_csv(zikv_path)
        plot_ma(df_zikv, "ZIKV MA Plot", OUT_DIR / "Bulk_MA_ZIKV.tiff")
    
    pass  # Execution logging removed for final release
    denv_path = BULK_DIR / "GSE279208/Results/Tables/Dengue_vs_healthy_All_Genes.csv"
    if denv_path.exists():
        df_denv = pd.read_csv(denv_path)
        plot_ma(df_denv, "DENV MA Plot", OUT_DIR / "Bulk_MA_DENV.tiff")

if __name__ == "__main__":
    main()
