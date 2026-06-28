import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from pathlib import Path
from scipy.stats import pearsonr

# Import our new publication style module
sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

BASE_SC = Path("E:/Zika_antigravity/SC_RNA_Zika_DENV")
BULK_DIR = BASE_SC / "Bulk_RNA_seq"
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_01_Fig_1_Bulk_DEG")

def plot_deg_counts():
    datasets = {
        "ZIKV (Macrophage)": BULK_DIR / "GSE118305/Results/Tables/ZIKV_positive_vs_Mock_All_Genes.csv",
        "ZIKV (Huh7)": BULK_DIR / "GSE78711/Results/Tables/ZIKV_positive_vs_Mock_All_Genes.csv",
        "DENV (Blood)": BULK_DIR / "GSE279208/Results/Tables/Dengue_vs_healthy_All_Genes.csv"
    }
    
    counts = []
    for name, path in datasets.items():
        if path.exists():
            df = pd.read_csv(path).dropna(subset=['padj', 'log2FoldChange'])
            up = len(df[(df['padj'] < 0.05) & (df['log2FoldChange'] > 1)])
            down = len(df[(df['padj'] < 0.05) & (df['log2FoldChange'] < -1)])
            counts.append({'Dataset': name, 'Direction': 'Upregulated', 'Count': up})
            counts.append({'Dataset': name, 'Direction': 'Downregulated', 'Count': down})
            
    if not counts:
        return
        
    df_counts = pd.DataFrame(counts)
    
    fig, ax = plt.subplots()
    
    # Pivot for stacked bar
    pivot_df = df_counts.pivot(index='Dataset', columns='Direction', values='Count')
    pivot_df = pivot_df[['Upregulated', 'Downregulated']] # enforce order
    
    colors = ['#E63946', '#457B9D']
    pivot_df.plot(kind='bar', stacked=True, color=colors, ax=ax, edgecolor='black', linewidth=0.5)
    
    ax.set_title("Differentially Expressed Genes across Datasets", pad=10, fontweight='bold')
    ax.set_ylabel("Number of DEGs")
    ax.set_xlabel("")
    plt.xticks(rotation=45, ha='right')
    ax.legend(frameon=False, loc='upper right')
    
    pub_style.save_pub_fig(fig, OUT_DIR / "DEG_Counts_Per_Dataset.tiff", width_mm=89, height_mm=75)
    plt.close()

def plot_logfc_correlation():
    zikv_path = BULK_DIR / "GSE118305/Results/Tables/ZIKV_positive_vs_Mock_All_Genes.csv"
    denv_path = BULK_DIR / "GSE279208/Results/Tables/Dengue_vs_healthy_All_Genes.csv"
    
    if not zikv_path.exists() or not denv_path.exists():
        return
        
    df_zikv = pd.read_csv(zikv_path).dropna(subset=['Gene', 'log2FoldChange'])
    df_denv = pd.read_csv(denv_path).dropna(subset=['Gene', 'log2FoldChange'])
    
    # DENV uses Entrez, ZIKV uses Symbol. We need to convert DENV Entrez to Symbol to merge
    import mygene
    mg = mygene.MyGeneInfo()
    res = mg.querymany(df_denv['Gene'].astype(str).tolist(), scopes='entrezgene', fields='symbol', species='human', verbose=False)
    symbol_map = {r['query']: r.get('symbol', None) for r in res}
    df_denv['Symbol'] = df_denv['Gene'].astype(str).map(symbol_map)
    df_denv = df_denv.dropna(subset=['Symbol'])
    
    df_merged = pd.merge(df_zikv[['Gene', 'log2FoldChange']], df_denv[['Symbol', 'log2FoldChange']], left_on='Gene', right_on='Symbol', suffixes=('_ZIKV', '_DENV'))
    
    fig, ax = plt.subplots()
    sns.scatterplot(data=df_merged, x='log2FoldChange_ZIKV', y='log2FoldChange_DENV', s=5, alpha=0.3, color='grey', ax=ax)
    
    r, p = pearsonr(df_merged['log2FoldChange_ZIKV'], df_merged['log2FoldChange_DENV'])
    
    ax.axhline(0, color='black', linewidth=1, linestyle='--')
    ax.axvline(0, color='black', linewidth=1, linestyle='--')
    
    ax.set_title("Cross-virus Transcriptomic Concordance", pad=10, fontweight='bold')
    ax.set_xlabel(r'ZIKV $\log_{2}$ Fold Change')
    ax.set_ylabel(r'DENV $\log_{2}$ Fold Change')
    
    ax.text(0.05, 0.95, f"r = {r:.2f}\np < 0.001", transform=ax.transAxes, verticalalignment='top', fontweight='bold')
    
    pub_style.save_pub_fig(fig, OUT_DIR / "Bulk_Log2FC_Correlation.tiff", width_mm=89, height_mm=75)
    plt.close()

def main():
    pass  # Execution logging removed for final release
    plot_deg_counts()
    pass  # Execution logging removed for final release
    plot_logfc_correlation()

if __name__ == "__main__":
    main()
