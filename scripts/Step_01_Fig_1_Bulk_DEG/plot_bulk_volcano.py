import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import mygene
import sys
from pathlib import Path

# Import our new publication style module
sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

BASE_SC = Path("E:/Zika_antigravity/SC_RNA_Zika_DENV")
BULK_DIR = BASE_SC / "Bulk_RNA_seq"
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_01_Fig_1_Bulk_DEG")

# Dengue uses Entrez IDs, Zika uses Symbols in the Gene column.
def plot_volcano(df, title, out_path, is_dengue=False):
    df = df.dropna(subset=['padj', 'log2FoldChange', 'Gene']).copy()
    
    # Calculate -log10(padj)
    # Clip tiny p-values to avoid inf
    df['padj'] = df['padj'].clip(lower=1e-300)
    df['nlog10'] = -np.log10(df['padj'])
    
    # Classify
    df['Sig'] = 'Not Significant'
    df.loc[(df['padj'] < 0.05) & (df['log2FoldChange'] > 1), 'Sig'] = 'Up'
    df.loc[(df['padj'] < 0.05) & (df['log2FoldChange'] < -1), 'Sig'] = 'Down'
    
    if is_dengue:
        # Convert Entrez to Symbol for top labels
        top_genes = df[df['Sig'] != 'Not Significant'].sort_values('nlog10', ascending=False).head(10)
        mg = mygene.MyGeneInfo()
        res = mg.querymany(top_genes['Gene'].astype(str).tolist(), scopes='entrezgene', fields='symbol', species='human', verbose=False)
        symbol_map = {r['query']: r.get('symbol', r['query']) for r in res}
        df['Symbol'] = df['Gene'].astype(str).map(symbol_map).fillna(df['Gene'])
    else:
        df['Symbol'] = df['Gene']

    # Colors
    palette = {'Up': '#E63946', 'Down': '#457B9D', 'Not Significant': '#E0E0E0'}
    
    fig, ax = plt.subplots()
    
    # Plot non-sig first to push to background
    sns.scatterplot(data=df[df['Sig']=='Not Significant'], x='log2FoldChange', y='nlog10', color=palette['Not Significant'], s=15, alpha=0.5, edgecolor=None, ax=ax)
    sns.scatterplot(data=df[df['Sig']=='Up'], x='log2FoldChange', y='nlog10', color=palette['Up'], s=25, alpha=0.8, edgecolor='black', linewidth=0.2, ax=ax)
    sns.scatterplot(data=df[df['Sig']=='Down'], x='log2FoldChange', y='nlog10', color=palette['Down'], s=25, alpha=0.8, edgecolor='black', linewidth=0.2, ax=ax)
    
    # Lines
    ax.axvline(1, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline(-1, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # Annotate top 5 Up and top 5 Down
    top_up = df[df['Sig'] == 'Up'].sort_values('nlog10', ascending=False).head(5)
    top_down = df[df['Sig'] == 'Down'].sort_values('nlog10', ascending=False).head(5)
    
    texts = []
    from adjustText import adjust_text
    
    for _, row in pd.concat([top_up, top_down]).iterrows():
        texts.append(ax.text(row['log2FoldChange'], row['nlog10'], row['Symbol'], fontsize=9, fontweight='bold'))
    
    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
        
    ax.set_title(title, pad=15, fontweight='bold')
    ax.set_xlabel(r'$\log_{2}$ Fold Change')
    ax.set_ylabel(r'$-\log_{10}$ adjusted p-value')
    
    # Custom Legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], marker='o', color='w', markerfacecolor=palette['Up'], markersize=8),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=palette['Down'], markersize=8),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=palette['Not Significant'], markersize=8)]
    ax.legend(custom_lines, ['Upregulated', 'Downregulated', 'Not Significant'], loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
    
    # Adjust layout so legend isn't cut off
    plt.subplots_adjust(right=0.75)
    pub_style.save_pub_fig(fig, out_path, width_mm=89, height_mm=75)
    plt.close()

def main():
    pass  # Execution logging removed for final release
    zikv_path = BULK_DIR / "GSE118305/Results/Tables/ZIKV_positive_vs_Mock_All_Genes.csv"
    if zikv_path.exists():
        df_zikv = pd.read_csv(zikv_path)
        plot_volcano(df_zikv, "ZIKV vs Mock (Bulk RNA-seq)", OUT_DIR / "Bulk_Volcano_ZIKV.tiff", is_dengue=False)
    
    pass  # Execution logging removed for final release
    denv_path = BULK_DIR / "GSE279208/Results/Tables/Dengue_vs_healthy_All_Genes.csv"
    if denv_path.exists():
        df_denv = pd.read_csv(denv_path)
        plot_volcano(df_denv, "DENV vs Healthy (Bulk RNA-seq)", OUT_DIR / "Bulk_Volcano_DENV.tiff", is_dengue=True)

if __name__ == "__main__":
    main()
