import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Nature SR Aesthetics
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

DATA_PATH = Path("E:/Zika_antigravity/SC_RNA_Zika_DENV/03_results/phase4_All_DEG/All_Significant_DEGs.csv")
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_03_Fig_3_scRNA_CellStates")

def plot_volcano(df, disease, out_path):
    fc_col = f'{disease}_log2FC'
    pval_col = f'{disease}_padj'
    
    if fc_col not in df.columns or pval_col not in df.columns:
        pass  # Execution logging removed for final release
        return
        
    sub_df = df.dropna(subset=[fc_col, pval_col, 'symbol']).copy()
    
    # Clip tiny p-values to avoid inf
    sub_df[pval_col] = sub_df[pval_col].clip(lower=1e-300)
    sub_df['nlog10'] = -np.log10(sub_df[pval_col])
    
    # Classify
    sub_df['Sig'] = 'Not Significant'
    sub_df.loc[(sub_df[pval_col] < 0.05) & (sub_df[fc_col] > 1), 'Sig'] = 'Up'
    sub_df.loc[(sub_df[pval_col] < 0.05) & (sub_df[fc_col] < -1), 'Sig'] = 'Down'
    
    # Colors
    palette = {'Up': '#E63946', 'Down': '#457B9D', 'Not Significant': '#E0E0E0'}
    
    fig, ax = plt.subplots(figsize=(7.08, 6.30)) # 180x160mm
    
    # Plot non-sig first to push to background
    sns.scatterplot(data=sub_df[sub_df['Sig']=='Not Significant'], x=fc_col, y='nlog10', color=palette['Not Significant'], s=15, alpha=0.5, edgecolor=None, ax=ax)
    sns.scatterplot(data=sub_df[sub_df['Sig']=='Up'], x=fc_col, y='nlog10', color=palette['Up'], s=25, alpha=0.8, edgecolor='black', linewidth=0.2, ax=ax)
    sns.scatterplot(data=sub_df[sub_df['Sig']=='Down'], x=fc_col, y='nlog10', color=palette['Down'], s=25, alpha=0.8, edgecolor='black', linewidth=0.2, ax=ax)
    
    # Lines
    ax.axvline(1, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline(-1, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # Annotate top 5 Up and top 5 Down
    top_up = sub_df[sub_df['Sig'] == 'Up'].sort_values('nlog10', ascending=False).head(5)
    top_down = sub_df[sub_df['Sig'] == 'Down'].sort_values('nlog10', ascending=False).head(5)
    
    texts = []
    from adjustText import adjust_text
    
    for _, row in pd.concat([top_up, top_down]).iterrows():
        texts.append(ax.text(row[fc_col], row['nlog10'], str(row['symbol']), fontsize=9, fontweight='bold'))
    
    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
        
    ax.set_title(f"High {disease} vs Bystander (scRNA-seq)", pad=15, fontweight='bold')
    ax.set_xlabel(r'$\log_{2}$ Fold Change')
    ax.set_ylabel(r'$-\log_{10}$ adjusted p-value')
    
    # Custom Legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], marker='o', color='w', markerfacecolor=palette['Up'], markersize=8),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=palette['Down'], markersize=8),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=palette['Not Significant'], markersize=8)]
    ax.legend(custom_lines, ['Upregulated', 'Downregulated', 'Not Significant'], loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
    
    plt.tight_layout()
    plt.subplots_adjust(right=0.75)
    plt.savefig(out_path, dpi=300, format='tiff', bbox_inches='tight')
    plt.close()
    pass  # Execution logging removed for final release

def main():
    pass  # Execution logging removed for final release
    df = pd.read_csv(DATA_PATH)
    
    pass  # Execution logging removed for final release
    plot_volcano(df, "ZIKV", OUT_DIR / "scRNA_Volcano_ZIKV.tiff")
    
    pass  # Execution logging removed for final release
    plot_volcano(df, "DENV", OUT_DIR / "scRNA_Volcano_DENV.tiff")

if __name__ == "__main__":
    main()
