import sys
from pathlib import Path
_root_dir = Path(__file__).resolve().parent.parent.parent.parent
if str(_root_dir) not in sys.path:
    sys.path.append(str(_root_dir))
from src import config

import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import sys
from pathlib import Path
from src.utils import pub_style
pub_style.apply_style()
FIG_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/publication_figures")

import seaborn as sns
import mygene
from pathlib import Path

# Nature SR Aesthetics
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 8

BASE_SC = config.ROOT_DIR
BULK_DIR = BASE_SC / "Step05b_bulk_differential_expression"
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_01_Fig_1_Bulk_DEG")

def get_bulk_df(disease):
    if disease == 'ZIKV':
        df = pd.read_csv(BULK_DIR / "bulk_GSE118305_ZIKV_macrophage_DE.csv")
        df['Symbol'] = df['SYM']
    else:
        df = pd.read_csv(BULK_DIR / "bulk_GSE279208_DENV_blood_DE.csv")
        df = df.dropna(subset=['SYM'])
        mg = mygene.MyGeneInfo()
        res = mg.querymany(df['SYM'].astype(str).tolist(), scopes='entrezgene', fields='symbol', species='human', verbose=False)
        symbol_map = {r['query']: r.get('symbol', r['query']) for r in res}
        df['Symbol'] = df['SYM'].astype(str).map(symbol_map).fillna(df['SYM'])
        
    df = df.dropna(subset=['padj', 'lfc', 'Symbol'])
    df['Symbol'] = df['Symbol'].astype(str).str.upper()
    # Deduplicate by picking the most significant
    df = df.sort_values('padj').drop_duplicates('Symbol')
    return df.set_index('Symbol')

def main():
    pass  # Execution logging removed for final release
    zikv_df = get_bulk_df('ZIKV')
    denv_df = get_bulk_df('DENV')
    
    # Filter for significant in at least one
    zikv_sig = zikv_df[(zikv_df['padj'] < 0.05) & (zikv_df['lfc'].abs() > 1)].index
    denv_sig = denv_df[(denv_df['padj'] < 0.05) & (denv_df['lfc'].abs() > 1)].index
    
    shared_sig = list(set(zikv_sig) & set(denv_sig))
    
    if not shared_sig:
        pass  # Execution logging removed for final release
        return
        
    pass  # Execution logging removed for final release
    
    # Take top 40 for clean visualization, sorted by average absolute log2FC
    combined = pd.DataFrame({'ZIKV_log2FC': zikv_df['lfc'], 'DENV_log2FC': denv_df['lfc']})
    combined = combined.loc[shared_sig].dropna()
    combined['mean_abs_fc'] = combined[['ZIKV_log2FC', 'DENV_log2FC']].abs().mean(axis=1)
    
    top_genes = combined.sort_values('mean_abs_fc', ascending=False).head(40).index
    plot_df = combined.loc[top_genes, ['ZIKV_log2FC', 'DENV_log2FC']]
    plot_df.columns = ['Zika (Bulk)', 'Dengue (Bulk)']
    
    # 180x160mm -> 7.08x6.30 inches
    fig, ax = plt.subplots(figsize=(7.08, 6.30))
    
    # coolwarm or RdBu_r is classic Nature
    sns.heatmap(plot_df, cmap='RdBu_r', center=0, annot=False, fmt=".1f", linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'log2 Fold Change', 'shrink': 0.8}, ax=ax)
                
    plt.title("Top Shared DEGs Across Infections", pad=15, fontweight='bold')
    plt.xlabel("")
    plt.ylabel("")
    
    plt.tight_layout()
    pub_style.save_pub_fig(plt.gcf(), OUT_DIR / "Bulk_Heatmap_TopShared.tiff", dpi=300, format='tiff')
    plt.close()
    pass  # Execution logging removed for final release

if __name__ == "__main__":
    main()
