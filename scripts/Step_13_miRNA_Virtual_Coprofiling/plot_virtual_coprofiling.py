import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import mygene

# AESTHETICS
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.linewidth'] = 1.2

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_13_miRNA_Virtual_Coprofiling")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def main():
    pass  # Execution logging removed for final release
    adata = sc.read_h5ad("E:/Zika_Enrichment/cell_state_analysis/adata_annotated.h5ad")

    # Define target groups based on previous network hub knowledge
    targets_3135b = ['AKNA', 'IFIT3', 'KCNK6', 'KLHDC7B', 'S1PR3', 'SLC35E4', 'STAR']
    targets_491_5p = ['AKNA', 'ERF', 'FOLR1', 'KRT15', 'PSD4']
    targets_140_3p = ['ABCD3', 'ERF', 'PI4K2B', 'S1PR3', 'TAP1']

    # Map symbols to Ensembl IDs since adata uses Ensembl
    pass  # Execution logging removed for final release
    mg = mygene.MyGeneInfo()
    all_symbols = list(set(targets_3135b + targets_491_5p + targets_140_3p))
    res = mg.querymany(all_symbols, scopes='symbol', fields='ensembl.gene', species='human', verbose=False)
    
    symbol_to_ensembl = {}
    for r in res:
        if 'ensembl' in r:
            ens = r['ensembl']
            if isinstance(ens, list):
                symbol_to_ensembl[r['query']] = ens[0]['gene']
            else:
                symbol_to_ensembl[r['query']] = ens['gene']

    # Convert lists to Ensembl IDs
    targets_3135b = [symbol_to_ensembl.get(g) for g in targets_3135b if symbol_to_ensembl.get(g)]
    targets_491_5p = [symbol_to_ensembl.get(g) for g in targets_491_5p if symbol_to_ensembl.get(g)]
    targets_140_3p = [symbol_to_ensembl.get(g) for g in targets_140_3p if symbol_to_ensembl.get(g)]

    # Filter valid genes in adata
    targets_3135b = [g for g in targets_3135b if g in adata.var_names]
    targets_491_5p = [g for g in targets_491_5p if g in adata.var_names]
    targets_140_3p = [g for g in targets_140_3p if g in adata.var_names]
    
    pass  # Execution logging removed for final release

    # 2. Score cells based on target expression (Module Score)
    pass  # Execution logging removed for final release
    if targets_3135b: sc.tl.score_genes(adata, gene_list=targets_3135b, score_name='miR_3135b_Target_Score')
    if targets_491_5p: sc.tl.score_genes(adata, gene_list=targets_491_5p, score_name='miR_491_5p_Target_Score')
    if targets_140_3p: sc.tl.score_genes(adata, gene_list=targets_140_3p, score_name='miR_140_3p_Target_Score')

    # Filter to only relevant states (High vs Bystander)
    adata_sub = adata[adata.obs['infection_state'].isin(['Bystander_ZIKV', 'High_ZIKV', 'Bystander_DENV', 'High_DENV'])].copy()

    # Ensure categorical order
    adata_sub.obs['infection_state'] = adata_sub.obs['infection_state'].astype('category')
    adata_sub.obs['infection_state'] = adata_sub.obs['infection_state'].cat.reorder_categories(['Bystander_DENV', 'High_DENV', 'Bystander_ZIKV', 'High_ZIKV'])

    # Prepare df for plotting
    scores_to_plot = []
    if targets_3135b: scores_to_plot.append('miR_3135b_Target_Score')
    if targets_491_5p: scores_to_plot.append('miR_491_5p_Target_Score')
    if targets_140_3p: scores_to_plot.append('miR_140_3p_Target_Score')
    
    df = adata_sub.obs[['infection_state'] + scores_to_plot].copy()

    # Plot Violin plots
    pass  # Execution logging removed for final release
    fig, axes = plt.subplots(1, len(scores_to_plot), figsize=(5 * len(scores_to_plot), 5))
    if len(scores_to_plot) == 1: axes = [axes]

    palette = ['#4e79a7', '#f28e2b', '#e15759', '#76b7b2']
    
    for i, score in enumerate(scores_to_plot):
        sns.violinplot(data=df, x='infection_state', y=score, ax=axes[i], palette=palette)
        axes[i].set_title(f"{score.replace('_Target_Score', '')} Target Expression", fontweight='bold')
        axes[i].set_ylabel('Target Gene Module Score')
        axes[i].set_xticklabels(['Bys_DENV', 'High_DENV', 'Bys_ZIKV', 'High_ZIKV'], rotation=45)

    plt.tight_layout()
    out_violin = OUT_DIR / "Virtual_Coprofiling_Violin.tiff"
    plt.savefig(out_violin, dpi=300)
    plt.close()
    pass  # Execution logging removed for final release

    # Aggregated dotplot for specific targets
    all_targets = list(set(targets_3135b + targets_491_5p + targets_140_3p))
    if all_targets:
        pass  # Execution logging removed for final release
        sc.pl.dotplot(adata_sub, var_names=all_targets, groupby='infection_state', show=False, cmap='Reds', standard_scale='var')
        out_dot = OUT_DIR / "Virtual_Coprofiling_Dotplot.tiff"
        plt.savefig(out_dot, dpi=300, bbox_inches='tight')
        plt.close()
        pass  # Execution logging removed for final release

    pass  # Execution logging removed for final release

if __name__ == "__main__":
    main()
