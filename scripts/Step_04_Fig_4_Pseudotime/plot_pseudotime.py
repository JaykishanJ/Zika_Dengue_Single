import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Nature SR Aesthetics
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

import seaborn as sns
from pathlib import Path
import scipy.sparse

DATA_PATH = Path("E:/Zika_Enrichment/cell_state_analysis/adata_trajectory.h5ad")
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_04_Fig_4_Pseudotime")

def plot_cell_cycle(adata):
    # Filter out Unknown
    df = adata.obs[adata.obs['infection_state'] != 'Unknown'].copy()
    
    # Calculate proportions
    props = df.groupby(['infection_state', 'phase']).size().unstack(fill_value=0)
    props = props.div(props.sum(axis=1), axis=0) * 100
    
    # Ensure standard order
    order = ['Bystander_DENV', 'Low_DENV', 'High_DENV', 'Bystander_ZIKV', 'Low_ZIKV', 'High_ZIKV']
    props = props.reindex([o for o in order if o in props.index])
    
    # 180x150mm -> 7.08x5.90 inches
    fig, ax = plt.subplots(figsize=(7.08, 5.90))
    props.plot(kind='bar', stacked=True, ax=ax, color=['#66c2a5', '#fc8d62', '#8da0cb'])
    plt.title('Cell Cycle Phase by Infection State', fontsize=14)
    plt.ylabel('Percentage of Cells (%)', fontsize=12)
    plt.xlabel('Infection State', fontsize=12)
    # Legend on the right side
    plt.legend(title='Cell Cycle Phase', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(OUT_DIR / "Cell_Cycle_Distribution.tiff", dpi=300, format='tiff', bbox_inches='tight')
    plt.close()

def plot_gene_trends(adata, genes):
    # Ensure genes exist in var_names
    valid_genes = [g for g in genes if g in adata.var_names]
    if not valid_genes:
        pass  # Execution logging removed for final release
        return
        
    for virus in ['ZIKV', 'DENV']:
        adata_sub = adata[adata.obs['infection_state'].str.contains(virus)].copy()
        adata_sub = adata_sub[~adata_sub.obs['pseudotime'].isna()]
        
        if len(adata_sub) < 10: continue
            
        fig, axes = plt.subplots(len(valid_genes), 1, figsize=(8, 3 * len(valid_genes)), sharex=True)
        if len(valid_genes) == 1: axes = [axes]
        
        # Sort cells by pseudotime
        order = np.argsort(adata_sub.obs['pseudotime'].values)
        pt = adata_sub.obs['pseudotime'].values[order]
        
        for ax, g in zip(axes, valid_genes):
            idx = adata_sub.var_names.get_loc(g)
            if scipy.sparse.issparse(adata_sub.X):
                expr = adata_sub.X[:, idx].toarray().flatten()
            else:
                expr = adata_sub.X[:, idx]
            expr = expr[order]
            
            # Smoothing (moving average)
            window = max(len(expr) // 20, 5)
            expr_smooth = pd.Series(expr).rolling(window=window, min_periods=1).mean().values
            
            ax.scatter(pt, expr, alpha=0.2, s=10, color='gray')
            ax.plot(pt, expr_smooth, color='red', linewidth=2)
            ax.set_ylabel(f"{g} Expression")
            ax.set_title(f"{g} Dynamics during {virus} Infection")
            
        axes[-1].set_xlabel("Infection Pseudotime")
        plt.tight_layout()
        plt.savefig(OUT_DIR / f"Pseudotime_Trends_{virus}.png", dpi=300)
        plt.close()

def main():
    pass  # Execution logging removed for final release
    adata = sc.read_h5ad(DATA_PATH)
    
    pass  # Execution logging removed for final release
    plot_cell_cycle(adata)
    
    # Try using symbols if they exist, otherwise try Ensembl mapping
    # Actually, we mapped symbols to ENSG in the first step! The var_names are ENSG.
    # Let's map our target genes back to ENSG just to extract them, OR just use the ENSG IDs.
    target_symbols = ['CCL4', 'MT2A', 'ISG15', 'CXCL10', 'MX1']
    import mygene
    mg = mygene.MyGeneInfo()
    res = mg.querymany(target_symbols, scopes='symbol', fields='ensembl.gene', species='human', verbose=False)
    ensg_ids = []
    symbol_map = {}
    for r in res:
        if 'ensembl' in r:
            ens = r['ensembl']
            if isinstance(ens, list):
                e_id = ens[0].get('gene')
            elif isinstance(ens, dict):
                e_id = ens.get('gene')
            ensg_ids.append(e_id)
            symbol_map[e_id] = r['query']
            
    pass  # Execution logging removed for final release
    valid_ensg = [e for e in ensg_ids if e in adata.var_names]
    
    for virus in ['ZIKV', 'DENV']:
        adata_sub = adata[adata.obs['infection_state'].str.contains(virus)].copy()
        adata_sub = adata_sub[~adata_sub.obs['pseudotime'].isna()]
        if len(adata_sub) < 10: continue
            
        fig, axes = plt.subplots(len(valid_ensg), 1, figsize=(8, 3 * len(valid_ensg)), sharex=True)
        if len(valid_ensg) == 1: axes = [axes]
        
        order = np.argsort(adata_sub.obs['pseudotime'].values)
        pt = adata_sub.obs['pseudotime'].values[order]
        
        for ax, g in zip(axes, valid_ensg):
            idx = adata_sub.var_names.get_loc(g)
            if scipy.sparse.issparse(adata_sub.X):
                expr = adata_sub.X[:, idx].toarray().flatten()
            else:
                expr = adata_sub.X[:, idx]
            expr = expr[order]
            
            window = max(len(expr) // 20, 5)
            expr_smooth = pd.Series(expr).rolling(window=window, min_periods=1).mean().values
            
            sym = symbol_map[g]
            ax.scatter(pt, expr, alpha=0.2, s=10, color='lightgray')
            ax.plot(pt, expr_smooth, color='firebrick', linewidth=2)
            ax.set_ylabel(f"Expression")
            ax.set_title(f"{sym} Dynamics ({virus})")
            
        axes[-1].set_xlabel("Infection Pseudotime")
        plt.tight_layout()
        plt.savefig(OUT_DIR / f"Pseudotime_Trends_{virus}.png", dpi=300)
        plt.close()
        
    pass  # Execution logging removed for final release

if __name__ == "__main__":
    main()
