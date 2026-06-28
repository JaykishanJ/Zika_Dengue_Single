import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

# Import our new publication style module
sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

import seaborn as sns
from pathlib import Path
import os
import scipy.sparse
import mygene

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_05_Fig_5_CellChat")

pass  # Execution logging removed for final release
core_pairs = [
    ('CXCL10', 'CXCR3'), ('CCL4', 'CCR5'), ('CCL5', 'CCR5'), ('CCL5', 'CCR1'), 
    ('CCL5', 'CCR3'), ('IFNB1', 'IFNAR1'), ('IFNB1', 'IFNAR2'), ('IFNA1', 'IFNAR1'), 
    ('IFNA2', 'IFNAR2'), ('IL6', 'IL6R'), ('IL6', 'IL6ST'), ('TNF', 'TNFRSF1A'),
    ('TNF', 'TNFRSF1B'), ('CXCL8', 'CXCR1'), ('CXCL8', 'CXCR2'), ('IL1B', 'IL1R1'), 
    ('IL1A', 'IL1R1'), ('FASLG', 'FAS'), ('VEGFA', 'FLT1'), ('VEGFA', 'KDR'), 
    ('CSF1', 'CSF1R'), ('CCL2', 'CCR2'), ('CXCL9', 'CXCR3'), ('CXCL11', 'CXCR3'),
    ('IL10', 'IL10RA'), ('IL10', 'IL10RB'), ('IL15', 'IL15RA'), ('IL18', 'IL18R1'),
    ('IL33', 'IL1RL1'), ('TGFB1', 'TGFBR1'), ('TGFB1', 'TGFBR2'), ('CCL3', 'CCR1'),
    ('CCL3', 'CCR5'), ('CX3CL1', 'CX3CR1'), ('MIF', 'CD74'), ('MIF', 'CXCR4'),
    ('MIF', 'CXCR2'), ('TNFSF10', 'TNFRSF10A'), ('TNFSF10', 'TNFRSF10B'),
    ('TNFSF13B', 'TNFRSF13B'), ('TNFSF13B', 'TNFRSF17'), ('TNFSF14', 'TNFRSF14'),
    ('CD40LG', 'CD40'), ('NAMPT', 'ITGA5'), ('APP', 'CD74'), ('CALR', 'LRP1')
]
df_lr = pd.DataFrame(core_pairs, columns=['Ligand', 'Receptor'])
pass  # Execution logging removed for final release

pass  # Execution logging removed for final release
adata = sc.read_h5ad("E:/Zika_Enrichment/cell_state_analysis/adata_trajectory.h5ad")

pass  # Execution logging removed for final release
mg = mygene.MyGeneInfo()
res = mg.querymany(adata.var_names.tolist(), scopes='ensembl.gene', fields='symbol', species='human', verbose=False)
ensg_to_symbol = {r['query']: r.get('symbol', r['query']) for r in res}
adata.var['Symbol'] = [ensg_to_symbol.get(g, g) for g in adata.var_names]

pass  # Execution logging removed for final release
dense_X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
df_expr = pd.DataFrame(dense_X, columns=adata.var['Symbol'].values)
df_expr['infection_state'] = adata.obs['infection_state'].values

mean_expr = df_expr.groupby('infection_state', observed=False).mean()

pass  # Execution logging removed for final release
results = []

for virus in ['DENV', 'ZIKV']:
    high_state = f'High_{virus}'
    bystander_state = f'Bystander_{virus}'
    
    if high_state not in mean_expr.index or bystander_state not in mean_expr.index:
        continue
        
    for _, row in df_lr.iterrows():
        ligand = row['Ligand']
        receptor = row['Receptor']
        
        if ligand in mean_expr.columns and receptor in mean_expr.columns:
            l_expr = mean_expr.loc[high_state, ligand]
            r_expr = mean_expr.loc[bystander_state, receptor]
            
            score = l_expr * r_expr
            
            # Since data is normalized log(1p)
            if score > 0.001: 
                results.append({
                    'Virus': virus,
                    'Direction': f'{high_state} -> {bystander_state}',
                    'LR_Pair': f'{ligand} - {receptor}',
                    'Score': score,
                    'Ligand_Expr': l_expr,
                    'Receptor_Expr': r_expr
                })

df_res = pd.DataFrame(results)

if df_res.empty:
    pass  # Execution logging removed for final release
else:
    df_res = df_res.sort_values('Score', ascending=False)
    
    top_pairs = df_res.head(30).copy()
    top_pairs.to_csv(OUT_DIR / "CellChat_Alternative_Top_LR_Pairs.csv", index=False)
    
    pass  # Execution logging removed for final release
    
    top_pairs['Normalized_Score'] = (top_pairs['Score'] - top_pairs['Score'].min()) / (top_pairs['Score'].max() - top_pairs['Score'].min() + 1e-9)
    top_pairs['Normalized_Score'] = top_pairs['Normalized_Score'] * 600 + 50
    
    # 180x170mm -> 7.08x6.69 inches
    fig, ax = plt.subplots(figsize=(7.08, 6.69))
    sns.scatterplot(
        data=top_pairs, 
        x='Direction', 
        y='LR_Pair', 
        size='Normalized_Score', 
        hue='Score', 
        palette='viridis', 
        sizes=(100, 700),
        ax=ax
    )
    plt.title("Paracrine Signaling (High Infection -> Bystander Cells)", fontsize=13, fontweight='bold', pad=20)
    plt.ylabel("Ligand - Receptor Pair", fontsize=11)
    plt.xlabel("Signaling Direction", fontsize=11)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Comm. Score")
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # We adjust layout to make sure legend fits
    plt.subplots_adjust(right=0.7)
    
    pub_style.save_pub_fig(fig, OUT_DIR / "CellChat_Alternative_DotPlot.tiff", width_mm=180, height_mm=170)
    plt.close()
    
    pass  # Execution logging removed for final release
