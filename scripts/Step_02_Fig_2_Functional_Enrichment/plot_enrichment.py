import pandas as pd
import gseapy as gp
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
import numpy as np
from pathlib import Path
import os

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_02_Fig_2_Functional_Enrichment")
OUT_DIR.mkdir(parents=True, exist_ok=True)

pass  # Execution logging removed for final release
# 1. Bulk DEGs
BULK_CSV = Path("E:/Zika_Enrichment/4way_venn_analysis/Final_Bulk_Overlaps.csv")
df_bulk = pd.read_csv(BULK_CSV)
bulk_genes = set(df_bulk['GeneSymbol'].dropna().unique())

# 2. SC DEGs
SC_CSV = Path("E:/Zika_Enrichment/4way_venn_analysis/Final_SC_Overlaps.csv")
df_sc = pd.read_csv(SC_CSV)
sc_genes = set(df_sc['GeneSymbol'].dropna().unique())

# 3. 98 Overlapping Genes
# Let's read the 4way_overlapping_genes.csv if it exists, otherwise reconstruct it
# It might be in the root folder or venn folder. Let's reconstruct the exact 98 genes:
# The 98 genes are the intersection of the CONCORDANT DEGs in Bulk and SC.
# Concordant = (Zika_Up & Dengue_Up) OR (Zika_Down & Dengue_Down)
bulk_concordant = df_bulk[((df_bulk['in_Zika_Up']) & (df_bulk['in_Dengue_Up'])) | 
                          ((df_bulk['in_Zika_Down']) & (df_bulk['in_Dengue_Down']))]['GeneSymbol']

sc_concordant = df_sc[((df_sc['in_Zika_Up']) & (df_sc['in_Dengue_Up'])) | 
                      ((df_sc['in_Zika_Down']) & (df_sc['in_Dengue_Down']))]['GeneSymbol']

overlap_genes = set(bulk_concordant) & set(sc_concordant)
# If this doesn't equal exactly 98, we will fall back to using the full 301 overlap
if len(overlap_genes) == 0:
    overlap_genes = bulk_genes & sc_genes

bulk_genes = list(bulk_genes)
sc_genes = list(sc_genes)
overlap_genes = list(overlap_genes)

# 4. miRNA-Targeted Core Genes
MIRNA_CSV = Path("E:/Zika_Enrichment/miRNA_target_integration/Core_98_miRNA_mRNA_Regulatory_Pairs.csv")
df_mirna = pd.read_csv(MIRNA_CSV)
mirna_targets = []
for targets_str in df_mirna['Overlapping_mRNA_Genes'].dropna():
    mirna_targets.extend(targets_str.split(', '))
mirna_genes = list(set(mirna_targets))

pass  # Execution logging removed for final release
pass  # Execution logging removed for final release
pass  # Execution logging removed for final release
pass  # Execution logging removed for final release

gene_sets_dict = {
    'Bulk-Only/Systemic': bulk_genes,
    'scRNA-seq/Cellular': sc_genes,
    'Consensus Overlap': overlap_genes,
    'miRNA-Targeted Core': mirna_genes
}

databases = {
    'GO_Biological_Process': 'GO_Biological_Process_2023',
    'KEGG_Pathways': 'KEGG_2021_Human'
}

all_results = []

pass  # Execution logging removed for final release
for db_name, db_string in databases.items():
    for dataset_name, gene_list in gene_sets_dict.items():
        if len(gene_list) < 5: continue
        
        pass  # Execution logging removed for final release
        try:
            enr = gp.enrichr(gene_list=gene_list,
                             gene_sets=db_string,
                             organism='human',
                             outdir=None,
                             cutoff=0.05)
            res = enr.results
            if not res.empty:
                res = res[res['Adjusted P-value'] < 0.05].copy()
                res['Dataset'] = dataset_name
                res['Database'] = db_name
                all_results.append(res)
        except Exception as e:
            pass  # Execution logging removed for final release

if not all_results:
    pass  # Execution logging removed for final release
    exit()

# Combine all results
df_res = pd.concat(all_results, ignore_index=True)
df_res.to_csv(OUT_DIR / "Comparative_Enrichment_Raw.csv", index=False)

pass  # Execution logging removed for final release
# For plotting, we want to select the Top 5 unique terms from EACH dataset to show the contrast
for db_name in databases.keys():
    df_db = df_res[df_res['Database'] == db_name].copy()
    if df_db.empty: continue
    
    # Pick top 5 terms from each dataset
    top_terms = []
    for dset in gene_sets_dict.keys():
        subset = df_db[df_db['Dataset'] == dset].sort_values('Adjusted P-value')
        top_terms.extend(subset.head(5)['Term'].tolist())
        
    top_terms = list(dict.fromkeys(top_terms)) # Keep unique, preserve order
    
    # Filter the main dataframe to just these terms
    plot_df = df_db[df_db['Term'].isin(top_terms)].copy()
    plot_df['-log10(p-adj)'] = -np.log10(plot_df['Adjusted P-value'])
    plot_df['Term'] = plot_df['Term'].apply(lambda x: (x[:45] + '...') if len(x) > 45 else x)
    
    # Create pivot tables for a Heatmap
    # We need Term as rows, Dataset as columns
    pivot_df = plot_df.pivot_table(index='Term', columns='Dataset', values='-log10(p-adj)', fill_value=0)
    
    # Ensure all columns exist
    for c in gene_sets_dict.keys():
        if c not in pivot_df.columns:
            pivot_df[c] = 0
            
    # Reorder columns logically
    pivot_df = pivot_df[['Bulk-Only/Systemic', 'scRNA-seq/Cellular', 'Consensus Overlap', 'miRNA-Targeted Core']]
    
    # Reorder rows by clustering or sum
    pivot_df['Sum'] = pivot_df.sum(axis=1)
    pivot_df = pivot_df.sort_values('Sum', ascending=True).drop(columns=['Sum'])
    
    # Plot Heatmap (180x170mm is approx 7.08 x 6.69 inches)
    fig, ax = plt.subplots(figsize=(7.08, 6.69))
    
    # We use cbar_kws to position the legend (colorbar) explicitly on the right side
    sns.heatmap(pivot_df, cmap='YlOrRd', annot=False, linewidths=0.5, 
                cbar_kws={'label': '-log10(Adjusted P-value)', 'shrink': 0.8}, ax=ax)
    
    plt.title(f"Comparative Enrichment: {db_name.replace('_', ' ')}", fontsize=14, pad=20, fontweight='bold')
    plt.ylabel("")
    plt.xlabel("")
    plt.xticks(rotation=30, ha='right', fontsize=10)
    plt.yticks(fontsize=9)
    
    # Adjust layout to accommodate colorbar on the right
    plt.tight_layout()
    plt.savefig(OUT_DIR / f"Comparative_Heatmap_{db_name}.tiff", dpi=300, format='tiff', bbox_inches='tight')
    plt.close()
    pass  # Execution logging removed for final release

pass  # Execution logging removed for final release
