import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from pathlib import Path

# Import our new publication style module
sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_02_Fig_2_Functional_Enrichment")

def get_core_signature():
    BULK_CSV = Path("E:/Zika_Enrichment/4way_venn_analysis/Final_Bulk_Overlaps.csv")
    SC_CSV = Path("E:/Zika_Enrichment/4way_venn_analysis/Final_SC_Overlaps.csv")
    
    if not BULK_CSV.exists() or not SC_CSV.exists():
        return []
        
    df_bulk = pd.read_csv(BULK_CSV)
    df_sc = pd.read_csv(SC_CSV)
    
    bulk_genes = set(df_bulk['GeneSymbol'].dropna().unique())
    sc_genes = set(df_sc['GeneSymbol'].dropna().unique())
    
    bulk_concordant = df_bulk[((df_bulk['in_Zika_Up']) & (df_bulk['in_Dengue_Up'])) | 
                              ((df_bulk['in_Zika_Down']) & (df_bulk['in_Dengue_Down']))]['GeneSymbol']
                              
    sc_concordant = df_sc[((df_sc['in_Zika_Up']) & (df_sc['in_Dengue_Up'])) | 
                          ((df_sc['in_Zika_Down']) & (df_sc['in_Dengue_Down']))]['GeneSymbol']
                          
    overlap_genes = set(bulk_concordant) & set(sc_concordant)
    if len(overlap_genes) == 0:
        overlap_genes = bulk_genes & sc_genes
        
    return list(overlap_genes)

def plot_dotplot(enr_res, title, out_path, top_n=15):
    res_df = enr_res.res2d.copy()
    if len(res_df) == 0:
        return
        
    # Clean up long term names
    res_df['Term'] = res_df['Term'].apply(lambda x: x.split('(GO:')[0].strip() if '(GO:' in x else x)
    res_df['Term'] = res_df['Term'].apply(lambda x: (x[:40] + '...') if len(x) > 43 else x)
    
    # Sort by adjusted p-value
    res_df = res_df.sort_values('Adjusted P-value').head(top_n)
    
    # Calculate GeneRatio (hits / signature size)
    res_df['Gene_Count'] = res_df['Overlap'].str.split('/').apply(lambda x: int(x[0]))
    total_genes = int(res_df['Overlap'].iloc[0].split('/')[1])
    res_df['GeneRatio'] = res_df['Gene_Count'] / total_genes
    
    # -log10 FDR for color
    res_df['nlog10_FDR'] = -np.log10(res_df['Adjusted P-value'].clip(lower=1e-50))
    
    # Reverse to plot most significant at the top
    res_df = res_df.iloc[::-1]
    
    fig, ax = plt.subplots(figsize=(4.72, 3.93)) # 120x100mm
    
    scatter = ax.scatter(x=res_df['GeneRatio'], y=res_df['Term'], 
                         s=res_df['Gene_Count']*15, c=res_df['nlog10_FDR'], 
                         cmap='viridis', alpha=0.8, edgecolors='black', linewidth=0.5)
                         
    ax.set_title(title, pad=15, fontweight='bold')
    ax.set_xlabel('Gene Ratio')
    ax.set_ylabel('')
    
    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label(r'$-\log_{10}$ FDR', rotation=270, labelpad=15)
    
    # Make room for labels
    plt.subplots_adjust(left=0.45, right=0.95)
    
    pub_style.save_pub_fig(fig, out_path, width_mm=120, height_mm=100)
    plt.close()

def main():
    core_genes = get_core_signature()
    if not core_genes:
        pass  # Execution logging removed for final release
        return
        
    pass  # Execution logging removed for final release
    
    # GO BP
    pass  # Execution logging removed for final release
    go_enr = gp.enrichr(gene_list=core_genes, gene_sets='GO_Biological_Process_2021', organism='human', cutoff=0.05)
    plot_dotplot(go_enr, "GO Biological Process Enrichment", OUT_DIR / "GO_BP_DotPlot.tiff")
    
    # KEGG
    pass  # Execution logging removed for final release
    kegg_enr = gp.enrichr(gene_list=core_genes, gene_sets='KEGG_2021_Human', organism='human', cutoff=0.05)
    plot_dotplot(kegg_enr, "KEGG Pathway Enrichment", OUT_DIR / "KEGG_DotPlot.tiff")

if __name__ == "__main__":
    main()
