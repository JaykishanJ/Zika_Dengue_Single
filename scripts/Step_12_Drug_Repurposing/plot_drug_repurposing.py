import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import gseapy as gp
import seaborn as sns

# Nature SR Aesthetics
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

BASE_DIR = Path("E:/Zika_antigravity/SC_RNA_Zika_DENV")
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_12_Drug_Repurposing")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def main():
    pass  # Execution logging removed for final release
    deg_df = pd.read_csv(BASE_DIR / "03_results/phase4_All_DEG/All_Significant_DEGs.csv")
    core_genes = deg_df[deg_df['significance_group'] == 'Both']['symbol'].dropna().unique().tolist()
    pass  # Execution logging removed for final release
    
    pass  # Execution logging removed for final release
    try:
        # DSigDB contains drug signatures from GEO
        enr = gp.enrichr(gene_list=core_genes, gene_sets='DSigDB', outdir=None)
        res = enr.results
    except Exception as e:
        pass  # Execution logging removed for final release
        return
        
    if res.empty:
        pass  # Execution logging removed for final release
        return
        
    # Filter for significance
    res = res[res['Adjusted P-value'] < 0.05].copy()
    if res.empty:
        pass  # Execution logging removed for final release
        return
        
    # Clean the terms (DSigDB terms often look like "drugname HL60 UP" or "drugname CTD 000...")
    res['Drug_Name'] = res['Term'].apply(lambda x: x.split(' ')[0].capitalize() if ' ' in x else x.capitalize())
    
    # Deduplicate by Drug Name, keeping the most significant
    res = res.sort_values('Adjusted P-value').drop_duplicates('Drug_Name').head(15)
    
    # Save the table
    out_csv = OUT_DIR / "Top_Candidate_Drugs.csv"
    res[['Drug_Name', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Genes']].to_csv(out_csv, index=False)
    pass  # Execution logging removed for final release
    
    # Plot
    pass  # Execution logging removed for final release
    # Calculate -log10 FDR
    res['nlog10_FDR'] = -np.log10(res['Adjusted P-value'])
    res = res.sort_values('nlog10_FDR', ascending=True) # Ascending for horizontal bar chart
    
    import numpy as np
    
    fig, ax = plt.subplots(figsize=(6, 7))
    bars = ax.barh(res['Drug_Name'], res['nlog10_FDR'], color='#457B9D')
    
    # Annotate with the overlapping genes count
    for bar, overlap in zip(bars, res['Overlap']):
        count = overlap.split('/')[0]
        ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height()/2, f"n={count}", 
                va='center', ha='left', fontsize=10, color='black')
                
    ax.axvline(-np.log10(0.05), color='red', linestyle='--', linewidth=1, label='FDR 0.05')
    
    ax.set_xlabel(r'$-\log_{10}$ (Adjusted p-value)', fontsize=12, fontweight='bold')
    ax.set_title('Top Drug Repurposing Candidates\nTargeting Core Consensus Signature', fontsize=14, fontweight='bold', pad=15)
    ax.legend(frameon=False)
    
    plt.tight_layout()
    out_file = OUT_DIR / "Drug_Repurposing_Barplot.tiff"
    plt.savefig(out_file, dpi=300, format='tiff', bbox_inches='tight')
    plt.close()
    pass  # Execution logging removed for final release

if __name__ == "__main__":
    main()
