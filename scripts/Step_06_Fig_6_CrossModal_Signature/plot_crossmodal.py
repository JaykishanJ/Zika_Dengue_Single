import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path

sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

from venn import venn
import mygene
from pathlib import Path
import os

BASE_SC = Path("E:/Zika_antigravity/SC_RNA_Zika_DENV")
SC_DEG_FILE = BASE_SC / "03_results/phase4_All_DEG/All_Significant_DEGs.csv"
BULK_DIR = BASE_SC / "Bulk_RNA_seq"

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_06_Fig_6_CrossModal_Signature")

bulk_datasets = {
    "ZIKV": [
        BULK_DIR / "GSE118305/Results/Tables/ZIKV_positive_vs_Mock_All_Genes.csv",
        BULK_DIR / "GSE78711/Results/Tables/GSE78711_Processed_DEG_ZIKV_vs_Mock_All_Genes.csv"
    ],
    "DENV": [
        BULK_DIR / "GSE279208/Results/Tables/Dengue_vs_healthy_All_Genes.csv"
    ]
}

def get_bulk_genes(disease, direction):
    files = bulk_datasets[disease]
    genes = set()
    for f in files:
        if not f.exists():
            continue
        df = pd.read_csv(f)
        sig_df = df.dropna(subset=['padj', 'log2FoldChange', 'Gene']).copy()
        sig_df = sig_df[sig_df['padj'] < 0.05]
        if direction == 'Up':
            sig_df = sig_df[sig_df['log2FoldChange'] > 0]
        else:
            sig_df = sig_df[sig_df['log2FoldChange'] < 0]
        
        if disease == 'DENV':
            entrez_ids = sig_df['Gene'].astype(str).tolist()
            mg = mygene.MyGeneInfo()
            res = mg.querymany(entrez_ids, scopes='entrezgene', fields='symbol', species='human', verbose=False)
            symbols = [r.get('symbol', '').upper() for r in res if r.get('symbol')]
            genes.update(symbols)
        else:
            sig_genes = sig_df['Gene'].astype(str).str.upper().tolist()
            genes.update(sig_genes)
    return genes

def get_sc_genes(disease, direction):
    sc_df = pd.read_csv(SC_DEG_FILE)
    col = f"{disease}_regulation"
    if col not in sc_df.columns:
        return set()
    genes = sc_df[sc_df[col] == direction]['symbol'].astype(str).str.upper().tolist()
    return set(genes)

def main():
    pass  # Execution logging removed for final release
    bulk_zika_up = get_bulk_genes('ZIKV', 'Up')
    bulk_zika_down = get_bulk_genes('ZIKV', 'Down')
    bulk_denv_up = get_bulk_genes('DENV', 'Up')
    bulk_denv_down = get_bulk_genes('DENV', 'Down')

    pass  # Execution logging removed for final release
    sc_zika_up = get_sc_genes('ZIKV', 'Up')
    sc_zika_down = get_sc_genes('ZIKV', 'Down')
    sc_denv_up = get_sc_genes('DENV', 'Up')
    sc_denv_down = get_sc_genes('DENV', 'Down')

    # Plot Bulk Venn
    bulk_labels = {
        "Upregulated\nDEGs in Zika": bulk_zika_up,
        "Downregulated\nDEGs in Zika": bulk_zika_down,
        "Upregulated\nDEGs in Dengue": bulk_denv_up,
        "Downregulated\nDEGs in Dengue": bulk_denv_down
    }
    
    fig, ax = plt.subplots(figsize=(10, 10))
    venn(bulk_labels, ax=ax)
    plt.title("Overlapping DEGs in Bulk RNA-seq (Zika vs Dengue)", fontsize=16)
    plt.savefig(OUT_DIR / "4way_Venn_Bulk_Only.png", dpi=300, bbox_inches='tight')
    plt.close()
    pass  # Execution logging removed for final release

    # Plot SC Venn
    sc_labels = {
        "Upregulated\nDEGs in Zika": sc_zika_up,
        "Downregulated\nDEGs in Zika": sc_zika_down,
        "Upregulated\nDEGs in Dengue": sc_denv_up,
        "Downregulated\nDEGs in Dengue": sc_denv_down
    }
    
    fig, ax = plt.subplots(figsize=(10, 10))
    venn(sc_labels, ax=ax)
    plt.title("Overlapping DEGs in Single-Cell RNA-seq (Zika vs Dengue)", fontsize=16)
    plt.savefig(OUT_DIR / "4way_Venn_SC_Only.png", dpi=300, bbox_inches='tight')
    plt.close()
    pass  # Execution logging removed for final release

    # Plot Union Venn
    union_labels = {
        "Upregulated\nDEGs in Zika": bulk_zika_up | sc_zika_up,
        "Downregulated\nDEGs in Zika": bulk_zika_down | sc_zika_down,
        "Upregulated\nDEGs in Dengue": bulk_denv_up | sc_denv_up,
        "Downregulated\nDEGs in Dengue": bulk_denv_down | sc_denv_down
    }
    
    # 180x140mm -> 7.08x5.51 inches
    fig, ax = plt.subplots(figsize=(7.08, 5.51))
    # We will use Pastel1
    venn(union_labels, ax=ax, cmap="Set3")
    
    # Custom Legend
    if ax.get_legend():
        ax.get_legend().remove()
    
    cmap = plt.get_cmap("Set3")
    handles = []
    from matplotlib.patches import Patch
    for i, label_name in enumerate(union_labels.keys()):
        patch = Patch(color=cmap(i), label=label_name.replace('\n', ' '))
        handles.append(patch)
        
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.05, 0.5), ncol=1, frameon=False, fontsize=10)
    
    plt.title("Cross-Modal Gene Signature Overlap", fontsize=14, fontweight='bold', pad=20)
    pub_style.save_pub_fig(plt.gcf(), OUT_DIR / "CrossModal_Signature.tiff", width_mm=180, height_mm=120, tight_layout=False)
    plt.close()
    pass  # Execution logging removed for final release

if __name__ == "__main__":
    main()
