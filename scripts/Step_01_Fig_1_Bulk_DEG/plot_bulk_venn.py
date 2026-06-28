import pandas as pd
import matplotlib.pyplot as plt

# Nature SR Aesthetics
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

from matplotlib.patches import Patch
from venn import venn
import mygene
from pathlib import Path

BASE_SC = Path("E:/Zika_antigravity/SC_RNA_Zika_DENV")
SC_DEG_FILE = BASE_SC / "03_results/phase4_All_DEG/All_Significant_DEGs.csv"
BULK_DIR = BASE_SC / "Bulk_RNA_seq"
OUT_DIR = Path("E:/Zika_Enrichment/4way_venn_analysis")

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

def plot_custom_venn(labels_dict, title, out_path):
    # 180mm x 160mm is approximately 7.08 x 6.30 inches
    fig, ax = plt.subplots(figsize=(7.08, 6.30))
    
    # We use Pastel1 or Set3 which gives nice pastel colors similar to the user's example
    cmap_name = "Set3"
    
    # Generate Venn diagram. legend_loc=None tries to suppress the built-in legend.
    venn(labels_dict, ax=ax, cmap=cmap_name)
    
    # Hide default legend if venn package creates one
    if ax.get_legend():
        ax.get_legend().remove()
        
    cmap = plt.get_cmap(cmap_name)
    handles = []
    # In the venn library, the colors are assigned to sets in the order they appear
    for i, label_name in enumerate(labels_dict.keys()):
        patch = Patch(color=cmap(i), label=label_name)
        handles.append(patch)
        
    # Put legend on the RIGHT SIDE for volcano/venn
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.05, 0.5), ncol=1, frameon=False, fontsize=10)
    
    plt.title(title, fontsize=14)
    # Adjust layout to accommodate legend on the right
    plt.tight_layout()
    plt.subplots_adjust(right=0.75) 
    plt.savefig(out_path, dpi=300, bbox_inches='tight', format='tiff')
    plt.close()

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

    bulk_labels = {
        "Zika Upregulated": bulk_zika_up,
        "Zika Downregulated": bulk_zika_down,
        "Dengue Upregulated": bulk_denv_up,
        "Dengue Downregulated": bulk_denv_down
    }
    
    sc_labels = {
        "Zika Upregulated": sc_zika_up,
        "Zika Downregulated": sc_zika_down,
        "Dengue Upregulated": sc_denv_up,
        "Dengue Downregulated": sc_denv_down
    }
    
    out_dir = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_01_Fig_1_Bulk_DEG")
    plot_custom_venn(bulk_labels, "Overlapping DEGs in Bulk RNA-seq (Zika vs Dengue)", out_dir / "Final_Venn_Bulk.tiff")
    pass  # Execution logging removed for final release
    
    # We also do the SC venn here or let the Step 6 handle it. We will save it here too.
    plot_custom_venn(sc_labels, "Overlapping DEGs in Single-Cell RNA-seq (Zika vs Dengue)", out_dir / "Final_Venn_SC.tiff")
    pass  # Execution logging removed for final release

if __name__ == "__main__":
    main()
