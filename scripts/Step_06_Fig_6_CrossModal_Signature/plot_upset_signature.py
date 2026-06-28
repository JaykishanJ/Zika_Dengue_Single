import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path
from upsetplot import from_contents, plot

sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_06_Fig_6_CrossModal_Signature")
BASE_SC = Path("E:/Zika_antigravity/SC_RNA_Zika_DENV")
SC_DEG_FILE = BASE_SC / "03_results/phase4_All_DEG/All_Significant_DEGs.csv"
BULK_DIR = BASE_SC / "Bulk_RNA_seq"

bulk_datasets = {
    "ZIKV": [
        BULK_DIR / "GSE118305/Results/Tables/ZIKV_positive_vs_Mock_All_Genes.csv",
        BULK_DIR / "GSE78711/Results/Tables/ZIKV_positive_vs_Mock_All_Genes.csv" 
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
            import mygene
            mg = mygene.MyGeneInfo()
            res = mg.querymany(entrez_ids, scopes='entrezgene', fields='symbol', species='human', verbose=False)
            symbols = [r.get('symbol', '').upper() for r in res if r.get('symbol')]
            genes.update(symbols)
        else:
            sig_genes = sig_df['Gene'].astype(str).str.upper().tolist()
            genes.update(sig_genes)
    return genes

def get_sc_genes(disease, direction):
    if not SC_DEG_FILE.exists():
        return set()
    sc_df = pd.read_csv(SC_DEG_FILE)
    
    if disease == 'ZIKV':
        df = sc_df.dropna(subset=['ZIKV_padj', 'ZIKV_log2FC', 'symbol']).copy()
        df = df[df['ZIKV_padj'] < 0.05]
        if direction == 'Up':
            df = df[df['ZIKV_log2FC'] > 0]
        else:
            df = df[df['ZIKV_log2FC'] < 0]
    else:
        df = sc_df.dropna(subset=['DENV_padj', 'DENV_log2FC', 'symbol']).copy()
        df = df[df['DENV_padj'] < 0.05]
        if direction == 'Up':
            df = df[df['DENV_log2FC'] > 0]
        else:
            df = df[df['DENV_log2FC'] < 0]
            
    return set(df['symbol'].str.upper().tolist())

def main():
    pass  # Execution logging removed for final release
    zika_bulk_up = get_bulk_genes("ZIKV", "Up")
    zika_bulk_down = get_bulk_genes("ZIKV", "Down")
    denv_bulk_up = get_bulk_genes("DENV", "Up")
    denv_bulk_down = get_bulk_genes("DENV", "Down")
    
    zika_sc_up = get_sc_genes("ZIKV", "Up")
    zika_sc_down = get_sc_genes("ZIKV", "Down")
    denv_sc_up = get_sc_genes("DENV", "Up")
    denv_sc_down = get_sc_genes("DENV", "Down")
    
    zika_up = zika_bulk_up | zika_sc_up
    zika_down = zika_bulk_down | zika_sc_down
    denv_up = denv_bulk_up | denv_sc_up
    denv_down = denv_bulk_down | denv_sc_down
    
    contents = {
        'ZIKV Up': zika_up,
        'ZIKV Down': zika_down,
        'DENV Up': denv_up,
        'DENV Down': denv_down
    }
    
    pass  # Execution logging removed for final release
    data = from_contents(contents)
    
    fig = plt.figure(figsize=(7.08, 5))
    plot(data, fig=fig, sort_by='cardinality', show_counts=True, element_size=40)
    plt.suptitle("Cross-Modal Gene Signature Overlap", fontweight='bold')
    
    pub_style.save_pub_fig(fig, OUT_DIR / "UpSet_CrossModal.tiff", width_mm=180, height_mm=120)
    plt.close()

if __name__ == "__main__":
    main()
