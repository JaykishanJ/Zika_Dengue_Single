import sys
from pathlib import Path
_root_dir = Path(__file__).resolve().parent.parent.parent.parent
if str(_root_dir) not in sys.path:
    sys.path.append(str(_root_dir))
from src import config

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path
from upsetplot import from_contents, plot
from src.utils import pub_style

# We intentionally do NOT apply_style() here because upsetplot 
# internally manages grid subplots and crashes when custom rcParams (like spines=False) are used.

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_06_Fig_6_CrossModal_Signature")
BASE_SC = config.ROOT_DIR
SC_DEG_FILE = BASE_SC / "Step06_crossmodal_validation/crossmodal_validated_progression_signature.csv"
BULK_DIR = BASE_SC / "Step05b_bulk_differential_expression"

bulk_datasets = {
    "ZIKV": [
        BULK_DIR / "bulk_GSE118305_ZIKV_macrophage_DE.csv",
        BULK_DIR / "GSE78711/Results/Tables/ZIKV_positive_vs_Mock_All_Genes.csv" 
    ],
    "DENV": [
        BULK_DIR / "bulk_GSE279208_DENV_blood_DE.csv"
    ]
}

def get_bulk_genes(disease, direction):
    files = bulk_datasets[disease]
    genes = set()
    for f in files:
        if not f.exists():
            continue
        df = pd.read_csv(f)
        sig_df = df.dropna(subset=['padj', 'lfc', 'Gene']).copy()
        sig_df = sig_df[sig_df['padj'] < 0.05]
        if direction == 'Up':
            sig_df = sig_df[sig_df['lfc'] > 0]
        else:
            sig_df = sig_df[sig_df['lfc'] < 0]
        
        if disease == 'DENV':
            entrez_ids = sig_df['SYM'].astype(str).tolist()
            import mygene
            mg = mygene.MyGeneInfo()
            res = mg.querymany(entrez_ids, scopes='entrezgene', fields='symbol', species='human', verbose=False)
            symbols = [r.get('symbol', '').upper() for r in res if r.get('symbol')]
            genes.update(symbols)
        else:
            sig_genes = sig_df['SYM'].astype(str).str.upper().tolist()
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
    plot(data, fig=fig, sort_by='cardinality', element_size=40)
    plt.suptitle("Cross-Modal Gene Signature Overlap", fontweight='bold')
    
    plt.savefig(OUT_DIR / "UpSet_CrossModal.tiff", dpi=300, format='tiff')
    plt.close()

if __name__ == "__main__":
    main()
