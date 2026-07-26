import sys
from pathlib import Path
_root_dir = Path(__file__).resolve().parent.parent.parent.parent
if str(_root_dir) not in sys.path:
    sys.path.append(str(_root_dir))
from src import config

import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from pathlib import Path
from scipy.stats import pearsonr
import mygene
from src.utils import pub_style

pub_style.apply_style()

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_06_Fig_6_CrossModal_Signature")
BASE_SC = config.ROOT_DIR

def main():
    pass  # Execution logging removed for final release
    SC_DEG_FILE = BASE_SC / "Step06_crossmodal_validation/crossmodal_validated_progression_signature.csv"
    if not SC_DEG_FILE.exists():
        return
        
    sc_df = pd.read_csv(SC_DEG_FILE)
    
    # ZIKV
    zikv_bulk_path = BASE_SC / "Step05b_bulk_differential_expression/bulk_GSE118305_ZIKV_macrophage_DE.csv"
    if zikv_bulk_path.exists():
        zikv_bulk = pd.read_csv(zikv_bulk_path).dropna(subset=['Gene', 'lfc', 'padj'])
        zikv_bulk = zikv_bulk[zikv_bulk['padj'] < 0.05]
        
        zikv_sc = sc_df.dropna(subset=['ZIKV_padj', 'ZIKV_log2FC', 'symbol']).copy()
        zikv_sc = zikv_sc[zikv_sc['ZIKV_padj'] < 0.05]
        zikv_sc['symbol'] = zikv_sc['symbol'].str.upper()
        
        merged_zikv = pd.merge(zikv_bulk, zikv_sc, left_on='Gene', right_on='symbol', suffixes=('_bulk', '_sc'))
        
        if len(merged_zikv) > 5:
            fig, ax = plt.subplots(figsize=(4.72, 3.93))
            sns.scatterplot(data=merged_zikv, x='lfc', y='ZIKV_log2FC', s=10, alpha=0.5, color='#457B9D', ax=ax)
            r, p = pearsonr(merged_zikv['lfc'], merged_zikv['ZIKV_log2FC'])
            
            ax.set_title("ZIKV Bulk vs scRNA-seq", fontweight='bold')
            ax.set_xlabel("Bulk log2FC")
            ax.set_ylabel("scRNA-seq log2FC")
            ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
            ax.axvline(0, color='black', linestyle='--', linewidth=0.5)
            ax.text(0.05, 0.95, f"r = {r:.2f}\np = {p:.2e}", transform=ax.transAxes, verticalalignment='top')
            
            pub_style.save_pub_fig(fig, OUT_DIR / "CrossModal_ZIKV_Correlation.tiff", width_mm=120, height_mm=100)
            plt.close()

    # DENV
    denv_bulk_path = BASE_SC / "Step05b_bulk_differential_expression/bulk_GSE279208_DENV_blood_DE.csv"
    if denv_bulk_path.exists():
        denv_bulk = pd.read_csv(denv_bulk_path).dropna(subset=['Gene', 'lfc', 'padj'])
        denv_bulk = denv_bulk[denv_bulk['padj'] < 0.05]
        
        denv_sc = sc_df.dropna(subset=['DENV_padj', 'DENV_log2FC', 'symbol']).copy()
        denv_sc = denv_sc[denv_sc['DENV_padj'] < 0.05]
        denv_sc['symbol'] = denv_sc['symbol'].str.upper()
        
        # Convert Entrez to Symbol for DENV
        mg = mygene.MyGeneInfo()
        res = mg.querymany(denv_bulk['Gene'].astype(str).tolist(), scopes='entrezgene', fields='symbol', species='human', verbose=False)
        symbol_map = {r['query']: r.get('symbol', '').upper() for r in res if r.get('symbol')}
        denv_bulk['Symbol'] = denv_bulk['Gene'].astype(str).map(symbol_map)
        
        merged_denv = pd.merge(denv_bulk, denv_sc, left_on='Symbol', right_on='symbol', suffixes=('_bulk', '_sc'))
        
        if len(merged_denv) > 5:
            fig, ax = plt.subplots(figsize=(4.72, 3.93))
            sns.scatterplot(data=merged_denv, x='lfc', y='DENV_log2FC', s=10, alpha=0.5, color='#E63946', ax=ax)
            r, p = pearsonr(merged_denv['lfc'], merged_denv['DENV_log2FC'])
            
            ax.set_title("DENV Bulk vs scRNA-seq", fontweight='bold')
            ax.set_xlabel("Bulk log2FC")
            ax.set_ylabel("scRNA-seq log2FC")
            ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
            ax.axvline(0, color='black', linestyle='--', linewidth=0.5)
            ax.text(0.05, 0.95, f"r = {r:.2f}\np = {p:.2e}", transform=ax.transAxes, verticalalignment='top')
            
            pub_style.save_pub_fig(fig, OUT_DIR / "CrossModal_DENV_Correlation.tiff", width_mm=120, height_mm=100)
            plt.close()

if __name__ == "__main__":
    main()
