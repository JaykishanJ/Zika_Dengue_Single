import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import scipy.sparse

def export_data():
    adata_path = Path("E:/Zika_Enrichment/cell_state_analysis/adata_trajectory.h5ad")
    out_dir = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_04_Fig_4_Pseudotime")
    
    print(f"Loading {adata_path}...")
    adata = sc.read_h5ad(adata_path)
    
    print("Reconstructing counts from log_norm layer...")
    log_norm = adata.layers['log_norm']
    
    # Exponentiate and round to get approximate integer counts
    if scipy.sparse.issparse(log_norm):
        counts = log_norm.copy()
        counts.data = np.round(np.expm1(counts.data)).astype(int)
        counts = counts.T
    else:
        counts = np.round(np.expm1(log_norm)).astype(int).T
        
    genes = adata.var_names
    cells = adata.obs_names
    
    print("Exporting counts...")
    df_counts = pd.DataFrame(counts.toarray() if scipy.sparse.issparse(counts) else counts, 
                             index=genes, columns=cells)
    df_counts.to_csv(out_dir / "tradeseq_counts.csv")
    
    print("Exporting metadata...")
    adata.obs.to_csv(out_dir / "tradeseq_metadata.csv")
    print("Done.")

if __name__ == "__main__":
    export_data()
