import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import mygene

# Nature SR Aesthetics
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# Paths
DATA_PATH = Path("E:/Zika_Enrichment/cell_state_analysis/adata_annotated.h5ad")
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_03_Fig_3_scRNA_CellStates")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def add_umap_arrows(ax):
    """Add Nature SR style UMAP corner arrows"""
    x0, y0 = 0.05, 0.05
    length = 0.15
    ax.annotate("", xy=(x0 + length, y0), xytext=(x0, y0),
                xycoords='axes fraction', arrowprops=dict(arrowstyle="->", lw=1.5, color='black'))
    ax.annotate("", xy=(x0, y0 + length), xytext=(x0, y0),
                xycoords='axes fraction', arrowprops=dict(arrowstyle="->", lw=1.5, color='black'))
    ax.text(x0 + length/2, y0 - 0.02, "UMAP_1", ha='center', va='top', transform=ax.transAxes, fontsize=10, fontweight='bold')
    ax.text(x0 - 0.02, y0 + length/2, "UMAP_2", ha='right', va='center', rotation=90, transform=ax.transAxes, fontsize=10, fontweight='bold')

def main():
    pass  # Execution logging removed for final release
    adata = sc.read_h5ad(DATA_PATH)
    
    if 'leiden' not in adata.obs:
        pass  # Execution logging removed for final release
        if 'neighbors' not in adata.uns:
            sc.pp.neighbors(adata)
        sc.tl.leiden(adata, resolution=0.5, key_added='leiden_unannotated', flavor='igraph', n_iterations=2, directed=False)
    else:
        adata.obs['leiden_unannotated'] = adata.obs['leiden']
        
    pass  # Execution logging removed for final release
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color='leiden_unannotated', ax=ax, show=False, 
               legend_loc='on data', frameon=False, size=20, alpha=0.8, title="UMAP (Unannotated Clusters)")
    add_umap_arrows(ax)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "UMAP_Unannotated_Clusters.tiff", dpi=300, bbox_inches='tight', format='tiff')
    plt.close()

    if 'infection_state' in adata.obs:
        pass  # Execution logging removed for final release
        fig, ax = plt.subplots(figsize=(7, 6))
        sc.pl.umap(adata, color='infection_state', ax=ax, show=False, 
                   legend_loc='on data', legend_fontsize=10, legend_fontoutline=2,
                   frameon=False, size=20, alpha=0.8, title="Annotated Cell States")
        add_umap_arrows(ax)
        plt.tight_layout()
        plt.savefig(OUT_DIR / "UMAP_Annotated_States.tiff", dpi=300, bbox_inches='tight', format='tiff')
        plt.close()
        
        pass  # Execution logging removed for final release
        sc.tl.rank_genes_groups(adata, 'infection_state', method='wilcoxon', key_added='rank_genes_infection')
        
        result = adata.uns['rank_genes_infection']
        groups = result['names'].dtype.names
        marker_dict = {}
        top_plot_genes = []
        dotplot_genes = []
        
        for group in groups:
            top_genes = result['names'][group][:5].tolist()
            marker_dict[group] = top_genes
            
            # For feature plots
            for g in top_genes:
                if g not in top_plot_genes and not g.startswith('__'):
                    top_plot_genes.append(g)
                    break
                    
            # For dotplot
            for g in result['names'][group][:3]:
                if g not in dotplot_genes and not g.startswith('__'):
                    dotplot_genes.append(g)
        
        # Map Ensembl IDs to Symbols using mygene
        all_to_map = list(set([g for genes in marker_dict.values() for g in genes] + top_plot_genes + dotplot_genes))
        pass  # Execution logging removed for final release
        mg = mygene.MyGeneInfo()
        res = mg.querymany(all_to_map, scopes='ensembl.gene', fields='symbol', species='human', verbose=False)
        symbol_map = {r['query']: r.get('symbol', r['query']) for r in res}
        
        # Make a copy of var_names, replace mapped ones, ensure uniqueness
        new_var_names = []
        seen = set()
        for v in adata.var_names:
            sym = symbol_map.get(v, v)
            if sym in seen:
                sym = f"{sym}_{v}" # Append Ensembl if duplicate
            seen.add(sym)
            new_var_names.append(sym)
        
        adata.var_names = new_var_names
        
        # Map our tracking lists to the new symbols
        mapped_marker_dict = {k: [symbol_map.get(g, g) for g in v] for k, v in marker_dict.items()}
        mapped_top_plot_genes = [symbol_map.get(g, g) for g in top_plot_genes]
        mapped_dotplot_genes = [symbol_map.get(g, g) for g in dotplot_genes]
        
        marker_df = pd.DataFrame({
            "Cell State": list(mapped_marker_dict.keys()),
            "Top Marker Genes": [", ".join(genes) for genes in mapped_marker_dict.values()]
        })
        marker_csv_path = OUT_DIR / "Marker_Genes_Table.csv"
        marker_df.to_csv(marker_csv_path, index=False)
        pass  # Execution logging removed for final release

        pass  # Execution logging removed for final release
        genes_to_plot = [g for g in mapped_top_plot_genes if g in adata.var_names][:4]
        if genes_to_plot:
            fig, axes = plt.subplots(1, len(genes_to_plot), figsize=(4 * len(genes_to_plot), 4))
            if len(genes_to_plot) == 1:
                axes = [axes]
                
            for i, gene in enumerate(genes_to_plot):
                sc.pl.umap(adata, color=gene, ax=axes[i], show=False, 
                           color_map='Purples', frameon=False, size=15, title=gene)
                add_umap_arrows(axes[i])
                
            plt.tight_layout()
            plt.savefig(OUT_DIR / "Validation_Feature_Plots.tiff", dpi=300, bbox_inches='tight', format='tiff')
            plt.close()

        pass  # Execution logging removed for final release
        # Recompute rank_genes_groups so it uses the new symbols internally if we wanted to show all,
        # but we are just passing var_names to dotplot which works fine.
        dp = sc.pl.dotplot(adata, var_names=mapped_dotplot_genes, groupby='infection_state', 
                           standard_scale='var', cmap='Purples', show=False, return_fig=True)
        dp.savefig(OUT_DIR / "DotPlot_Example.tiff", dpi=300, bbox_inches='tight', format='tiff')
        pass  # Execution logging removed for final release

if __name__ == "__main__":
    main()
