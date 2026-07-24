import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import sys
from pathlib import Path

sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_07_Fig_7_miRNA_Targets")
DATA_FILE = OUT_DIR / "Core_98_miRNA_mRNA_Regulatory_Pairs.csv"

def main():
    if not DATA_FILE.exists():
        pass  # Execution logging removed for final release
        return
        
    df = pd.read_csv(DATA_FILE)
    df['Target_Count'] = df['Overlapping_mRNA_Genes'].apply(lambda x: len(str(x).split(',')) if pd.notna(x) else 0)
    top_mirnas = df.sort_values('Target_Count', ascending=False).head(3)
    
    G = nx.DiGraph()
    
    # Layer 1: miRNAs
    # Layer 2: mRNA Targets (top 15 unique targets across these 3)
    # Layer 3: Pathways (simplified)
    
    pathways = {
        'Antiviral Response': ['OAS1', 'OAS2', 'IFIT1', 'ISG15', 'MX1', 'RSAD2'],
        'Inflammation': ['CXCL10', 'CCL5', 'IL6', 'TNF', 'CXCL8'],
        'Apoptosis': ['CASP3', 'CASP8', 'FAS', 'BAX']
    }
    
    all_targets = set()
    for _, row in top_mirnas.iterrows():
        mirna = row['miRNA']
        G.add_node(mirna, layer=0, type='miRNA')
        
        targets = [t.strip() for t in str(row['Overlapping_mRNA_Genes']).split(',') if t.strip()]
        for t in targets:
            all_targets.add(t)
            G.add_node(t, layer=1, type='mRNA')
            G.add_edge(mirna, t)
            
    # Add pathway edges
    for t in all_targets:
        for p_name, p_genes in pathways.items():
            if t in p_genes:
                G.add_node(p_name, layer=2, type='Pathway')
                G.add_edge(t, p_name)
                
    # Remove mRNAs that don't connect to pathways just to keep the plot clean
    nodes_to_remove = [n for n, attr in G.nodes(data=True) if attr['type'] == 'mRNA' and G.out_degree(n) == 0]
    G.remove_nodes_from(nodes_to_remove)
    
    fig, ax = plt.subplots(figsize=(7.08, 5))
    
    pos = nx.multipartite_layout(G, subset_key="layer", align='horizontal')
    
    # Rotate layout 90 degrees
    pos = {node: (y, -x) for node, (x, y) in pos.items()}
    
    # Draw nodes
    colors = []
    for n, attr in G.nodes(data=True):
        if attr['type'] == 'miRNA': colors.append('#E63946')
        elif attr['type'] == 'mRNA': colors.append('#457B9D')
        else: colors.append('#2A9D8F')
        
    nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=500, ax=ax)
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5, edge_color='gray', arrows=True, arrowsize=10, ax=ax)
    
    # Custom labels to avoid overlap
    labels = {n: n for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight='bold')
    
    ax.set_title("Regulatory Cascade: miRNA \u2192 mRNA \u2192 Pathway", fontweight='bold')
    
    pub_style.save_pub_fig(fig, OUT_DIR / "miRNA_Target_Cascade_Network.tiff", width_mm=180, height_mm=120)
    plt.close()

if __name__ == "__main__":
    main()
