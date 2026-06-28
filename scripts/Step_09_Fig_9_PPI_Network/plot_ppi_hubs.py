import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from pathlib import Path

sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_09_Fig_9_PPI_Network")
DATA_FILE = Path("E:/Zika_Enrichment/Publication_Standard_PPI_Genes.csv")

def main():
    if not DATA_FILE.exists():
        pass  # Execution logging removed for final release
        return
        
    df = pd.read_csv(DATA_FILE)
    query_genes = df['GeneSymbol'].tolist()
    
    pass  # Execution logging removed for final release
    string_api_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "\r".join(query_genes),
        "species": 9606,  
        "required_score": 900,  
        "add_nodes": 10  
    }
    
    response = requests.post(string_api_url, data=params)
    if response.status_code != 200:
        pass  # Execution logging removed for final release
        return
        
    interactions = response.json()
    
    pass  # Execution logging removed for final release
    G = nx.Graph()
    for inter in interactions:
        nodeA = inter['preferredName_A']
        nodeB = inter['preferredName_B']
        G.add_edge(nodeA, nodeB)
        
    G.remove_nodes_from(list(nx.isolates(G)))
    
    degrees = dict(G.degree())
    df_deg = pd.DataFrame(list(degrees.items()), columns=['Gene', 'Degree'])
    df_sorted = df_deg.sort_values('Degree', ascending=False).head(15)
    
    pass  # Execution logging removed for final release
    fig, ax = plt.subplots(figsize=(4.72, 3.93))
    
    sns.barplot(data=df_sorted, x='Degree', y='Gene', palette='Reds_r', edgecolor='black', linewidth=0.5, ax=ax)
    
    ax.set_title("Top Protein-Protein Interaction Hubs", fontweight='bold', pad=10)
    ax.set_xlabel("Node Degree (STRING Confidence > 0.9)")
    ax.set_ylabel("")
    
    pub_style.save_pub_fig(fig, OUT_DIR / "PPI_Network_Hubs.tiff", width_mm=120, height_mm=100)
    plt.close()

if __name__ == "__main__":
    main()
