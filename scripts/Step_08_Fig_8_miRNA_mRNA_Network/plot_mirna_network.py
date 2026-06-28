import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Nature SR Aesthetics
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

from pathlib import Path
import math

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_08_Fig_8_miRNA_mRNA_Network")
OUT_DIR.mkdir(parents=True, exist_ok=True)

INPUT_CSV = Path("E:/Zika_Enrichment/miRNA_target_integration/Core_98_miRNA_mRNA_Regulatory_Pairs.csv")

pass  # Execution logging removed for final release
df = pd.read_csv(INPUT_CSV)

# Create an undirected graph
G = nx.Graph()

# Keep track of node attributes
node_types = {}  # 'miRNA' or 'mRNA'
mirna_virus = {} # 'DENV', 'ZIKV', or 'SHARED'

pass  # Execution logging removed for final release
for idx, row in df.iterrows():
    mirna = row['miRNA']
    virus = row['Virus']
    targets_str = str(row['Overlapping_mRNA_Genes'])
    if pd.isna(targets_str) or targets_str.strip() == "":
        continue
        
    # Handle shared miRNAs (some miRNAs might exist in both sets)
    if mirna in mirna_virus:
        if mirna_virus[mirna] != virus:
            mirna_virus[mirna] = 'SHARED'
    else:
        mirna_virus[mirna] = virus
        
    G.add_node(mirna)
    node_types[mirna] = 'miRNA'
    
    targets = targets_str.split(', ')
    for t in targets:
        t = t.strip()
        if not t: continue
        G.add_node(t)
        node_types[t] = 'mRNA'
        G.add_edge(mirna, t)

pass  # Execution logging removed for final release

# Calculate Degree for sizing
degrees = dict(G.degree())

# Assign colors and sizes
node_colors = []
node_sizes = []
labels = {}

for node in G.nodes():
    deg = degrees[node]
    
    # Base size formula (scale up)
    size = (deg * 150) + 300
    
    if node_types[node] == 'miRNA':
        v = mirna_virus.get(node, 'DENV')
        if v == 'DENV':
            node_colors.append('#E63946') # Red
        elif v == 'ZIKV':
            node_colors.append('#457B9D') # Blue
        else:
            node_colors.append('#9B5DE5') # Purple (Shared)
        
        # Only label highly connected miRNAs to avoid clutter
        if deg >= 3:
            labels[node] = node
        else:
            labels[node] = ""
            
    else: # mRNA
        node_colors.append('#F4A261') # Gold/Orange
        node_sizes.append((deg * 200) + 600) # mRNAs slightly larger base
        labels[node] = node # Label all mRNA targets
        continue
        
    node_sizes.append(size)

pass  # Execution logging removed for final release
# Kamada-Kawai or Spring layout usually works best for bipartite-like networks
pos = nx.spring_layout(G, k=0.3, iterations=100, seed=42)

pass  # Execution logging removed for final release
# 180x170mm -> 7.08x6.69 inches
fig, ax = plt.subplots(figsize=(7.08, 6.69))

# Draw edges
nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='#B0BEC5', width=1.0, ax=ax)

# Draw nodes
nx.draw_networkx_nodes(G, pos, 
                       node_color=node_colors, 
                       node_size=node_sizes, 
                       alpha=0.9, 
                       edgecolors='white', 
                       linewidths=0.5,
                       ax=ax)

# Draw labels with background box for readability
for node, (x, y) in pos.items():
    if labels.get(node, ""):
        font_size = 5 if node_types[node] == 'miRNA' else 7
        font_weight = 'normal' if node_types[node] == 'miRNA' else 'bold'
        ax.text(x, y, s=labels[node], fontsize=font_size, fontweight=font_weight,
                 horizontalalignment='center', verticalalignment='center',
                 bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', boxstyle='round,pad=0.1'))

# Create custom legend
import matplotlib.lines as mlines
dengue_patch = mlines.Line2D([], [], color='#E63946', marker='o', linestyle='None', markersize=8, label='Dengue miRNA')
zika_patch = mlines.Line2D([], [], color='#457B9D', marker='o', linestyle='None', markersize=8, label='Zika miRNA')
shared_patch = mlines.Line2D([], [], color='#9B5DE5', marker='o', linestyle='None', markersize=8, label='Shared miRNA')
mrna_patch = mlines.Line2D([], [], color='#F4A261', marker='o', linestyle='None', markersize=8, label='Target mRNA Gene')

ax.legend(handles=[dengue_patch, zika_patch, shared_patch, mrna_patch], 
           loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=8, frameon=False, title='Node Types', title_fontsize=9)

plt.title("miRNA-mRNA Master Regulatory Network", fontsize=12, pad=10, fontweight='bold')
plt.axis('off')

# Save
out_file = OUT_DIR / "miRNA_mRNA_Regulatory_Network.tiff"
plt.tight_layout()
plt.subplots_adjust(right=0.75) 
plt.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white', format='tiff')
plt.close()

pass  # Execution logging removed for final release

pass  # Execution logging removed for final release
import numpy as np
from collections import Counter

degree_sequence = sorted((d for n, d in G.degree()), reverse=True)
degree_counts = Counter(degree_sequence)
degrees = list(degree_counts.keys())
counts = list(degree_counts.values())

fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(degrees, counts, color='#E63946', alpha=0.8)

# Calculate a linear fit in log-log space (proxy for power-law)
log_deg = np.log10(degrees)
log_cnt = np.log10(counts)
slope, intercept = np.polyfit(log_deg, log_cnt, 1)

# Plot the fit
fit_counts = 10**(intercept + slope * log_deg)
ax.plot(degrees, fit_counts, color='black', linestyle='--', label=f'Fit (slope: {slope:.2f})')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Degree (log scale)', fontsize=10, fontweight='bold')
ax.set_ylabel('Frequency (log scale)', fontsize=10, fontweight='bold')
ax.set_title('Network Degree Distribution', fontsize=12, fontweight='bold')
ax.legend(frameon=False)

plt.tight_layout()
deg_out = OUT_DIR / "Network_Degree_Distribution.tiff"
plt.savefig(deg_out, dpi=300, format='tiff')
plt.close()
pass  # Execution logging removed for final release

