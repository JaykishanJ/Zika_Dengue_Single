import pandas as pd
from pathlib import Path

OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_07_Fig_7_miRNA_Targets")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 1. Load the 98 core genes
pass  # Execution logging removed for final release
core_df = pd.read_csv("E:/Zika_antigravity/SC_RNA_Zika_DENV/03_results/phase3_shared_degs/target_DEGs_annotated.csv")
core_genes = set(core_df['symbol'].dropna())
pass  # Execution logging removed for final release

# 2. Load the miRNA target datasets
pass  # Execution logging removed for final release
zika_targets_df = pd.read_csv("E:/Zika_Enrichment/System_zika_mirna/02_literature_resources/zika_set_targets.csv")
zika_targets_df['Virus'] = 'ZIKV'

dengue_targets_df = pd.read_csv("E:/Zika_Enrichment/System_zika_mirna/02_literature_resources/dengue_full_set_targets.csv")
dengue_targets_df['Virus'] = 'DENV'

# Combine them for processing
all_targets_df = pd.concat([zika_targets_df, dengue_targets_df], ignore_index=True)
pass  # Execution logging removed for final release

# 3. Intersect with 98 core genes
pass  # Execution logging removed for final release
overlapping_df = all_targets_df[all_targets_df['Target Gene'].isin(core_genes)].copy()
pass  # Execution logging removed for final release

# 4. Format the output: list of miRNA with respective mRNA overlapping gene
pass  # Execution logging removed for final release
# Group by miRNA and Virus, then collect unique target genes
results = []
grouped = overlapping_df.groupby(['miRNA', 'Virus'])
for (mirna, virus), group in grouped:
    # Get unique overlapping mRNA genes targeted by this miRNA
    targets = sorted(list(set(group['Target Gene'])))
    target_count = len(targets)
    targets_str = ", ".join(targets)
    results.append({
        'Virus': virus,
        'miRNA': mirna,
        'Overlapping_mRNA_Count': target_count,
        'Overlapping_mRNA_Genes': targets_str
    })

# Convert to DataFrame and sort by number of targeted genes
final_df = pd.DataFrame(results)
final_df = final_df.sort_values(by='Overlapping_mRNA_Count', ascending=False)

# Save to CSV
out_path = OUT_DIR / "Core_98_miRNA_mRNA_Regulatory_Pairs.csv"
final_df.to_csv(out_path, index=False)
pass  # Execution logging removed for final release

pass  # Execution logging removed for final release
import matplotlib.pyplot as plt

# Nature SR Aesthetics
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

from venn import venn
from matplotlib.patches import Patch

zika_targets_core = set(overlapping_df[overlapping_df['Virus'] == 'ZIKV']['Target Gene'])
dengue_targets_core = set(overlapping_df[overlapping_df['Virus'] == 'DENV']['Target Gene'])

labels = {
    "ZIKV miRNA Targets": zika_targets_core,
    "DENV miRNA Targets": dengue_targets_core
}

# 180x160mm -> 7.08x6.30 inches
fig, ax = plt.subplots(figsize=(7.08, 6.30))
venn(labels, ax=ax, cmap="Set3")

if ax.get_legend():
    ax.get_legend().remove()

cmap = plt.get_cmap("Set3")
handles = [Patch(color=cmap(i), label=lbl) for i, lbl in enumerate(labels.keys())]
ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.05, 0.5), ncol=1, frameon=False, fontsize=10)

plt.title("Core DEGs Targeted by Viral miRNAs", fontsize=14)
plt.tight_layout()
plt.subplots_adjust(right=0.75) 
plt.savefig(OUT_DIR / "miRNA_Target_Venn.tiff", dpi=300, format='tiff', bbox_inches='tight')
plt.close()
pass  # Execution logging removed for final release

# Print a quick summary of the top 5 miRNAs regulating the core signature
pass  # Execution logging removed for final release
pass  # Execution logging removed for final release
