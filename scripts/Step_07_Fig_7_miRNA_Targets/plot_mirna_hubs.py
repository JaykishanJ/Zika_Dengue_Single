import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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
    
    # Sort by the number of targets
    df['Target_Count'] = df['Overlapping_mRNA_Genes'].apply(lambda x: len(str(x).split(',')) if pd.notna(x) else 0)
    df_sorted = df.sort_values('Target_Count', ascending=False).head(15)
    
    fig, ax = plt.subplots(figsize=(4.72, 3.93))
    
    sns.barplot(data=df_sorted, x='Target_Count', y='miRNA', palette='Blues_r', edgecolor='black', linewidth=0.5, ax=ax)
    
    ax.set_title("Top miRNA Regulatory Hubs", fontweight='bold', pad=10)
    ax.set_xlabel("Number of Targeted Core Genes")
    ax.set_ylabel("")
    
    pub_style.save_pub_fig(fig, OUT_DIR / "Top_miRNA_Hubs.tiff", width_mm=120, height_mm=100)
    plt.close()

if __name__ == "__main__":
    main()
