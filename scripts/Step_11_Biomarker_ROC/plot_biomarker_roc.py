import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_curve, auc
import mygene

sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))
import pub_style

pub_style.apply_style()

BASE_DIR = Path("E:/Zika_antigravity/SC_RNA_Zika_DENV")
OUT_DIR = Path("E:/Zika_Enrichment/Publication_Pipeline/Step_11_Biomarker_ROC")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def main():
    pass  # Execution logging removed for final release
    deg_df = pd.read_csv(BASE_DIR / "03_results/phase4_All_DEG/All_Significant_DEGs.csv")
    core_genes = deg_df[deg_df['significance_group'] == 'Both']['symbol'].dropna().unique().tolist()
    pass  # Execution logging removed for final release
    
    pass  # Execution logging removed for final release
    mg = mygene.MyGeneInfo()
    res = mg.querymany(core_genes, scopes='symbol', fields='entrezgene', species='human', verbose=False)
    entrez_ids = [str(r.get('entrezgene')) for r in res if 'entrezgene' in r]
    
    pass  # Execution logging removed for final release
    expr_df = pd.read_csv(BASE_DIR / "Bulk_RNA_seq/GSE279208/Results/Counts/GSE279208_VST_Normalized_Matrix.csv", index_col=0)
    pheno_df = pd.read_csv(BASE_DIR / "Bulk_RNA_seq/GSE279208/Results/Annotation/GSE279208_Final_Phenotype.csv")
    
    # Ensure expression matrix indices are strings
    expr_df.index = expr_df.index.astype(str)
    
    # Filter for core genes that exist in the validation cohort
    available_genes = list(set(entrez_ids) & set(expr_df.index))
    pass  # Execution logging removed for final release
    
    X_full = expr_df.loc[available_genes].T
    
    # Align phenotypes
    pheno_df.set_index('SampleID', inplace=True)
    pheno_df = pheno_df.loc[X_full.index]
    
    # Create binary target: healthy=0, Dengue/Severe_Dengue=1
    y = (pheno_df['Group'] != 'healthy').astype(int).values
    X = X_full.values
    
    pass  # Execution logging removed for final release
    loo = LeaveOneOut()
    y_true = []
    y_prob = []
    
    clf = LogisticRegression(penalty='l2', max_iter=1000, random_state=42)
    
    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        
        clf.fit(X_train, y_train)
        y_prob.append(clf.predict_proba(X_test)[0, 1])
        y_true.append(y_test[0])
        
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    roc_auc = auc(fpr, tpr)
    
    pass  # Execution logging removed for final release
    
    # Plot ROC Curve
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(fpr, tpr, color='#E63946', lw=2.5, label=f'Core Signature (AUC = {roc_auc:.3f})')
    ax.plot([0, 1], [0, 1], color='gray', lw=1.5, linestyle='--')
    
    ax.set_xlim([-0.02, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate', fontsize=12, fontweight='bold')
    ax.set_ylabel('True Positive Rate', fontsize=12, fontweight='bold')
    ax.set_title('Diagnostic Performance of Core Signature\nValidation Cohort (GSE279208)', fontsize=14, fontweight='bold', pad=15)
    ax.legend(loc="lower right", frameon=False, fontsize=12)
    
    pub_style.save_pub_fig(fig, OUT_DIR / "Biomarker_ROC_Curve.tiff", width_mm=120, height_mm=120)
    plt.close()

if __name__ == "__main__":
    main()
