import glob
import re
from pathlib import Path

scripts_dir = 'd:/Somenath_Zika_work/Somenath_Zika_work/scripts'
py_files = glob.glob(scripts_dir + '/**/*.py', recursive=True)

for f in py_files:
    with open(f, 'r', encoding='utf-8') as file:
        content = file.read()
    
    # sys.path replacements to find utils
    content = content.replace('sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))', 'sys.path.append(str(Path(__file__).resolve().parent.parent.parent / "utils"))')
    content = content.replace('if str(Path("E:/Zika_Enrichment/Publication_Pipeline")) not in sys.path:\n    sys.path.append(str(Path("E:/Zika_Enrichment/Publication_Pipeline")))', 'if str(Path(__file__).resolve().parent.parent.parent / "utils") not in sys.path:\n    sys.path.append(str(Path(__file__).resolve().parent.parent.parent / "utils"))')

    # Replace absolute path bases with dynamic ROOT
    content = content.replace('Path("E:/Zika_Enrichment/Publication_Pipeline")', 'Path(__file__).resolve().parent.parent.parent')
    content = content.replace('Path("E:/Zika_antigravity/SC_RNA_Zika_DENV")', 'Path(__file__).resolve().parent.parent.parent')
    
    # Some scripts just use strings:
    content = content.replace('"E:/Zika_Enrichment/Publication_Pipeline"', 'str(Path(__file__).resolve().parent.parent.parent)')
    content = content.replace('"E:/Zika_antigravity/SC_RNA_Zika_DENV"', 'str(Path(__file__).resolve().parent.parent.parent)')
    
    # Fix paths
    content = content.replace('final_publication_figures', 'publication_figures')
    content = content.replace('Bulk_RNA_seq', 'Step05b_bulk_differential_expression')
    content = content.replace('03_results/phase4_All_DEG/All_Significant_DEGs.csv', 'Step03_differential_expression/DE_High_vs_Mock_DENV.csv')
    
    # specific bulk logic fixes (bulk DE output format changed from log2FoldChange/padj/Gene to lfc/padj/SYM/dir)
    # ZIKV_positive_vs_Mock_All_Genes.csv -> bulk_GSE118305_ZIKV_macrophage_DE.csv
    content = content.replace('GSE118305/Results/Tables/ZIKV_positive_vs_Mock_All_Genes.csv', 'bulk_GSE118305_ZIKV_macrophage_DE.csv')
    content = content.replace('GSE279208/Results/Tables/Dengue_vs_healthy_All_Genes.csv', 'bulk_GSE279208_DENV_blood_DE.csv')
    
    if 'log2FoldChange' in content:
        content = content.replace("'log2FoldChange'", "'lfc'")
        content = content.replace('"log2FoldChange"', '"lfc"')
    if "df['Symbol'] = df['Gene']" in content:
        content = content.replace("df['Symbol'] = df['Gene']", "df['Symbol'] = df['SYM']")
    if "dropna(subset=['Gene'])" in content:
        content = content.replace("dropna(subset=['Gene'])", "dropna(subset=['SYM'])")
        
    content = content.replace("df['Gene']", "df['SYM']")
    
    with open(f, 'w', encoding='utf-8') as file:
        file.write(content)

print("Updated all scripts.")
