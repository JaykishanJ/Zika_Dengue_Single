import glob

scripts_dir = 'd:/Somenath_Zika_work/Somenath_Zika_work/scripts'
py_files = glob.glob(scripts_dir + '/**/*.py', recursive=True)

for f in py_files:
    with open(f, 'r', encoding='utf-8') as file:
        content = file.read()
    
    # Fix single quotes
    content = content.replace("'E:/Zika_Enrichment/Publication_Pipeline'", "str(Path(__file__).resolve().parent.parent.parent)")
    content = content.replace("'E:/Zika_antigravity/SC_RNA_Zika_DENV'", "str(Path(__file__).resolve().parent.parent.parent)")
    
    # Fix bbox_inches
    content = content.replace(", bbox_inches='tight'", "")
    content = content.replace(", bbox_inches=\"tight\"", "")

    # Fix biomarker/drug repurposing logic if significance_group is missing
    content = content.replace("deg_df[deg_df['significance_group'] == 'Both']['symbol']", "deg_df[deg_df['tier'] == 1]['symbol']")
    content = content.replace("All_Significant_DEGs.csv", "crossmodal_validated_progression_signature.csv")
    content = content.replace("Step03_differential_expression/DE_High_vs_Mock_DENV.csv", "Step06_crossmodal_validation/crossmodal_validated_progression_signature.csv")

    with open(f, 'w', encoding='utf-8') as file:
        file.write(content)
