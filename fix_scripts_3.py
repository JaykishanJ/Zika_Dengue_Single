import glob
import re

scripts_dir = 'd:/Somenath_Zika_work/Somenath_Zika_work/scripts'
py_files = glob.glob(scripts_dir + '/**/*.py', recursive=True)

for f in py_files:
    with open(f, 'r', encoding='utf-8') as file:
        content = file.read()
    
    # 1. Add numpy if missing
    if 'import numpy' not in content:
        # Add it after pandas or at the top
        content = "import numpy as np\n" + content
        
    # 2. Fix the raw counts path for biomarker ROC (and any others)
    # The old bad path is likely "Step05b_bulk_differential_expression/GSE279208/Results/Counts/GSE279208_VST_Normalized_Matrix.csv"
    # or similar. I'll replace the prefix.
    content = content.replace("Step05b_bulk_differential_expression/GSE279208/Results/Counts/GSE279208_VST_Normalized_Matrix.csv", 
                              "raw_data/bulk/GSE279208/GSE279208_VST_Normalized_Matrix.csv")
    content = content.replace("Bulk_RNA_seq/GSE279208/Results/Counts/GSE279208_VST_Normalized_Matrix.csv", 
                              "raw_data/bulk/GSE279208/GSE279208_VST_Normalized_Matrix.csv")
                              
    # 3. Fix save_pub_fig facecolor argument
    # Some scripts have `save_pub_fig(..., facecolor='white')`
    content = re.sub(r",\s*facecolor=['\"][^'\"]+['\"]", "", content)

    # 4. Any missing imports (just in case)
    if 'pd.' in content and 'import pandas' not in content:
        content = "import pandas as pd\n" + content

    with open(f, 'w', encoding='utf-8') as file:
        file.write(content)

print("Applied final fixes to visualization scripts!")
