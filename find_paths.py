import glob
import re

py_files = glob.glob('d:/Somenath_Zika_work/Somenath_Zika_work/scripts/**/*.py', recursive=True)
paths_found = set()

for f in py_files:
    with open(f, 'r', encoding='utf-8') as file:
        content = file.read()
        paths = re.findall(r'Path\([\"\']([^\"\']+)[\"\']\)', content)
        for p in paths:
            paths_found.add(p)

for p in sorted(list(paths_found)):
    print(p)
