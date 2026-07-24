# Somenath_Zika_work — step-by-step reproducibility workspace

> **START HERE:**
> - **`PROJECT_REPORT.txt`** — the full write-up: objective, methods, findings, limitations (read this for the science).
> - **`Step08_figures/FIGURES_README.md`** — description of every figure.
> - This README (below) — how the pipeline is organised and how to re-run it.


Each analysis step is a **self-contained folder with its own `.py` script**. Run them **one by one, in order** — each reads the previous step's output and the raw data, and writes its own results. This regenerates the entire convergent Zika/Dengue progression analysis from scratch.

**Fully self-contained:** all raw data is included in `raw_data/` inside this folder (single-cell GSE110496.zip, the 3 bulk cohorts, and the Ensembl→symbol map). Nothing outside this folder is needed — zip it and anyone can reproduce the whole analysis. Only Steps 06–07 need internet (mygene / Enrichr).

```
raw_data/
├── single_cell_GSE110496/GSE110496.zip   (2,260 counts + SOFT)
├── reference/ensembl_to_symbol.csv
└── bulk/  GSE118305 (macrophage) + GSE78711 (neural) + GSE279208 (blood)
```

## How to run (one by one)
```
conda activate zika_sc      # scanpy, pydeseq2, leidenalg, igraph, gseapy, mygene, scipy
cd Step00_build_data          && python Step00_build_data.py
cd ../Step01_qc               && python Step01_qc.py
cd ../Step02_clustering       && python Step02_clustering.py
cd ../Step03_differential_expression && python Step03_differential_expression.py
cd ../Step04_core_signature   && python Step04_core_signature.py
cd ../Step05_progression      && python Step05_progression.py
cd ../Step05b_bulk_differential_expression && python Step05b_bulk_DE.py           # needs internet (mygene)
cd ../Step06_crossmodal_validation   && python Step06_crossmodal_validation.py
cd ../Step07_enrichment       && python Step07_enrichment.py                      # needs internet (Enrichr)
cd ../Step08_figures          && python Step08_figures.py
```
Each script uses paths **relative to its own location**, so the folder is portable.

## Steps & outputs
| Step | Script | Produces |
|---|---|---|
| 00 | build_data | `adata_raw.h5ad` (2,260 cells, from GSE110496.zip) |
| 01 | qc | `adata_qc.h5ad`, QC violin, summary |
| 02 | clustering | `adata_processed.h5ad`, UMAPs |
| 03 | differential_expression | DE tables (High vs Mock, batch-matched) + volcanoes |
| 04 | core_signature | 12-gene convergent core + heatmap |
| 05 | progression | 513 convergent progression genes + trends |
| 05b | **bulk_differential_expression** | **bulk RNA-seq DE (computed in-house) for all 3 cohorts** |
| 06 | crossmodal_validation | 30 common genes (Tier 1), reading Step 05b's bulk DE |
| 07 | enrichment | pathway enrichment (up/down) + barplot |
| 08 | figures | publication/ + supplementary/ (600 dpi PNG) |

## Reproducibility check (this run vs the original ForProfessor package)
Every step reproduced the original **exactly**:

| Step | This workspace | Original package | Match |
|---|---|---|---|
| Build | 2,260 × 60,716 | 2,260 × 60,716 | ✅ |
| QC | 2,260 → 2,069 cells | 2,260 → 2,069 | ✅ |
| Clustering | 11 Leiden clusters | 11 | ✅ |
| DE | DENV 134 (70↑/64↓), ZIKV 248 (184↑/64↓) | identical | ✅ |
| Convergent core | 12 genes (10↑, 2↓) — **identical gene set** | 12 | ✅ |
| Progression | 513 genes (273↑, 240↓) — **identical gene set** | 513 | ✅ |
| Cross-modal | GOLD 2, Tier 1 = 30 (**identical gene set**), Tier 2 = 138 | same | ✅ |
| Enrichment | 168 up-terms (**identical**), 10 down | same | ✅ |

Gene-set identity was verified programmatically (core, progression, Tier-1 common, enrichment terms all identical). **The pipeline is fully reproducible.**
