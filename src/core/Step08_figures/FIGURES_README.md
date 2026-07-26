# Figures — descriptions (all 600 dpi PNG)

## Publication figures
- **Figure1_scRNA_landscape** — UMAP of the 2,069 Huh7 single cells, coloured by condition
  (Mock/DENV/ZIKV), infection state (Mock/Bystander/Low/High), log viral load, and Leiden cluster.
- **Figure2_scRNA_volcanoes** — single-cell differential expression, High-infection vs Mock
  (batch-matched), for DENV and ZIKV. Red = up, blue = down (padj<0.05, |log2FC|>0.58).
- **Figure3_convergent_core_heatmap** — the 12-gene single-cell convergent core (log2FC in
  DENV vs ZIKV); genes significant and same-direction in both viruses.
- **Figure4_progression_trends** — expression vs per-cell viral load for the top convergent
  progression genes, in both viruses (the ER-stress/ATF4 program rising with viral load).
- **Figure5_crossmodal_evidence_heatmap** — the 30 common genes (rows) x 5 datasets (columns);
  colour = direction (red up / blue down / white ns). Shows DNAJC3 & HYOU1 confirmed in every
  dataset; the macrophage column is mostly ns (that dataset suppresses transcription).
- **Figure6_pathway_enrichment** — enriched pathways of the up-signature (Unfolded Protein
  Response, ER protein processing, mTORC1, IRE1/XBP1, tRNA aminoacylation).
- **Figure7_bulk_volcanoes** — the three bulk RNA-seq DE analyses we computed (ZIKV macrophage,
  ZIKV neural, DENV blood); red = up, blue = down.
- **Figure8_venn_singlecell_vs_bulk** — 3-set Venn of single-cell convergent progression genes
  (513) vs ZIKV bulk DEGs vs DENV bulk DEGs. Centre = 70 (single-cell genes DE in both viruses'
  bulk); 30 of those are direction-concordant = the common signature.
- **Figure9_CrossModal_Signature_Venn** — 4-set Venn (ZIKV up/down x DENV up/down; each set =
  single-cell + bulk union). Convergent UP = 178, convergent DOWN = 149.

## Supplementary figures
- **FigureS1_QC_violin** — QC metrics (genes/cell, counts, % mito) with filter thresholds (2,260 -> 2,069 cells).
- **FigureS2_progression_concordance** — Spearman rho(DENV) vs rho(ZIKV) vs viral load for the 513 genes.
- **FigureS3_DEG_counts_per_dataset** — number of DE genes in single-cell (DENV, ZIKV) and each bulk cohort.
- **FigureS4_crossmodal_tiers** — distribution of the 513 genes across validation tiers (1-4).

## IMPORTANT — note on the "common gene" numbers (they are different views, all correct)
Different figures show different-sized "convergent" sets because they use different strictness:
- **30 genes** (Figure 5) = the rigorous common signature: significant + SAME direction in both
  viruses across single-cell AND bulk (a symmetric rule). This is the headline result.
- **70** (Figure 8, Venn centre) = single-cell convergent genes that are DE in both viruses' bulk,
  IGNORING direction. 30 of these 70 are direction-concordant.
- **178 up / 149 down** (Figure 9) = the loosest, union-based view (up in either modality of both
  viruses), matching the standard cross-modal-Venn method. Includes discordant/contradictory genes.

For the scientific write-up (objective, methods, findings) see **../PROJECT_REPORT.txt**.
