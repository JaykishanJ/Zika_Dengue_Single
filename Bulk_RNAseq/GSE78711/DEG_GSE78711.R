################################################################################
# FINAL PUBLICATION-READY SCRIPT
# RNA-SEQ PROCESSED DIFFERENTIAL EXPRESSION ANALYSIS
#
# Dataset: GSE78711
#
# Study:
# Zika Virus Targets Human Cortical Neural Precursors and Attenuates Their Growth
#
# Important:
# GEO contains 8 samples across 2 GEO matrix/platform objects:
#   Mock1-1, Mock2-1, ZIKV1-1, ZIKV2-1
#   Mock1-2, Mock2-2, ZIKV1-2, ZIKV2-2
#
# Supplementary file:
#   GSE78711_gene.exp.all.txt.gz
#
# This file is NOT an 8-sample raw count matrix.
# It is already a processed Mock vs ZIKV DEG summary table with columns:
#   gene, sample_1, sample_2, value_1, value_2,
#   log2.fold_change., p_value, significant
#
# Therefore:
#   1. This script correctly builds phenotype for all 8 GEO samples.
#   2. It analyzes the processed DEG table directly.
#   3. It does NOT falsely treat value_1/value_2/logFC/p_value as samples.
#   4. PCA/sample correlation/sample clustering is NOT performed here.
#
# Publication-ready outputs:
#   - DEG summary table
#   - P-value histogram
#   - FDR histogram
#   - log2FC density plot
#   - DEG count barplot
#   - Mock vs ZIKV expression scatter
#   - MA plot
#   - Volcano plot
#   - EnhancedVolcano plot
#   - Top up/down DEG barplot
#   - Top DEG lollipop plot
#   - Top 50 DEG summary heatmap
#   - Ranked log2FC profile
#   - GO BP/MF/CC enrichment
#   - KEGG enrichment
#   - Publication-ready PDF + 600 dpi PNG outputs
################################################################################

rm(list = ls())
gc()

# ==============================================================================
# 0. USER SETTINGS
# ==============================================================================

WORKDIR <- "E:/zika_virus/GSE78711"
GSE_ID  <- "GSE78711"

setwd(WORKDIR)

dir.create("Results", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Plots", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/QC", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Tables", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Annotation", showWarnings = FALSE, recursive = TRUE)

cat("\n============================================================\n")
cat("Processed RNA-seq DEG Analysis:", GSE_ID, "\n")
cat("Working directory:", getwd(), "\n")
cat("============================================================\n\n")


# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================

cat("[1] Loading required libraries...\n")

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "readxl",
  "scales",
  "RColorBrewer"
)

bioc_pkgs <- c(
  "GEOquery",
  "Biobase",
  "EnhancedVolcano",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot"
)

install_if_missing <- function(pkgs, bioc = FALSE) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing missing package:", pkg, "\n")
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg)
      }
    }
  }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(readxl)
  library(scales)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
})

theme_set(theme_bw(base_size = 14))

cat("✓ Libraries loaded successfully.\n\n")


# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

safe_read_table <- function(file) {
  
  cat("Trying file:", file, "\n")
  
  if (grepl("\\.csv(\\.gz)?$", file, ignore.case = TRUE)) {
    
    if (grepl("\\.gz$", file, ignore.case = TRUE)) {
      con <- gzfile(file, "rt")
      on.exit(close(con), add = TRUE)
      x <- read.csv(
        con,
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      x <- read.csv(
        file,
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
    
  } else {
    
    if (grepl("\\.gz$", file, ignore.case = TRUE)) {
      con <- gzfile(file, "rt")
      on.exit(close(con), add = TRUE)
      x <- read.table(
        con,
        header = TRUE,
        sep = "\t",
        quote = "",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      x <- read.table(
        file,
        header = TRUE,
        sep = "\t",
        quote = "",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
  
  return(x)
}


clean_sample_names <- function(x) {
  x <- gsub("^X", "", x)
  x <- gsub("\\.bam$|\\.fastq$|\\.fq$|\\.txt$|\\.gz$", "", x, ignore.case = TRUE)
  x <- gsub("\\.", "-", x)
  x <- gsub("_$", "", x)
  x <- gsub("[^A-Za-z0-9_-]", "", x)
  x
}


make_safe_pdf <- function(filename, expr, width = 8, height = 6) {
  
  pdf(filename, width = width, height = height)
  
  tryCatch(
    expr,
    error = function(e) {
      plot.new()
      text(
        0.5,
        0.5,
        paste("Plot failed:\n", e$message),
        cex = 1
      )
      message("Plot failed: ", filename, " | ", e$message)
    },
    finally = dev.off()
  )
}


save_pub_plot <- function(plot_object, filename_base, width = 8, height = 6, dpi = 600) {
  
  pdf_file <- paste0(filename_base, ".pdf")
  png_file <- paste0(filename_base, ".png")
  
  ggsave(
    filename = pdf_file,
    plot = plot_object,
    width = width,
    height = height,
    device = cairo_pdf
  )
  
  ggsave(
    filename = png_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi
  )
}


pub_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6)
  )


# ==============================================================================
# 3. DOWNLOAD AND COMBINE GEO METADATA FROM BOTH MATRIX OBJECTS
# ==============================================================================

cat("[2] Downloading GEO metadata...\n")

options(timeout = 300000)

gset <- getGEO(GSE_ID, GSEMatrix = TRUE, getGPL = FALSE)

cat("Number of GEO matrix objects found:", length(gset), "\n")

metadata_list <- lapply(seq_along(gset), function(i) {
  p <- pData(gset[[i]])
  p$Platform <- annotation(gset[[i]])
  p$GEO_matrix_index <- i
  p
})

metadata <- dplyr::bind_rows(metadata_list)

write.csv(
  metadata,
  "Results/Annotation/All_GEO_Metadata_8_Total_Samples.csv",
  row.names = FALSE
)

cat("Total GEO samples after combining all matrix objects:", nrow(metadata), "\n\n")

phenotype <- data.frame(
  geo_accession = metadata$geo_accession,
  title = metadata$title,
  Platform = metadata$Platform,
  GEO_matrix_index = metadata$GEO_matrix_index,
  stringsAsFactors = FALSE
)

phenotype$SampleID <- clean_sample_names(phenotype$title)

phenotype$Condition <- dplyr::case_when(
  grepl("mock", phenotype$title, ignore.case = TRUE) ~ "Mock",
  grepl("zikv|zika", phenotype$title, ignore.case = TRUE) ~ "ZIKV",
  TRUE ~ NA_character_
)

phenotype$Replicate_Set <- dplyr::case_when(
  grepl("-1$", phenotype$title) ~ "Set_1",
  grepl("-2$", phenotype$title) ~ "Set_2",
  TRUE ~ NA_character_
)

phenotype$Condition <- factor(
  phenotype$Condition,
  levels = c("Mock", "ZIKV")
)

write.csv(
  phenotype,
  "Results/Annotation/GSE78711_Final_8_Sample_Phenotype.csv",
  row.names = FALSE
)

cat("Final 8-sample phenotype table:\n")
print(phenotype)

cat("\nPhenotype summary:\n")
print(table(phenotype$Condition, useNA = "ifany"))

cat("\nReplicate-set summary:\n")
print(table(phenotype$Replicate_Set, phenotype$Condition, useNA = "ifany"))

if (nrow(phenotype) != 8) {
  warning("Expected 8 GEO samples, but detected: ", nrow(phenotype))
}

if (any(is.na(phenotype$Condition))) {
  stop(
    "Some samples could not be assigned to Mock/ZIKV. Check Results/Annotation/GSE78711_Final_8_Sample_Phenotype.csv"
  )
}

cat("\n✓ GEO phenotype correctly built from both platforms.\n\n")


# ==============================================================================
# 4. DOWNLOAD SUPPLEMENTARY FILES
# ==============================================================================

cat("[3] Downloading supplementary files...\n")

supp_files <- getGEOSuppFiles(
  GSE_ID,
  makeDirectory = TRUE,
  baseDir = getwd()
)

cat("\nSupplementary files reported by GEO:\n")
print(rownames(supp_files))
cat("\n")

supp_dir <- file.path(getwd(), GSE_ID)

all_files <- list.files(
  supp_dir,
  recursive = TRUE,
  full.names = TRUE
)

write.csv(
  data.frame(
    File = all_files,
    Basename = basename(all_files)
  ),
  "Results/Annotation/All_Supplementary_Files_Found.csv",
  row.names = FALSE
)

cat("All supplementary files found:\n")
print(basename(all_files))
cat("\n")


# ==============================================================================
# 5. OPTIONAL: READ TABLE S2A READ SUMMARY XLSX
# ==============================================================================

cat("[4] Reading Table S2A read summary if available...\n")

xlsx_file <- all_files[
  grepl("GSE78711_Table.S2A.reads_R1.xlsx$", basename(all_files), ignore.case = TRUE)
]

if (length(xlsx_file) > 0) {
  
  xlsx_file <- xlsx_file[1]
  cat("Found XLSX file:", xlsx_file, "\n")
  
  sheets <- readxl::excel_sheets(xlsx_file)
  
  for (sh in sheets) {
    
    cat("Reading sheet:", sh, "\n")
    
    tmp <- tryCatch(
      readxl::read_excel(xlsx_file, sheet = sh),
      error = function(e) NULL
    )
    
    if (!is.null(tmp)) {
      
      out_name <- paste0(
        "Results/Annotation/GSE78711_Table_S2A_",
        gsub("[^A-Za-z0-9]", "_", sh),
        ".csv"
      )
      
      write.csv(
        as.data.frame(tmp),
        out_name,
        row.names = FALSE
      )
    }
  }
  
  cat("✓ Table S2A was saved into Results/Annotation.\n")
  cat("Note: This file is a read/sample information table, not necessarily a gene-level count matrix.\n\n")
  
} else {
  
  cat("No Table S2A XLSX file found.\n\n")
}


# ==============================================================================
# 6. FIND PROCESSED DEG TABLE
# ==============================================================================

cat("[5] Finding processed DEG table...\n")

deg_candidates <- all_files[
  grepl("GSE78711_gene.exp.all.txt.gz$", basename(all_files), ignore.case = TRUE)
]

if (length(deg_candidates) == 0) {
  stop("Could not find GSE78711_gene.exp.all.txt.gz")
}

deg_file <- deg_candidates[1]

cat("Selected processed DEG table:\n")
cat(deg_file, "\n\n")


# ==============================================================================
# 7. LOAD PROCESSED DEG TABLE
# ==============================================================================

cat("[6] Loading processed DEG table...\n")

deg_raw <- safe_read_table(deg_file)

cat("Raw processed DEG table dimensions:\n")
cat("Rows:", nrow(deg_raw), "\n")
cat("Columns:", ncol(deg_raw), "\n\n")

cat("Column names:\n")
print(colnames(deg_raw))
cat("\n")

cat("First rows:\n")
print(head(deg_raw))
cat("\n")

write.csv(
  data.frame(Column_Names = colnames(deg_raw)),
  "Results/Annotation/Processed_DEG_Table_Column_Names.csv",
  row.names = FALSE
)

expected_cols <- c(
  "gene",
  "sample_1",
  "sample_2",
  "value_1",
  "value_2",
  "log2.fold_change.",
  "p_value",
  "significant"
)

missing_cols <- setdiff(expected_cols, colnames(deg_raw))

if (length(missing_cols) > 0) {
  stop(
    "Processed DEG table does not have expected columns. Missing: ",
    paste(missing_cols, collapse = ", ")
  )
}

cat("✓ Confirmed: file is a processed Mock vs ZIKV DEG summary table.\n")
cat("  It will NOT be treated as sample-level expression data.\n\n")


# ==============================================================================
# 8. FORMAT DEG TABLE
# ==============================================================================

cat("[7] Formatting processed DEG table...\n")

deg <- deg_raw %>%
  dplyr::rename(
    Gene = gene,
    Group_1 = sample_1,
    Group_2 = sample_2,
    Mock_expression = value_1,
    ZIKV_expression = value_2,
    log2FoldChange = log2.fold_change.,
    pvalue = p_value,
    Original_significant_call = significant
  )

deg$Mock_expression <- as.numeric(deg$Mock_expression)
deg$ZIKV_expression <- as.numeric(deg$ZIKV_expression)
deg$log2FoldChange <- as.numeric(deg$log2FoldChange)
deg$pvalue <- as.numeric(deg$pvalue)

deg$padj <- p.adjust(deg$pvalue, method = "BH")

deg <- deg %>%
  dplyr::mutate(
    Regulation = dplyr::case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >= 1  ~ "Up_in_ZIKV",
      !is.na(padj) & padj < 0.05 & log2FoldChange <= -1 ~ "Down_in_ZIKV",
      TRUE ~ "Not_significant"
    )
  ) %>%
  dplyr::arrange(padj, pvalue)

write.csv(
  deg,
  "Results/Tables/GSE78711_Processed_DEG_ZIKV_vs_Mock_All_Genes.csv",
  row.names = FALSE
)

deg_p005 <- deg %>%
  dplyr::filter(!is.na(pvalue)) %>%
  dplyr::filter(pvalue < 0.05)

deg_fdr005 <- deg %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::filter(padj < 0.05)

deg_fdr005_logfc1 <- deg %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) >= 1)

deg_up <- deg_fdr005_logfc1 %>%
  dplyr::filter(log2FoldChange >= 1)

deg_down <- deg_fdr005_logfc1 %>%
  dplyr::filter(log2FoldChange <= -1)

write.csv(
  deg_p005,
  "Results/Tables/GSE78711_DEG_pvalue_0.05.csv",
  row.names = FALSE
)

write.csv(
  deg_fdr005,
  "Results/Tables/GSE78711_DEG_FDR_0.05.csv",
  row.names = FALSE
)

write.csv(
  deg_fdr005_logfc1,
  "Results/Tables/GSE78711_DEG_FDR_0.05_log2FC_1.csv",
  row.names = FALSE
)

write.csv(
  deg_up,
  "Results/Tables/GSE78711_DEG_UP_in_ZIKV_FDR_0.05_log2FC_1.csv",
  row.names = FALSE
)

write.csv(
  deg_down,
  "Results/Tables/GSE78711_DEG_DOWN_in_ZIKV_FDR_0.05_log2FC_1.csv",
  row.names = FALSE
)

summary_table <- data.frame(
  Dataset = GSE_ID,
  Contrast = "ZIKV_vs_Mock",
  Total_genes = nrow(deg),
  Pvalue_0.05 = nrow(deg_p005),
  FDR_0.05 = nrow(deg_fdr005),
  FDR_0.05_log2FC_1 = nrow(deg_fdr005_logfc1),
  Up_in_ZIKV = nrow(deg_up),
  Down_in_ZIKV = nrow(deg_down)
)

write.csv(
  summary_table,
  "Results/Tables/GSE78711_DEG_Summary.csv",
  row.names = FALSE
)

cat("DEG summary:\n")
print(summary_table)
cat("\n")


# ==============================================================================
# 9. PUBLICATION-READY DEG PLOTS
# ==============================================================================

cat("[8] Generating publication-ready DEG plots...\n")

deg <- deg %>%
  dplyr::mutate(
    minusLog10P = -log10(pvalue),
    minusLog10FDR = -log10(padj),
    Mean_expression = rowMeans(
      cbind(Mock_expression, ZIKV_expression),
      na.rm = TRUE
    ),
    log2Mean_expression = log2(Mean_expression + 1),
    Direction = dplyr::case_when(
      padj < 0.05 & log2FoldChange >= 1  ~ "Up in ZIKV",
      padj < 0.05 & log2FoldChange <= -1 ~ "Down in ZIKV",
      TRUE ~ "Not significant"
    )
  )

deg$minusLog10P[is.infinite(deg$minusLog10P)] <- max(
  deg$minusLog10P[is.finite(deg$minusLog10P)],
  na.rm = TRUE
)

deg$minusLog10FDR[is.infinite(deg$minusLog10FDR)] <- max(
  deg$minusLog10FDR[is.finite(deg$minusLog10FDR)],
  na.rm = TRUE
)

deg$Direction <- factor(
  deg$Direction,
  levels = c("Up in ZIKV", "Down in ZIKV", "Not significant")
)

deg_colors <- c(
  "Up in ZIKV" = "#D73027",
  "Down in ZIKV" = "#4575B4",
  "Not significant" = "grey75"
)

# Automatically label top genes by FDR
top_label <- deg %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 20)

# Optional biologically relevant ZIKV/interferon showcase genes
showcase_genes <- c(
  "IFNB1", "IFIT1", "IFIT2", "IFIT3",
  "ISG15", "MX1", "OAS1", "OAS2", "OAS3",
  "CXCL10", "DDX58", "IFIH1", "STAT1",
  "STAT2", "IRF7", "RSAD2", "IFI44", "IFI44L"
)

showcase_label <- deg %>%
  dplyr::filter(Gene %in% showcase_genes)

label_df <- dplyr::bind_rows(top_label, showcase_label) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

write.csv(
  label_df,
  "Results/Tables/GSE78711_Genes_Labelled_In_Publication_Plots.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9.1 DEG count barplot
# ------------------------------------------------------------------------------

deg_count_df <- deg %>%
  dplyr::count(Direction) %>%
  dplyr::mutate(Direction = factor(Direction, levels = levels(deg$Direction)))

p_deg_count <- ggplot(deg_count_df, aes(x = Direction, y = n, fill = Direction)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.4, fontface = "bold", size = 4) +
  scale_fill_manual(values = deg_colors) +
  labs(
    title = "Differentially expressed genes",
    subtitle = "FDR < 0.05 and |log2FC| ≥ 1",
    x = NULL,
    y = "Number of genes"
  ) +
  pub_theme +
  theme(legend.position = "none")

save_pub_plot(
  p_deg_count,
  "Results/Plots/01_Publication_DEG_Count_Barplot",
  width = 7,
  height = 6
)

# ------------------------------------------------------------------------------
# 9.2 P-value distribution
# ------------------------------------------------------------------------------

p_pvalue_hist <- ggplot(deg, aes(x = pvalue)) +
  geom_histogram(
    bins = 60,
    fill = "#4DAF4A",
    color = "white"
  ) +
  labs(
    title = "P-value distribution",
    subtitle = "Processed ZIKV vs Mock DEG table",
    x = "P-value",
    y = "Gene count"
  ) +
  pub_theme

save_pub_plot(
  p_pvalue_hist,
  "Results/QC/02_Publication_Pvalue_Distribution",
  width = 7,
  height = 5
)

# ------------------------------------------------------------------------------
# 9.3 FDR distribution
# ------------------------------------------------------------------------------

p_fdr_hist <- ggplot(deg, aes(x = padj)) +
  geom_histogram(
    bins = 60,
    fill = "#377EB8",
    color = "white"
  ) +
  labs(
    title = "FDR distribution",
    subtitle = "Benjamini-Hochberg adjusted p-values",
    x = "Adjusted p-value / FDR",
    y = "Gene count"
  ) +
  pub_theme

save_pub_plot(
  p_fdr_hist,
  "Results/QC/03_Publication_FDR_Distribution",
  width = 7,
  height = 5
)

# ------------------------------------------------------------------------------
# 9.4 log2 fold-change density
# ------------------------------------------------------------------------------

p_logfc_density <- ggplot(deg, aes(x = log2FoldChange)) +
  geom_density(fill = "#984EA3", alpha = 0.45, linewidth = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.5) +
  labs(
    title = "Distribution of log2 fold change",
    subtitle = "Dashed lines indicate |log2FC| = 1",
    x = "log2 fold change",
    y = "Density"
  ) +
  pub_theme

save_pub_plot(
  p_logfc_density,
  "Results/QC/04_Publication_log2FC_Density",
  width = 7,
  height = 5
)

# ------------------------------------------------------------------------------
# 9.5 Mock vs ZIKV expression scatter
# ------------------------------------------------------------------------------

scatter_df <- deg %>%
  dplyr::mutate(
    log2_Mock_expression = log2(Mock_expression + 1),
    log2_ZIKV_expression = log2(ZIKV_expression + 1)
  )

label_scatter_df <- label_df %>%
  dplyr::mutate(
    log2_Mock_expression = log2(Mock_expression + 1),
    log2_ZIKV_expression = log2(ZIKV_expression + 1)
  )

p_scatter <- ggplot(
  scatter_df,
  aes(
    x = log2_Mock_expression,
    y = log2_ZIKV_expression,
    color = Direction
  )
) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.7) +
  scale_color_manual(values = deg_colors) +
  ggrepel::geom_text_repel(
    data = label_scatter_df,
    aes(label = Gene),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  labs(
    title = "Mock vs ZIKV expression summary",
    subtitle = "Processed group-level expression values",
    x = "log2(Mock expression + 1)",
    y = "log2(ZIKV expression + 1)",
    color = "Regulation"
  ) +
  pub_theme

save_pub_plot(
  p_scatter,
  "Results/Plots/05_Publication_Mock_vs_ZIKV_Expression_Scatter",
  width = 8,
  height = 7
)

# ------------------------------------------------------------------------------
# 9.6 MA plot
# ------------------------------------------------------------------------------

p_ma <- ggplot(
  deg,
  aes(
    x = log2Mean_expression,
    y = log2FoldChange,
    color = Direction
  )
) +
  geom_point(alpha = 0.75, size = 1.7) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", linewidth = 0.6) +
  scale_color_manual(values = deg_colors) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = Gene),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  labs(
    title = "MA plot: ZIKV vs Mock",
    subtitle = "Mean expression versus log2 fold change",
    x = "log2 mean expression",
    y = "log2 fold change",
    color = "Regulation"
  ) +
  pub_theme

save_pub_plot(
  p_ma,
  "Results/Plots/06_Publication_MA_Plot_ZIKV_vs_Mock",
  width = 8,
  height = 7
)

# ------------------------------------------------------------------------------
# 9.7 Volcano plot
# ------------------------------------------------------------------------------

p_volcano <- ggplot(
  deg,
  aes(
    x = log2FoldChange,
    y = minusLog10FDR,
    color = Direction
  )
) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.7) +
  scale_color_manual(values = deg_colors) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = Gene),
    size = 3.2,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  labs(
    title = "Volcano plot: ZIKV vs Mock",
    subtitle = "FDR < 0.05 and |log2FC| ≥ 1",
    x = "log2 fold change",
    y = "-log10 FDR",
    color = "Regulation"
  ) +
  pub_theme

save_pub_plot(
  p_volcano,
  "Results/Plots/07_Publication_Volcano_ZIKV_vs_Mock",
  width = 9,
  height = 7
)

# ------------------------------------------------------------------------------
# 9.8 EnhancedVolcano
# ------------------------------------------------------------------------------

make_safe_pdf(
  "Results/Plots/08_Publication_EnhancedVolcano_ZIKV_vs_Mock.pdf",
  {
    print(
      EnhancedVolcano(
        deg,
        lab = deg$Gene,
        x = "log2FoldChange",
        y = "padj",
        selectLab = label_df$Gene,
        title = "GSE78711: ZIKV vs Mock",
        subtitle = "Processed GEO DEG summary table",
        caption = "Cutoff: FDR < 0.05 and |log2FC| ≥ 1",
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = 2.4,
        labSize = 3.5,
        colAlpha = 0.85,
        legendPosition = "bottom",
        drawConnectors = TRUE,
        widthConnectors = 0.5
      )
    )
  },
  width = 10,
  height = 8
)

png(
  "Results/Plots/08_Publication_EnhancedVolcano_ZIKV_vs_Mock.png",
  width = 10,
  height = 8,
  units = "in",
  res = 600
)

print(
  EnhancedVolcano(
    deg,
    lab = deg$Gene,
    x = "log2FoldChange",
    y = "padj",
    selectLab = label_df$Gene,
    title = "GSE78711: ZIKV vs Mock",
    subtitle = "Processed GEO DEG summary table",
    caption = "Cutoff: FDR < 0.05 and |log2FC| ≥ 1",
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.4,
    labSize = 3.5,
    colAlpha = 0.85,
    legendPosition = "bottom",
    drawConnectors = TRUE,
    widthConnectors = 0.5
  )
)

dev.off()

# ------------------------------------------------------------------------------
# 9.9 Top upregulated and downregulated genes barplot
# ------------------------------------------------------------------------------

top_up <- deg %>%
  dplyr::filter(Direction == "Up in ZIKV") %>%
  dplyr::arrange(desc(log2FoldChange)) %>%
  dplyr::slice_head(n = 15)

top_down <- deg %>%
  dplyr::filter(Direction == "Down in ZIKV") %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::slice_head(n = 15)

top_bar_df <- dplyr::bind_rows(top_down, top_up) %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = Gene[order(log2FoldChange)])
  )

write.csv(
  top_bar_df,
  "Results/Tables/GSE78711_Top_Up_Down_Genes_For_Barplot.csv",
  row.names = FALSE
)

p_top_bar <- ggplot(
  top_bar_df,
  aes(
    x = Gene,
    y = log2FoldChange,
    fill = Direction
  )
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = deg_colors) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  labs(
    title = "Top differentially expressed genes",
    subtitle = "Top 15 upregulated and top 15 downregulated genes",
    x = NULL,
    y = "log2 fold change",
    fill = "Regulation"
  ) +
  pub_theme

save_pub_plot(
  p_top_bar,
  "Results/Plots/09_Publication_Top_Up_Down_DEG_Barplot",
  width = 8,
  height = 9
)

# ------------------------------------------------------------------------------
# 9.10 Top 30 DEG lollipop plot
# ------------------------------------------------------------------------------

top_lollipop_df <- deg %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = rev(Gene))
  )

write.csv(
  top_lollipop_df,
  "Results/Tables/GSE78711_Top30_DEG_For_Lollipop.csv",
  row.names = FALSE
)

p_lollipop <- ggplot(
  top_lollipop_df,
  aes(
    x = minusLog10FDR,
    y = Gene,
    color = Direction
  )
) +
  geom_segment(
    aes(x = 0, xend = minusLog10FDR, y = Gene, yend = Gene),
    linewidth = 0.7,
    color = "grey60"
  ) +
  geom_point(size = 3.2) +
  scale_color_manual(values = deg_colors) +
  labs(
    title = "Top 30 genes ranked by FDR",
    subtitle = "Higher value indicates stronger statistical significance",
    x = "-log10 FDR",
    y = NULL,
    color = "Regulation"
  ) +
  pub_theme

save_pub_plot(
  p_lollipop,
  "Results/Plots/10_Publication_Top30_DEG_Lollipop",
  width = 8,
  height = 9
)

# ------------------------------------------------------------------------------
# 9.11 Top 50 DEG summary heatmap
# ------------------------------------------------------------------------------

top_heatmap_genes <- deg %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::pull(Gene)

heat_df <- deg %>%
  dplyr::filter(Gene %in% top_heatmap_genes) %>%
  dplyr::select(Gene, Mock_expression, ZIKV_expression, log2FoldChange, padj)

heat_mat <- as.matrix(heat_df[, c("Mock_expression", "ZIKV_expression")])
rownames(heat_mat) <- heat_df$Gene

heat_mat <- log2(heat_mat + 1)

heat_mat_z <- t(scale(t(heat_mat)))
heat_mat_z[is.na(heat_mat_z)] <- 0

annotation_row <- data.frame(
  Regulation = deg$Direction[match(rownames(heat_mat_z), deg$Gene)]
)

rownames(annotation_row) <- rownames(heat_mat_z)

ann_colors <- list(
  Regulation = c(
    "Up in ZIKV" = "#D73027",
    "Down in ZIKV" = "#4575B4",
    "Not significant" = "grey75"
  )
)

make_safe_pdf(
  "Results/Plots/11_Publication_Top50_DEG_Summary_Heatmap.pdf",
  {
    pheatmap(
      heat_mat_z,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_row = annotation_row,
      annotation_colors = ann_colors,
      color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
      border_color = NA,
      fontsize_row = 7,
      fontsize_col = 12,
      main = "Top 50 genes: Mock vs ZIKV summary expression"
    )
  },
  width = 7,
  height = 11
)

png(
  "Results/Plots/11_Publication_Top50_DEG_Summary_Heatmap.png",
  width = 7,
  height = 11,
  units = "in",
  res = 600
)

pheatmap(
  heat_mat_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 12,
  main = "Top 50 genes: Mock vs ZIKV summary expression"
)

dev.off()

# ------------------------------------------------------------------------------
# 9.12 Ranked gene statistic plot
# ------------------------------------------------------------------------------

rank_df <- deg %>%
  dplyr::arrange(desc(log2FoldChange)) %>%
  dplyr::mutate(Rank = dplyr::row_number())

p_rank <- ggplot(
  rank_df,
  aes(
    x = Rank,
    y = log2FoldChange,
    color = Direction
  )
) +
  geom_point(alpha = 0.75, size = 1.4) +
  geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed")) +
  scale_color_manual(values = deg_colors) +
  labs(
    title = "Ranked log2 fold-change profile",
    subtitle = "Genes ordered from highest to lowest log2FC",
    x = "Gene rank",
    y = "log2 fold change",
    color = "Regulation"
  ) +
  pub_theme

save_pub_plot(
  p_rank,
  "Results/Plots/12_Publication_Ranked_log2FC_Profile",
  width = 8,
  height = 6
)

cat("✓ Publication-ready DEG plots completed.\n\n")


# ==============================================================================
# 10. PUBLICATION-READY FUNCTIONAL ENRICHMENT
# ==============================================================================

cat("[9] Running GO/KEGG enrichment using significant DEGs...\n")

main_deg <- deg %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) >= 1)

up_deg <- main_deg %>%
  dplyr::filter(log2FoldChange >= 1)

down_deg <- main_deg %>%
  dplyr::filter(log2FoldChange <= -1)

write.csv(
  main_deg,
  "Results/Tables/GSE78711_Main_DEG_Used_For_Enrichment.csv",
  row.names = FALSE
)

run_enrichment <- function(gene_vector, prefix_name) {
  
  gene_vector <- unique(gene_vector)
  gene_vector <- gene_vector[!is.na(gene_vector)]
  
  if (length(gene_vector) < 5) {
    cat("Skipping enrichment for", prefix_name, ": fewer than 5 genes.\n")
    return(NULL)
  }
  
  if (mean(grepl("^ENSG", gene_vector)) > 0.5) {
    from_type <- "ENSEMBL"
    gene_clean <- sub("\\..*$", "", gene_vector)
  } else {
    from_type <- "SYMBOL"
    gene_clean <- gene_vector
  }
  
  cat("Running enrichment for:", prefix_name, "\n")
  cat("Detected gene ID type:", from_type, "\n")
  
  entrez_conversion <- tryCatch(
    {
      clusterProfiler::bitr(
        gene_clean,
        fromType = from_type,
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      )
    },
    error = function(e) {
      cat("Gene conversion failed for", prefix_name, ":", e$message, "\n")
      data.frame()
    }
  )
  
  write.csv(
    entrez_conversion,
    paste0("Results/Tables/", prefix_name, "_GeneID_to_ENTREZID.csv"),
    row.names = FALSE
  )
  
  if (nrow(entrez_conversion) == 0) {
    cat("No ENTREZ IDs found for", prefix_name, "\n")
    return(NULL)
  }
  
  entrez_ids <- unique(entrez_conversion$ENTREZID)
  
  ego_bp <- tryCatch(
    {
      enrichGO(
        gene = entrez_ids,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
    },
    error = function(e) {
      cat("GO BP failed for", prefix_name, ":", e$message, "\n")
      NULL
    }
  )
  
  ego_mf <- tryCatch(
    {
      enrichGO(
        gene = entrez_ids,
        OrgDb = org.Hs.eg.db,
        ont = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
    },
    error = function(e) {
      cat("GO MF failed for", prefix_name, ":", e$message, "\n")
      NULL
    }
  )
  
  ego_cc <- tryCatch(
    {
      enrichGO(
        gene = entrez_ids,
        OrgDb = org.Hs.eg.db,
        ont = "CC",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
    },
    error = function(e) {
      cat("GO CC failed for", prefix_name, ":", e$message, "\n")
      NULL
    }
  )
  
  ekegg <- tryCatch(
    {
      enrichKEGG(
        gene = entrez_ids,
        organism = "hsa",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
    },
    error = function(e) {
      cat("KEGG failed for", prefix_name, ":", e$message, "\n")
      NULL
    }
  )
  
  enrichment_list <- list(
    GO_BP = ego_bp,
    GO_MF = ego_mf,
    GO_CC = ego_cc,
    KEGG = ekegg
  )
  
  for (nm in names(enrichment_list)) {
    
    obj <- enrichment_list[[nm]]
    
    if (is.null(obj)) {
      next
    }
    
    res_df <- as.data.frame(obj)
    
    write.csv(
      res_df,
      paste0("Results/Tables/", prefix_name, "_", nm, "_Enrichment.csv"),
      row.names = FALSE
    )
    
    if (nrow(res_df) > 0) {
      
      make_safe_pdf(
        paste0("Results/Plots/", prefix_name, "_", nm, "_Dotplot.pdf"),
        {
          print(
            dotplot(
              obj,
              showCategory = 15,
              title = paste0(prefix_name, " ", nm, " enrichment")
            ) +
              pub_theme
          )
        },
        width = 11,
        height = 7
      )
      
      png(
        paste0("Results/Plots/", prefix_name, "_", nm, "_Dotplot.png"),
        width = 11,
        height = 7,
        units = "in",
        res = 600
      )
      
      print(
        dotplot(
          obj,
          showCategory = 15,
          title = paste0(prefix_name, " ", nm, " enrichment")
        ) +
          pub_theme
      )
      
      dev.off()
      
      make_safe_pdf(
        paste0("Results/Plots/", prefix_name, "_", nm, "_Barplot.pdf"),
        {
          print(
            barplot(
              obj,
              showCategory = 15,
              title = paste0(prefix_name, " ", nm, " enrichment")
            )
          )
        },
        width = 11,
        height = 7
      )
      
      png(
        paste0("Results/Plots/", prefix_name, "_", nm, "_Barplot.png"),
        width = 11,
        height = 7,
        units = "in",
        res = 600
      )
      
      print(
        barplot(
          obj,
          showCategory = 15,
          title = paste0(prefix_name, " ", nm, " enrichment")
        )
      )
      
      dev.off()
    }
  }
  
  return(enrichment_list)
}

enrich_all <- run_enrichment(
  gene_vector = main_deg$Gene,
  prefix_name = "GSE78711_All_Significant_DEG"
)

enrich_up <- run_enrichment(
  gene_vector = up_deg$Gene,
  prefix_name = "GSE78711_Up_in_ZIKV_DEG"
)

enrich_down <- run_enrichment(
  gene_vector = down_deg$Gene,
  prefix_name = "GSE78711_Down_in_ZIKV_DEG"
)

cat("✓ Enrichment analysis completed.\n\n")


# ==============================================================================
# 11. SAVE SESSION INFO
# ==============================================================================

sink("Results/Annotation/sessionInfo.txt")
print(sessionInfo())
sink()


# ==============================================================================
# 12. FINISH
# ==============================================================================

cat("\n============================================================\n")
cat("WORKFLOW COMPLETED SUCCESSFULLY\n")
cat("Important note:\n")
cat("  GEO total samples = 8.\n")
cat("  These 8 samples are split across 2 GEO matrix/platform objects.\n")
cat("  This script combines both matrix objects for phenotype annotation.\n")
cat("  GSE78711_gene.exp.all.txt.gz is a processed DEG summary table.\n")
cat("  Therefore, DESeq2/limma sample-level analysis is NOT run here.\n")
cat("  PCA/sample clustering/sample correlation are NOT performed because no valid\n")
cat("  sample-level expression matrix is present in this processed DEG table.\n")
cat("  For true 8-sample DESeq2 analysis, download and reprocess raw SRA reads.\n")
cat("\nMain output folders:\n")
cat("  Results/QC\n")
cat("  Results/Plots\n")
cat("  Results/Tables\n")
cat("  Results/Annotation\n")
cat("\nMain publication-ready plots:\n")
cat("  Results/Plots/01_Publication_DEG_Count_Barplot.pdf/png\n")
cat("  Results/Plots/05_Publication_Mock_vs_ZIKV_Expression_Scatter.pdf/png\n")
cat("  Results/Plots/06_Publication_MA_Plot_ZIKV_vs_Mock.pdf/png\n")
cat("  Results/Plots/07_Publication_Volcano_ZIKV_vs_Mock.pdf/png\n")
cat("  Results/Plots/08_Publication_EnhancedVolcano_ZIKV_vs_Mock.pdf/png\n")
cat("  Results/Plots/09_Publication_Top_Up_Down_DEG_Barplot.pdf/png\n")
cat("  Results/Plots/10_Publication_Top30_DEG_Lollipop.pdf/png\n")
cat("  Results/Plots/11_Publication_Top50_DEG_Summary_Heatmap.pdf/png\n")
cat("  Results/Plots/12_Publication_Ranked_log2FC_Profile.pdf/png\n")
cat("============================================================\n")

