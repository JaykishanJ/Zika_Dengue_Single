################################################################################
# FINAL CLEAN SCRIPT
# RNA-SEQ DIFFERENTIAL EXPRESSION ANALYSIS
# Dataset: GSE118305
#
# Study:
# Genome-wide RNA-seq and epigenetic analysis of human monocyte-derived
# macrophage responses to Zika virus infection
#
# Important:
# GEO contains 60 total samples, but these include RNA-seq + ChIP-seq samples.
# This script analyzes only the RNA-seq expression/count matrix:
# GSE118305_RNAseq_HMDM_ZIKV_2018.txt.gz
#
# RNA-seq groups are extracted from matrix column names:
# Mock
# ZIKV_negative / bystander / 4G2neg
# ZIKV_positive / infected / 4G2pos
################################################################################

rm(list = ls())
gc()

# ==============================================================================
# 0. USER SETTINGS
# ==============================================================================

WORKDIR <- "E:/zika_virus/GSE118305"
GSE_ID  <- "GSE118305"

setwd(WORKDIR)

dir.create("Results", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Plots", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/QC", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Tables", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Annotation", showWarnings = FALSE, recursive = TRUE)

cat("\n============================================================\n")
cat("RNA-seq DESeq2 Analysis:", GSE_ID, "\n")
cat("Working directory:", getwd(), "\n")
cat("============================================================\n\n")


# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================

cat("[1] Loading required libraries...\n")

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ashr)
  library(AnnotationDbi)
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
  x <- gsub("\\.", "-", x)
  x <- gsub("_$", "", x)
  x
}


make_safe_pdf <- function(filename, expr, width = 8, height = 6) {
  pdf(filename, width = width, height = height)
  tryCatch(
    expr,
    error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Plot failed:\n", e$message), cex = 1)
      message("Plot failed: ", filename, " | ", e$message)
    },
    finally = dev.off()
  )
}


# ==============================================================================
# 3. DOWNLOAD GEO METADATA
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
  "Results/Annotation/All_GEO_Metadata_60_Total_Samples.csv",
  row.names = FALSE
)

cat("Total GEO samples:", nrow(metadata), "\n")
cat("Note: These include RNA-seq and ChIP-seq/epigenetic samples.\n\n")


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

extract_dir <- file.path(getwd(), paste0(GSE_ID, "_extracted"))
dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)

tar_files <- all_files[grepl("\\.tar$", all_files, ignore.case = TRUE)]

if (length(tar_files) > 0) {
  for (tf in tar_files) {
    cat("Extracting TAR:", tf, "\n")
    untar(tf, exdir = extract_dir)
  }
}

zip_files <- all_files[grepl("\\.zip$", all_files, ignore.case = TRUE)]

if (length(zip_files) > 0) {
  for (zf in zip_files) {
    cat("Extracting ZIP:", zf, "\n")
    unzip(zf, exdir = extract_dir)
  }
}

all_files <- unique(c(
  all_files,
  list.files(extract_dir, recursive = TRUE, full.names = TRUE)
))

write.csv(
  data.frame(
    File = all_files,
    Basename = basename(all_files)
  ),
  "Results/Annotation/All_Supplementary_Files_Found.csv",
  row.names = FALSE
)

cat("\nAll files found after download/extraction:\n")
print(basename(all_files))
cat("\n")


# ==============================================================================
# 5. FIND RNA-SEQ MATRIX FILE
# ==============================================================================

cat("[4] Finding RNA-seq matrix file...\n")

rnaseq_candidates <- all_files[
  grepl("RNAseq|rnaseq|RNA_seq|rna_seq", basename(all_files), ignore.case = TRUE) &
    grepl("\\.(txt|tsv|csv|gz)$", basename(all_files), ignore.case = TRUE) &
    !grepl("series_matrix|family|soft|chip|bed|bigwig|bw|narrowPeak|broadPeak|peaks|regions", basename(all_files), ignore.case = TRUE)
]

cat("Candidate RNA-seq files:\n")
print(rnaseq_candidates)
cat("\n")

if (length(rnaseq_candidates) == 0) {
  stop(
    "No RNA-seq matrix file found. Check Results/Annotation/All_Supplementary_Files_Found.csv"
  )
}

# Prefer the known GSE118305 RNA-seq matrix file
preferred_file <- rnaseq_candidates[
  grepl("GSE118305_RNAseq_HMDM_ZIKV_2018", basename(rnaseq_candidates), ignore.case = TRUE)
]

if (length(preferred_file) >= 1) {
  rnaseq_file <- preferred_file[1]
} else {
  rnaseq_file <- rnaseq_candidates[1]
}

cat("Selected RNA-seq matrix file:\n")
cat(rnaseq_file, "\n\n")


# ==============================================================================
# 6. LOAD RNA-SEQ MATRIX
# ==============================================================================

cat("[5] Loading RNA-seq matrix...\n")

raw_df <- safe_read_table(rnaseq_file)

cat("Raw RNA-seq matrix dimensions:\n")
cat("Rows:", nrow(raw_df), "\n")
cat("Columns:", ncol(raw_df), "\n\n")

cat("First 15 column names:\n")
print(head(colnames(raw_df), 15))
cat("\n")

write.csv(
  data.frame(Column_Names = colnames(raw_df)),
  "Results/Annotation/Selected_RNAseq_Matrix_Column_Names.csv",
  row.names = FALSE
)


# ==============================================================================
# 7. FORMAT HOMER ANALYZEREPEATS MATRIX
# ==============================================================================

cat("[6] Formatting HOMER analyzeRepeats RNA-seq matrix...\n")

counts_df <- raw_df
colnames(counts_df) <- clean_sample_names(colnames(counts_df))

# This matrix is HOMER analyzeRepeats-style:
# Column 1  = Transcript/RepeatID
# Column 2  = chr
# Column 3  = start
# Column 4  = end
# Column 5  = strand
# Column 6  = Length
# Column 7  = gene_name
# Column 8  = Annotation/Divergence
# Column 9  = gene_id
# Column 10 = gene_type
# Column 11 onward = RNA-seq expression/count columns

if (ncol(counts_df) < 11) {
  stop("Matrix has fewer than 11 columns. This does not look like expected HOMER output.")
}

annotation_cols <- 1:10
expr_cols <- 11:ncol(counts_df)

gene_annotation <- counts_df[, annotation_cols, drop = FALSE]
expr_df <- counts_df[, expr_cols, drop = FALSE]

# Save annotation
write.csv(
  gene_annotation,
  "Results/Annotation/Gene_Annotation_From_RNAseq_Matrix.csv",
  row.names = FALSE
)

# Use gene_name as primary gene ID if available
if ("gene_name" %in% colnames(gene_annotation)) {
  gene_ids <- as.character(gene_annotation$gene_name)
} else {
  gene_ids <- as.character(gene_annotation[[1]])
}

fallback_gene_ids <- as.character(gene_annotation[[1]])

gene_ids[is.na(gene_ids) | gene_ids == ""] <- fallback_gene_ids[is.na(gene_ids) | gene_ids == ""]
gene_ids <- make.unique(gene_ids)

rownames(expr_df) <- gene_ids

cat("Expression columns before sample filtering:\n")
print(colnames(expr_df))
cat("\n")

# Remove pooled columns because pooled samples are not independent biological replicates
expr_df <- expr_df[, !grepl("pooled", colnames(expr_df), ignore.case = TRUE), drop = FALSE]

cat("Expression columns after removing pooled samples:\n")
print(colnames(expr_df))
cat("\n")


# ==============================================================================
# 8. BUILD RNA-SEQ PHENOTYPE FROM MATRIX COLUMN NAMES
# ==============================================================================

cat("[7] Building RNA-seq phenotype from matrix column names...\n")

sample_info <- data.frame(
  SampleID = colnames(expr_df),
  stringsAsFactors = FALSE
)

sample_info$Condition <- dplyr::case_when(
  grepl("mock", sample_info$SampleID, ignore.case = TRUE) ~ "Mock",
  grepl("4G2neg|ZIKVneg|ZIKVnegative|negative", sample_info$SampleID, ignore.case = TRUE) ~ "ZIKV_negative",
  grepl("4G2pos|ZIKVpos|ZIKVpositive|positive", sample_info$SampleID, ignore.case = TRUE) ~ "ZIKV_positive",
  TRUE ~ NA_character_
)

sample_info$Timepoint <- dplyr::case_when(
  grepl("12h", sample_info$SampleID, ignore.case = TRUE) ~ "12h",
  grepl("18h", sample_info$SampleID, ignore.case = TRUE) ~ "18h",
  grepl("24h", sample_info$SampleID, ignore.case = TRUE) ~ "24h",
  TRUE ~ NA_character_
)

sample_info$Patient <- stringr::str_extract(sample_info$SampleID, "p[0-9]+")
sample_info$Patient <- ifelse(is.na(sample_info$Patient), sample_info$SampleID, sample_info$Patient)

sample_info <- sample_info %>%
  dplyr::filter(!is.na(Condition))

rownames(sample_info) <- sample_info$SampleID

sample_info$Condition <- factor(
  sample_info$Condition,
  levels = c("Mock", "ZIKV_negative", "ZIKV_positive")
)

sample_info$Timepoint <- factor(
  sample_info$Timepoint,
  levels = c("12h", "18h", "24h")
)

write.csv(
  sample_info,
  "Results/Annotation/RNAseq_Sample_Design_From_Matrix_Columns.csv",
  row.names = FALSE
)

cat("RNA-seq phenotype summary from matrix columns:\n")
print(table(sample_info$Condition))
cat("\n")

cat("RNA-seq timepoint summary:\n")
print(table(sample_info$Timepoint, sample_info$Condition))
cat("\n")

if (nrow(sample_info) < 3) {
  stop("Too few RNA-seq samples were detected from matrix column names.")
}


# ==============================================================================
# 9. CREATE COUNT MATRIX
# ==============================================================================

cat("[8] Creating DESeq2 count matrix...\n")

expr_df <- expr_df[, rownames(sample_info), drop = FALSE]

counts_data <- as.data.frame(expr_df, check.names = FALSE)

for (j in seq_len(ncol(counts_data))) {
  counts_data[[j]] <- suppressWarnings(as.numeric(as.character(counts_data[[j]])))
}

if (any(is.na(counts_data))) {
  bad_cols <- colnames(counts_data)[colSums(is.na(counts_data)) > 0]
  stop(
    "NA values created during numeric conversion. Problem columns: ",
    paste(bad_cols, collapse = ", ")
  )
}

counts_data <- as.matrix(counts_data)
counts_data <- round(counts_data)
storage.mode(counts_data) <- "integer"

# Remove bad/empty genes
valid_genes <- !is.na(rownames(counts_data)) & rownames(counts_data) != ""
counts_data <- counts_data[valid_genes, , drop = FALSE]

counts_data <- counts_data[rowSums(counts_data) > 0, , drop = FALSE]

# Remove duplicated rownames if any remain
if (any(duplicated(rownames(counts_data)))) {
  rownames(counts_data) <- make.unique(rownames(counts_data))
}

rna_metadata <- sample_info[colnames(counts_data), , drop = FALSE]

if (!all(colnames(counts_data) == rownames(rna_metadata))) {
  stop("Sample order mismatch between count matrix and RNA metadata.")
}

write.csv(
  as.data.frame(counts_data[1:min(20, nrow(counts_data)), , drop = FALSE]),
  "Results/Annotation/Count_Matrix_First20Genes_Check.csv"
)

write.csv(
  data.frame(
    SampleID = rownames(rna_metadata),
    Condition = rna_metadata$Condition,
    Timepoint = rna_metadata$Timepoint,
    Patient = rna_metadata$Patient
  ),
  "Results/Annotation/Final_RNAseq_Sample_Design.csv",
  row.names = FALSE
)

cat("✓ RNA-seq count matrix loaded successfully.\n")
cat("Genes:", nrow(counts_data), "\n")
cat("Samples:", ncol(counts_data), "\n")
cat("Final phenotype summary:\n")
print(table(rna_metadata$Condition))
cat("\n")


# ==============================================================================
# 10. MAIN RNA-SEQ ANALYSIS SUBSET
# ==============================================================================

cat("[9] Preparing main analysis subset...\n")

# For a clean biological contrast, use 24h samples.
# Mock samples are 24h. ZIKV_negative and ZIKV_positive also have 24h samples.
# This avoids mixing 12h, 18h, and 24h in one simple design.

main_samples <- rownames(rna_metadata)[
  rna_metadata$Timepoint == "24h" |
    rna_metadata$Condition == "Mock"
]

counts_main <- counts_data[, main_samples, drop = FALSE]
meta_main <- rna_metadata[main_samples, , drop = FALSE]

# Drop unused levels
meta_main$Condition <- droplevels(meta_main$Condition)
meta_main$Timepoint <- droplevels(meta_main$Timepoint)

cat("Main analysis sample summary:\n")
print(table(meta_main$Condition))
cat("\n")

write.csv(
  meta_main,
  "Results/Annotation/Main_Analysis_24h_Sample_Design.csv",
  row.names = FALSE
)

if (any(table(meta_main$Condition) < 2)) {
  warning("One or more groups have fewer than 2 samples. Interpret DESeq2 results carefully.")
}


# ==============================================================================
# 11. DESEQ2 ANALYSIS
# ==============================================================================

cat("[10] Running DESeq2...\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_main,
  colData = meta_main,
  design = ~ Condition
)

keep <- rowSums(counts(dds)) >= 15
dds <- dds[keep, ]

cat("Genes after low-count filtering:", nrow(dds), "\n")

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)

cat("✓ DESeq2 completed.\n\n")


# ==============================================================================
# 12. QC PLOTS
# ==============================================================================

cat("[11] Generating QC plots...\n")

sample_colors <- c(
  "Mock" = "#4DAF4A",
  "ZIKV_negative" = "#FBBC05",
  "ZIKV_positive" = "#E41A1C"
)

make_safe_pdf(
  "Results/QC/01_Dispersion_Profiles.pdf",
  {
    plotDispEsts(dds, main = "DESeq2 Dispersion Estimates")
  },
  width = 7,
  height = 6
)

make_safe_pdf(
  "Results/QC/02_Normalization_Diagnostics.pdf",
  {
    par(mfrow = c(1, 2))
    
    boxplot(
      log2(counts(dds, normalized = FALSE) + 1),
      outline = FALSE,
      las = 2,
      main = "Raw Counts",
      col = "#EA4335",
      ylab = "log2(count + 1)"
    )
    
    boxplot(
      assay(vsd),
      outline = FALSE,
      las = 2,
      main = "VST Normalized",
      col = "#4285F4",
      ylab = "VST expression"
    )
    
    par(mfrow = c(1, 1))
  },
  width = 12,
  height = 6
)

pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_color_manual(values = sample_colors) +
  labs(
    title = paste0("PCA: ", GSE_ID, " RNA-seq 24h"),
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_bw(base_size = 14)

group_counts <- table(pca_data$Condition)

if (all(group_counts >= 3)) {
  p_pca <- p_pca + stat_ellipse(level = 0.95, linewidth = 1, linetype = 2)
}

ggsave(
  "Results/QC/03_PCA_Analysis.pdf",
  p_pca,
  width = 8,
  height = 6
)

sample_cor <- cor(assay(vsd))

annotation_col <- data.frame(
  Condition = colData(vsd)$Condition
)

rownames(annotation_col) <- colnames(vsd)

make_safe_pdf(
  "Results/QC/04_Sample_Correlation_Heatmap.pdf",
  {
    pheatmap(
      sample_cor,
      annotation_col = annotation_col,
      annotation_colors = list(Condition = sample_colors),
      border_color = NA,
      main = "Sample Correlation Heatmap"
    )
  },
  width = 10,
  height = 10
)

cat("✓ QC plots completed.\n\n")


# ==============================================================================
# 13. DIFFERENTIAL EXPRESSION CONTRASTS
# ==============================================================================

cat("[12] Running differential expression contrasts...\n")

contrast_list <- list(
  ZIKV_positive_vs_ZIKV_negative = c("Condition", "ZIKV_positive", "ZIKV_negative"),
  ZIKV_positive_vs_Mock          = c("Condition", "ZIKV_positive", "Mock"),
  ZIKV_negative_vs_Mock          = c("Condition", "ZIKV_negative", "Mock")
)

all_contrast_results <- list()
summary_table <- data.frame()

for (contrast_name in names(contrast_list)) {
  
  cat("\nContrast:", contrast_name, "\n")
  
  contrast_vec <- contrast_list[[contrast_name]]
  
  res_raw <- results(dds, contrast = contrast_vec)
  
  res_shrunk <- lfcShrink(
    dds,
    contrast = contrast_vec,
    res = res_raw,
    type = "ashr"
  )
  
  deg_all <- as.data.frame(res_shrunk)
  deg_all$Gene <- rownames(deg_all)
  deg_all$Contrast <- contrast_name
  
  deg_all <- deg_all %>%
    dplyr::select(
      Gene,
      Contrast,
      baseMean,
      log2FoldChange,
      lfcSE,
      pvalue,
      padj
    ) %>%
    dplyr::arrange(padj)
  
  sig_deg <- deg_all %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)
  
  up_deg <- sig_deg %>%
    dplyr::filter(log2FoldChange > 1)
  
  down_deg <- sig_deg %>%
    dplyr::filter(log2FoldChange < -1)
  
  write.csv(
    deg_all,
    paste0("Results/Tables/", contrast_name, "_All_Genes.csv"),
    row.names = FALSE
  )
  
  write.csv(
    sig_deg,
    paste0("Results/Tables/", contrast_name, "_Significant_DEGs_padj0.05_log2FC1.csv"),
    row.names = FALSE
  )
  
  write.csv(
    up_deg,
    paste0("Results/Tables/", contrast_name, "_Upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    down_deg,
    paste0("Results/Tables/", contrast_name, "_Downregulated.csv"),
    row.names = FALSE
  )
  
  cat("Total significant DEGs:", nrow(sig_deg), "\n")
  cat("Upregulated:", nrow(up_deg), "\n")
  cat("Downregulated:", nrow(down_deg), "\n")
  
  summary_table <- rbind(
    summary_table,
    data.frame(
      Contrast = contrast_name,
      Total_DEG = nrow(sig_deg),
      Upregulated = nrow(up_deg),
      Downregulated = nrow(down_deg)
    )
  )
  
  all_contrast_results[[contrast_name]] <- deg_all
  
  make_safe_pdf(
    paste0("Results/QC/", contrast_name, "_Pvalue_Distribution.pdf"),
    {
      hist(
        deg_all$pvalue,
        breaks = 50,
        main = paste0("P-value Distribution: ", contrast_name),
        xlab = "P-value",
        ylab = "Frequency",
        col = "darkgreen",
        border = "white"
      )
    },
    width = 8,
    height = 6
  )
  
  make_safe_pdf(
    paste0("Results/Plots/", contrast_name, "_Enhanced_Volcano.pdf"),
    {
      print(
        EnhancedVolcano(
          deg_all,
          lab = deg_all$Gene,
          x = "log2FoldChange",
          y = "padj",
          title = contrast_name,
          subtitle = "DESeq2 with ashr log2FC shrinkage",
          pCutoff = 0.05,
          FCcutoff = 1,
          pointSize = 2.5,
          labSize = 3.5,
          colAlpha = 0.85,
          legendPosition = "bottom"
        )
      )
    },
    width = 10,
    height = 8
  )
  
  top_genes <- sig_deg %>%
    dplyr::arrange(padj) %>%
    head(50) %>%
    dplyr::pull(Gene)
  
  if (length(top_genes) >= 2) {
    
    heatmap_matrix <- assay(vsd)[top_genes, , drop = FALSE]
    
    make_safe_pdf(
      paste0("Results/Plots/", contrast_name, "_Top50_DEG_Heatmap.pdf"),
      {
        pheatmap(
          heatmap_matrix,
          scale = "row",
          annotation_col = annotation_col,
          annotation_colors = list(Condition = sample_colors),
          clustering_method = "ward.D2",
          show_colnames = FALSE,
          fontsize_row = 7,
          color = colorRampPalette(c("#4285F4", "white", "#EA4335"))(100),
          main = paste0("Top 50 DEGs: ", contrast_name)
        )
      },
      width = 12,
      height = 10
    )
  }
  
  top_6_targets <- sig_deg %>%
    dplyr::arrange(padj) %>%
    head(6) %>%
    dplyr::pull(Gene)
  
  if (length(top_6_targets) >= 1) {
    
    plot_violin_data <- as.data.frame(assay(vsd)[top_6_targets, , drop = FALSE]) %>%
      tibble::rownames_to_column("Gene") %>%
      tidyr::pivot_longer(
        -Gene,
        names_to = "SampleID",
        values_to = "Normalized_Expression"
      ) %>%
      dplyr::left_join(
        data.frame(
          SampleID = rownames(meta_main),
          Condition = meta_main$Condition
        ),
        by = "SampleID"
      )
    
    p_violin <- ggplot(
      plot_violin_data,
      aes(x = Condition, y = Normalized_Expression, fill = Condition)
    ) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.1, size = 1.5, alpha = 0.8) +
      facet_wrap(~ Gene, scales = "free_y") +
      scale_fill_manual(values = sample_colors) +
      labs(
        title = paste0("Top Candidate Gene Expression: ", contrast_name),
        y = "VST normalized expression",
        x = ""
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    ggsave(
      paste0("Results/Plots/", contrast_name, "_Top6_Gene_Violins.pdf"),
      p_violin,
      width = 10,
      height = 7
    )
  }
}

combined_results <- dplyr::bind_rows(all_contrast_results)

write.csv(
  combined_results,
  "Results/Tables/Combined_All_Contrasts_DESeq2_Results.csv",
  row.names = FALSE
)

write.csv(
  summary_table,
  "Results/Tables/DEG_Summary_All_Contrasts.csv",
  row.names = FALSE
)

cat("\n✓ All contrasts completed.\n\n")


# ==============================================================================
# 14. FUNCTIONAL ENRICHMENT FOR MAIN CONTRAST
# ==============================================================================

cat("[13] Running GO/KEGG enrichment for main contrast...\n")

main_contrast <- "ZIKV_positive_vs_ZIKV_negative"

main_deg <- all_contrast_results[[main_contrast]] %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

if (nrow(main_deg) == 0) {
  
  cat("No significant DEGs for enrichment in main contrast.\n")
  
} else {
  
  gene_ids <- main_deg$Gene
  
  if (mean(grepl("^ENSG", gene_ids)) > 0.5) {
    from_type <- "ENSEMBL"
    gene_ids_clean <- sub("\\..*$", "", gene_ids)
  } else {
    from_type <- "SYMBOL"
    gene_ids_clean <- gene_ids
  }
  
  cat("Detected gene ID type:", from_type, "\n")
  
  entrez_conversion <- tryCatch(
    {
      clusterProfiler::bitr(
        gene_ids_clean,
        fromType = from_type,
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      )
    },
    error = function(e) {
      cat("Gene conversion failed:", e$message, "\n")
      data.frame()
    }
  )
  
  write.csv(
    entrez_conversion,
    "Results/Tables/MainContrast_GeneID_to_ENTREZID.csv",
    row.names = FALSE
  )
  
  if (nrow(entrez_conversion) == 0) {
    
    cat("No genes converted to ENTREZID. Skipping enrichment.\n")
    
  } else {
    
    entrez_ids <- unique(entrez_conversion$ENTREZID)
    
    ego <- tryCatch(
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
        cat("GO enrichment failed:", e$message, "\n")
        NULL
      }
    )
    
    if (!is.null(ego)) {
      
      go_result <- as.data.frame(ego)
      
      write.csv(
        go_result,
        "Results/Tables/MainContrast_GO_BP_Enrichment.csv",
        row.names = FALSE
      )
      
      if (nrow(go_result) > 0) {
        make_safe_pdf(
          "Results/Plots/MainContrast_GO_BP_Dotplot.pdf",
          {
            print(dotplot(ego, showCategory = 15, title = "GO BP: ZIKV Positive vs Negative"))
          },
          width = 11,
          height = 7
        )
      }
    }
    
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
        cat("KEGG enrichment failed:", e$message, "\n")
        NULL
      }
    )
    
    if (!is.null(ekegg)) {
      
      kegg_result <- as.data.frame(ekegg)
      
      write.csv(
        kegg_result,
        "Results/Tables/MainContrast_KEGG_Enrichment.csv",
        row.names = FALSE
      )
      
      if (nrow(kegg_result) > 0) {
        
        ekegg_readable <- tryCatch(
          {
            setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
          },
          error = function(e) ekegg
        )
        
        make_safe_pdf(
          "Results/Plots/MainContrast_KEGG_Dotplot.pdf",
          {
            print(dotplot(ekegg_readable, showCategory = 15, title = "KEGG: ZIKV Positive vs Negative"))
          },
          width = 11,
          height = 7
        )
      }
    }
  }
}


# ==============================================================================
# 15. FINISH
# ==============================================================================

cat("\n============================================================\n")
cat("WORKFLOW COMPLETED SUCCESSFULLY\n")
cat("Important note:\n")
cat("  GEO total samples = 60, including RNA-seq and ChIP-seq samples.\n")
cat("  This script analyzes RNA-seq matrix samples only.\n")
cat("Output folders:\n")
cat("  Results/QC\n")
cat("  Results/Plots\n")
cat("  Results/Tables\n")
cat("  Results/Annotation\n")
cat("============================================================\n")

