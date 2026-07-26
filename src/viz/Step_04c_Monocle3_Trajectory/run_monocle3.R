# run_monocle3.R
# Script to run Monocle3 trajectory analysis on the extracted Zika/DENV dataset

print("Loading libraries...")
suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
})

# Handle working directory seamlessly
# Hardcoding to the absolute path of the pipeline to avoid relative path errors in R
setwd("e:/Zika_Enrichment/Publication_Pipeline")
base_dir <- "."
out_dir <- file.path(base_dir, "Step_04c_Monocle3_Trajectory")

# Ensure the output directory exists
if (!dir.exists(out_dir)) dir.create(out_dir, showWarnings = FALSE)

print("Loading data...")
# Load the counts matrix (genes x cells)
counts_file <- file.path(base_dir, "Step_04b_tradeSeq_Branch_Analysis/tradeseq_counts.csv")
meta_file <- file.path(base_dir, "Step_04b_tradeSeq_Branch_Analysis/tradeseq_metadata.csv")

if (!file.exists(counts_file) || !file.exists(meta_file)) {
  stop("Could not find the input CSV files. Please check your working directory.")
}

counts_df <- read.csv(counts_file, row.names=1)
counts_mat <- as.matrix(counts_df)
counts_sparse <- Matrix(counts_mat, sparse=TRUE)

# Load cell metadata
meta <- read.csv(meta_file, row.names=1)
# Ensure columns align
meta <- meta[colnames(counts_sparse), , drop=FALSE]

# Create gene metadata dataframe
gene_meta <- data.frame(gene_short_name = rownames(counts_sparse), row.names = rownames(counts_sparse))

print("Creating cell_data_set object...")
cds <- new_cell_data_set(counts_sparse,
                         cell_metadata = meta,
                         gene_metadata = gene_meta)

print("Preprocessing data (PCA)...")
cds <- preprocess_cds(cds, num_dim = 50)

print("Reducing dimensions (UMAP)...")
cds <- reduce_dimension(cds, reduction_method="UMAP")

print("Clustering cells...")
cds <- cluster_cells(cds)

print("Learning trajectory graph...")
cds <- learn_graph(cds)

print("Identifying root cells for pseudotime ordering...")
# Find the cells belonging to the Control / Bystander populations to set as the root
root_cells <- rownames(meta)[meta$condition == "Control" | meta$condition == "Mock"]
if (length(root_cells) == 0) {
  root_cells <- rownames(meta)[grepl("Bystander", meta$infection_state)]
}

print(paste("Found", length(root_cells), "root cells. Ordering cells..."))
cds <- order_cells(cds, root_cells = root_cells)

# Define figure output directory
fig_dir <- file.path(base_dir, "final_publication_figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, showWarnings = FALSE)

print("Plotting trajectories...")
# Panel A: Pseudotime
p1 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           show_trajectory_graph = FALSE,
           cell_size = 0.75) +
  theme(legend.position = "right") +
  ggtitle("Pseudotime Trajectory")

ggsave(file.path(fig_dir, "Monocle3_Pseudotime.pdf"), plot = p1, width = 6, height = 5)

# Panel B: Condition
p2 <- plot_cells(cds,
           color_cells_by = "condition",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           show_trajectory_graph = FALSE,
           cell_size = 0.75) +
  theme(legend.position = "right") +
  ggtitle("Condition Trajectory")

ggsave(file.path(fig_dir, "Monocle3_Condition.pdf"), plot = p2, width = 6, height = 5)

# Panel C: Infection States
p3 <- plot_cells(cds,
           color_cells_by = "infection_state",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           show_trajectory_graph = FALSE,
           cell_size = 0.75) +
  theme(legend.position = "right") +
  ggtitle("Infection States")

ggsave(file.path(fig_dir, "Monocle3_InfectionState.pdf"), plot = p3, width = 6, height = 5)

print("Exporting pseudotime results...")
pt_df <- data.frame(cell_id = rownames(meta), monocle3_pseudotime = pseudotime(cds))
write.csv(pt_df, file.path(out_dir, "monocle3_pseudotime_results.csv"), row.names=FALSE)

print("Monocle3 Analysis Complete!")
