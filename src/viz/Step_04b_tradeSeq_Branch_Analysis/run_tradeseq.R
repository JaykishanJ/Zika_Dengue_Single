library(tradeSeq)
library(SingleCellExperiment)

setwd("E:/Zika_Enrichment/Publication_Pipeline")

print("Loading data...")
counts <- read.csv("Step_04_Fig_4_Pseudotime/tradeseq_counts.csv", row.names=1)
meta <- read.csv("Step_04_Fig_4_Pseudotime/tradeseq_metadata.csv", row.names=1)

counts_mat <- as.matrix(counts)

# We have 2 lineages: DENV and ZIKV
cell_weights <- matrix(0, nrow=nrow(meta), ncol=2)
rownames(cell_weights) <- rownames(meta)
colnames(cell_weights) <- c("DENV_Lineage", "ZIKV_Lineage")

# Control cells belong to both lineages
cell_weights[meta$condition == "Control", 1] <- 1
cell_weights[meta$condition == "Control", 2] <- 1

# DENV cells belong to Lineage 1
cell_weights[meta$condition == "DENV", 1] <- 1

# ZIKV cells belong to Lineage 2
cell_weights[meta$condition == "ZIKV", 2] <- 1

pseudotime_mat <- matrix(meta$pseudotime, nrow=nrow(meta), ncol=2)
rownames(pseudotime_mat) <- rownames(meta)
colnames(pseudotime_mat) <- c("DENV_Lineage", "ZIKV_Lineage")

print("Filtering genes and cells...")
# Filter cells with 0 counts
cell_sums <- colSums(counts_mat)
counts_mat <- counts_mat[, cell_sums > 0]
meta <- meta[cell_sums > 0, ]
cell_weights <- cell_weights[cell_sums > 0, ]
pseudotime_mat <- pseudotime_mat[cell_sums > 0, ]

# Filter genes: only keep the top 2000 most highly variable genes to speed up computation
# Calculate variance of each gene across all cells
gene_vars <- apply(counts_mat, 1, var)
# Get the indices of the top 2000 genes by variance
top_genes <- order(gene_vars, decreasing=TRUE)[1:2000]
counts_mat <- counts_mat[top_genes, ]
# Also ensure these top 2000 genes have at least some expression across 10 cells
keep <- rowSums(counts_mat > 0) >= 10
counts_mat <- counts_mat[keep, ]

print("Fitting GAM... (This may take a few minutes)")
set.seed(42)
sce <- fitGAM(counts = counts_mat,
              pseudotime = pseudotime_mat,
              cellWeights = cell_weights,
              nknots = 6, 
              verbose = TRUE)

print("Running patternTest to identify branch-specific genes...")
patternRes <- patternTest(sce)

patternRes$FDR <- p.adjust(patternRes$pvalue, method="fdr")
sig_genes <- sum(patternRes$FDR < 0.05, na.rm=TRUE)

print(paste(">>> EXACT BRANCH-SPECIFIC GENE COUNT:", sig_genes, "<<<"))

write.csv(patternRes, "Step_04_Fig_4_Pseudotime/tradeSeq_branch_specific_genes.csv")
print("Results exported to tradeSeq_branch_specific_genes.csv")
