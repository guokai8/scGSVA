#!/usr/bin/env Rscript
# =============================================================================
# Comprehensive comparison: scGSVA vs competing single-cell gene set scoring
# tools (AUCell, UCell, singscore, GSVA)
#
# Output: results/ folder with comparison tables and figures
# =============================================================================

library(devtools)
load_all("/Users/bioguo/Library/Mobile Documents/com~apple~CloudDocs/Work/Packages/scGSVA")
library(ggplot2)
library(tidyr)
library(dplyr)

results_dir <- "/Users/bioguo/Library/Mobile Documents/com~apple~CloudDocs/Work/Packages/scGSVA/results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Prepare data and gene sets
# =============================================================================
cat("=== Loading data ===\n")
data(pbmcs)
hsko <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")

# Extract expression matrix
expr_mat <- as.matrix(Seurat::GetAssayData(pbmcs, assay = "RNA", layer = "counts"))
cat(sprintf("Expression matrix: %d genes x %d cells\n", nrow(expr_mat), ncol(expr_mat)))

# Build gene set list (named list of character vectors)
annot_df <- as.data.frame(hsko)
gene_sets <- split(annot_df$GeneID, annot_df$Annot)
# Make gene sets unique and filter by size
gene_sets <- lapply(gene_sets, unique)
gene_sets <- gene_sets[vapply(gene_sets, function(gs) {
    k <- length(intersect(gs, rownames(expr_mat)))
    k >= 5 & k <= 500
}, logical(1))]
cat(sprintf("Gene sets: %d (after size filter 5-500)\n", length(gene_sets)))

# Pre-intersect gene sets with expression genes (needed for UCell/AUCell)
gene_sets_intersected <- lapply(gene_sets, function(gs) intersect(gs, rownames(expr_mat)))
gene_sets_intersected <- gene_sets_intersected[vapply(gene_sets_intersected, length, integer(1)) >= 5]
cat(sprintf("Gene sets (intersected with expr): %d\n", length(gene_sets_intersected)))

# =============================================================================
# 2. Run all tools
# =============================================================================

# --- scGSVA (ssGSEA, custom C++) ---
cat("\n--- scGSVA ssGSEA (C++ backend) ---\n")
t_scgsva_ssgsea <- system.time({
    res_scgsva_ssgsea <- scgsva(pbmcs, hsko, method = "ssgsea", cores = 1, use.original = FALSE)
})

# --- scGSVA (GSVA, custom C++) ---
cat("--- scGSVA GSVA (C++ backend) ---\n")
t_scgsva_gsva <- system.time({
    res_scgsva_gsva <- scgsva(pbmcs, hsko, method = "gsva", cores = 1, use.original = FALSE)
})

# --- Original GSVA ssGSEA ---
cat("--- Original GSVA ssGSEA ---\n")
t_gsva_ssgsea <- system.time({
    res_gsva_ssgsea <- scgsva(pbmcs, hsko, method = "ssgsea", cores = 1, use.original = TRUE)
})

# --- Original GSVA ---
cat("--- Original GSVA ---\n")
t_gsva_gsva <- system.time({
    res_gsva_gsva <- scgsva(pbmcs, hsko, method = "gsva", cores = 1, use.original = TRUE)
})

# --- AUCell ---
cat("--- AUCell ---\n")
library(AUCell)
t_aucell <- system.time({
    rankings_aucell <- AUCell_buildRankings(expr_mat, plotStats = FALSE, verbose = FALSE)
    auc_scores <- AUCell_calcAUC(gene_sets_intersected, rankings_aucell, verbose = FALSE)
    aucell_mat <- getAUC(auc_scores)
})

# --- UCell ---
cat("--- UCell ---\n")
library(UCell)
t_ucell <- system.time({
    ucell_scores <- ScoreSignatures_UCell(expr_mat, features = gene_sets_intersected,
                                          ncores = 1, maxRank = nrow(expr_mat))
})

# --- singscore ---
cat("--- singscore ---\n")
library(singscore)
library(SummarizedExperiment)
t_singscore <- system.time({
    # singscore needs a SummarizedExperiment and GeneSetCollection
    se <- SummarizedExperiment(assays = list(counts = expr_mat))
    ranked_se <- rankGenes(se)
    # Score each gene set
    singscore_list <- lapply(names(gene_sets_intersected), function(gs_name) {
        gs_genes <- gene_sets_intersected[[gs_name]]
        if (length(gs_genes) < 2) return(NULL)
        tryCatch({
            res <- simpleScore(ranked_se, upSet = gs_genes)
            score <- res$TotalScore
            names(score) <- rownames(res)  # cell barcodes
            score
        }, error = function(e) NULL)
    })
    names(singscore_list) <- names(gene_sets_intersected)
    singscore_list <- singscore_list[!vapply(singscore_list, is.null, logical(1))]
    singscore_mat <- do.call(rbind, singscore_list)
})

# =============================================================================
# 3. Compile benchmark results
# =============================================================================
cat("\n=== Benchmark Results ===\n")
benchmark_tools <- data.frame(
    Tool = c("scGSVA ssGSEA", "scGSVA GSVA", "GSVA ssGSEA", "GSVA (original)",
             "AUCell", "UCell", "singscore"),
    Time_sec = round(c(
        t_scgsva_ssgsea["elapsed"], t_scgsva_gsva["elapsed"],
        t_gsva_ssgsea["elapsed"], t_gsva_gsva["elapsed"],
        t_aucell["elapsed"], t_ucell["elapsed"], t_singscore["elapsed"]
    ), 3),
    stringsAsFactors = FALSE
)
benchmark_tools$Speedup_vs_GSVA_ssGSEA <- round(
    t_gsva_ssgsea["elapsed"] / pmax(benchmark_tools$Time_sec, 0.001), 1
)
print(benchmark_tools)
write.csv(benchmark_tools, file.path(results_dir, "benchmark_all_tools.csv"), row.names = FALSE)

# =============================================================================
# 4. Cross-tool correlation analysis
# =============================================================================
cat("\n=== Cross-tool Correlations ===\n")

# Extract score matrices (pathways x cells)
# scGSVA stores as cells x pathways data.frame; transpose to pathways x cells
scgsva_ssgsea_df <- res_scgsva_ssgsea@gsva
scgsva_gsva_df   <- res_scgsva_gsva@gsva

# Normalize pathway names: scGSVA uses dots, others use spaces
# Create lookup: canonical name (spaces) -> scGSVA name (dots)
scgsva_path_names <- colnames(scgsva_ssgsea_df)
canonical_from_scgsva <- gsub("\\.", " ", scgsva_path_names)
scgsva_name_map <- setNames(scgsva_path_names, canonical_from_scgsva)

# AUCell and singscore use spaces directly
# UCell adds "_UCell" suffix to space-based names
aucell_paths <- rownames(aucell_mat)
ucell_paths_raw <- sub("_UCell$", "", colnames(ucell_scores))
singscore_paths <- rownames(singscore_mat)

# Find common pathways (using canonical space-based names)
common_paths <- Reduce(intersect, list(
    canonical_from_scgsva,
    aucell_paths,
    ucell_paths_raw,
    singscore_paths
))
cat(sprintf("Common pathways across all tools: %d\n", length(common_paths)))

if (length(common_paths) == 0) {
    cat("WARNING: No common pathways found. Attempting pairwise comparisons only.\n")
    # Use scGSVA vs AUCell common paths as fallback
    common_paths <- intersect(canonical_from_scgsva, aucell_paths)
    cat(sprintf("scGSVA vs AUCell common: %d\n", length(common_paths)))
}

# Build per-pathway correlation matrix
tool_names <- c("scGSVA_ssGSEA", "scGSVA_GSVA", "AUCell", "UCell", "singscore")
n_tools <- length(tool_names)
cor_matrix <- matrix(NA, nrow = n_tools, ncol = n_tools,
                     dimnames = list(tool_names, tool_names))

# Get common cells
cells <- colnames(expr_mat)

# Extract aligned matrices (pathways x cells) using canonical (space-based) names
get_scores <- function(tool, paths = common_paths) {
    switch(tool,
        "scGSVA_ssGSEA" = {
            scgsva_cols <- scgsva_name_map[paths]
            as.matrix(t(scgsva_ssgsea_df[cells, scgsva_cols, drop = FALSE]))
        },
        "scGSVA_GSVA" = {
            scgsva_cols <- scgsva_name_map[paths]
            as.matrix(t(scgsva_gsva_df[cells, scgsva_cols, drop = FALSE]))
        },
        "AUCell"    = as.matrix(aucell_mat[paths, cells, drop = FALSE]),
        "UCell"     = {
            ucell_cols <- paste0(paths, "_UCell")
            t(as.matrix(ucell_scores[cells, ucell_cols, drop = FALSE]))
        },
        "singscore" = singscore_mat[paths, cells, drop = FALSE]
    )
}

# Compute mean per-pathway Pearson correlation between all tool pairs
cor_df <- data.frame()
for (i in seq_len(n_tools)) {
    for (j in seq_len(n_tools)) {
        mat_i <- get_scores(tool_names[i])
        mat_j <- get_scores(tool_names[j])
        cors <- sapply(seq_len(nrow(mat_i)), function(k) {
            cor(mat_i[k, ], mat_j[k, ], use = "complete.obs")
        })
        mean_cor <- mean(cors, na.rm = TRUE)
        cor_matrix[i, j] <- round(mean_cor, 4)
        if (i < j) {
            cor_df <- rbind(cor_df, data.frame(
                Tool_A = tool_names[i], Tool_B = tool_names[j],
                Mean_Cor = round(mean_cor, 4),
                Min_Cor = round(min(cors, na.rm = TRUE), 4),
                Max_Cor = round(max(cors, na.rm = TRUE), 4),
                stringsAsFactors = FALSE
            ))
        }
    }
}

cat("\nMean per-pathway correlation matrix:\n")
print(cor_matrix)
write.csv(cor_matrix, file.path(results_dir, "cross_tool_correlation_matrix.csv"))
write.csv(cor_df, file.path(results_dir, "cross_tool_pairwise_correlations.csv"), row.names = FALSE)

# =============================================================================
# 5. Feature comparison table
# =============================================================================
cat("\n=== Feature Comparison Table ===\n")
feature_df <- data.frame(
    Feature = c(
        "Language",
        "Backend",
        "Single-cell native",
        "Seurat integration",
        "AnnData/Scanpy support",
        "Methods available",
        "Annotation builder",
        "Visualization built-in",
        "Statistical testing",
        "Spatial transcriptomics",
        "Species coverage",
        "Cross-platform (R+Python)",
        "GSVA-equivalent output",
        "Parallel execution"
    ),
    scGSVA = c(
        "R + C++",
        "Rcpp/C++ (raw pointers, R math)",
        "Yes (per-cell scoring)",
        "Yes (native)",
        "Yes (via pygsva)",
        "ssGSEA, GSVA, PLAGE, z-score, UCell",
        "Yes (KEGG/GO/MSigDB, all species)",
        "Yes (6+ plot types)",
        "Yes (findPathway, sigPathway)",
        "Yes (spatialFeaturePlot)",
        "All (via AnnotationHub)",
        "Yes (scGSVA + pygsva)",
        "Yes (r=1.0 validated)",
        "Yes (BiocParallel chunks)"
    ),
    GSVA = c(
        "R",
        "R + C (kernel_estimation.c)",
        "Yes (v2.0+)",
        "Manual extraction",
        "No",
        "ssGSEA, GSVA, PLAGE, z-score",
        "No",
        "No",
        "No",
        "No",
        "N/A (user provides sets)",
        "No",
        "Reference",
        "Yes (BiocParallel)"
    ),
    AUCell = c(
        "R",
        "R",
        "Yes",
        "Manual extraction",
        "No (pySCENIC separate)",
        "AUC-based only",
        "No",
        "Threshold plots",
        "No",
        "No",
        "N/A",
        "No",
        "No (different method)",
        "No"
    ),
    UCell = c(
        "R",
        "R (rank-based)",
        "Yes",
        "Yes (AddModuleScore_UCell)",
        "No",
        "Mann-Whitney U only",
        "No",
        "No",
        "No",
        "No",
        "N/A",
        "No",
        "No (different method)",
        "Yes (multicore)"
    ),
    singscore = c(
        "R",
        "R",
        "Yes",
        "Manual extraction",
        "No",
        "Rank-based scoring only",
        "No",
        "Yes (landscape plots)",
        "Yes (permutation test)",
        "No",
        "N/A",
        "No",
        "No (different method)",
        "No"
    ),
    stringsAsFactors = FALSE
)
write.csv(feature_df, file.path(results_dir, "feature_comparison_table.csv"), row.names = FALSE)
cat("Feature comparison table saved.\n")

# =============================================================================
# 6. Figures
# =============================================================================
cat("\n=== Generating Figures ===\n")

# --- Figure A: Runtime comparison (all tools, bar plot) ---
benchmark_tools$Tool <- factor(benchmark_tools$Tool,
    levels = benchmark_tools$Tool[order(benchmark_tools$Time_sec)])

p_runtime <- ggplot(benchmark_tools, aes(x = Tool, y = Time_sec, fill = Tool)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = paste0(Time_sec, "s")), hjust = -0.1, size = 3) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Runtime Comparison: Single-Cell Gene Set Scoring Tools",
         subtitle = sprintf("PBMC dataset (%d genes x %d cells, %d KEGG sets)",
                            nrow(expr_mat), ncol(expr_mat), length(gene_sets)),
         x = "", y = "Time (seconds)") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "tool_runtime_comparison.pdf"), p_runtime, width = 10, height = 5)
ggsave(file.path(results_dir, "tool_runtime_comparison.png"), p_runtime, width = 10, height = 5, dpi = 300)

# --- Figure B: Cross-tool correlation heatmap ---
cor_long <- as.data.frame(as.table(cor_matrix))
names(cor_long) <- c("Tool_A", "Tool_B", "Correlation")
cor_long$Correlation <- as.numeric(cor_long$Correlation)

p_corheat <- ggplot(cor_long, aes(x = Tool_A, y = Tool_B, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0.5, limits = c(0, 1)) +
    labs(title = "Cross-Tool Score Correlation (Mean Per-Pathway Pearson r)",
         x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "cross_tool_correlation_heatmap.pdf"), p_corheat, width = 7, height = 6)
ggsave(file.path(results_dir, "cross_tool_correlation_heatmap.png"), p_corheat, width = 7, height = 6, dpi = 300)

# --- Figure C: Scatter plots (scGSVA ssGSEA vs each tool) ---
pdf(file.path(results_dir, "scatter_scgsva_vs_tools.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

others <- c("AUCell", "UCell", "singscore", "scGSVA_GSVA")
ref_mat <- get_scores("scGSVA_ssGSEA")

for (other in others) {
    other_mat <- get_scores(other)
    ref_vals <- as.vector(ref_mat)
    other_vals <- as.vector(other_mat)
    valid <- is.finite(ref_vals) & is.finite(other_vals)
    r <- round(cor(ref_vals[valid], other_vals[valid]), 4)
    plot(ref_vals[valid], other_vals[valid],
         pch = 16, cex = 0.4, col = rgb(0.2, 0.4, 0.8, 0.3),
         xlab = "scGSVA ssGSEA", ylab = other,
         main = paste0("scGSVA ssGSEA vs ", other, " (r = ", r, ")"))
    abline(lm(other_vals[valid] ~ ref_vals[valid]), col = "red", lwd = 2, lty = 2)
}
dev.off()

png(file.path(results_dir, "scatter_scgsva_vs_tools.png"), width = 1200, height = 800, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (other in others) {
    other_mat <- get_scores(other)
    ref_vals <- as.vector(ref_mat)
    other_vals <- as.vector(other_mat)
    valid <- is.finite(ref_vals) & is.finite(other_vals)
    r <- round(cor(ref_vals[valid], other_vals[valid]), 4)
    plot(ref_vals[valid], other_vals[valid],
         pch = 16, cex = 0.4, col = rgb(0.2, 0.4, 0.8, 0.3),
         xlab = "scGSVA ssGSEA", ylab = other,
         main = paste0("scGSVA ssGSEA vs ", other, " (r = ", r, ")"))
    abline(lm(other_vals[valid] ~ ref_vals[valid]), col = "red", lwd = 2, lty = 2)
}
dev.off()

# --- Figure D: Per-pathway correlation distribution (scGSVA ssGSEA vs each tool) ---
pathway_cors <- data.frame()
for (other in c("AUCell", "UCell", "singscore")) {
    ref_mat_tmp <- get_scores("scGSVA_ssGSEA")
    other_mat <- get_scores(other)
    cors <- sapply(seq_len(nrow(ref_mat_tmp)), function(k) {
        cor(ref_mat_tmp[k, ], other_mat[k, ], use = "complete.obs")
    })
    pathway_cors <- rbind(pathway_cors, data.frame(
        Tool = other, Correlation = cors, stringsAsFactors = FALSE
    ))
}

p_cordist <- ggplot(pathway_cors, aes(x = Correlation, fill = Tool)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    facet_wrap(~Tool, ncol = 1) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Per-Pathway Correlation: scGSVA ssGSEA vs Other Tools",
         x = "Pearson Correlation (per pathway)", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "pathway_cor_scgsva_vs_tools.pdf"), p_cordist, width = 7, height = 8)
ggsave(file.path(results_dir, "pathway_cor_scgsva_vs_tools.png"), p_cordist, width = 7, height = 8, dpi = 300)

# =============================================================================
# 7. Big data comparison (2000 genes x 500 cells)
# =============================================================================
cat("\n=== Big Data Comparison (2000 genes x 500 cells) ===\n")
set.seed(42)
counts <- as.matrix(Seurat::GetAssayData(pbmcs, assay = "RNA", layer = "counts"))
n_genes_target <- 2000
n_cells_target <- 500

gene_idx <- sample(nrow(counts), n_genes_target, replace = TRUE)
cell_idx <- sample(ncol(counts), n_cells_target, replace = TRUE)
big_counts <- counts[gene_idx, cell_idx]
rownames(big_counts) <- paste0("Gene", seq_len(n_genes_target))
colnames(big_counts) <- paste0("Cell", seq_len(n_cells_target))
big_counts <- big_counts + matrix(rpois(n_genes_target * n_cells_target, lambda = 0.5),
                                   nrow = n_genes_target)

big_seu <- Seurat::CreateSeuratObject(counts = big_counts)

# Create mapped gene sets for big data
orig_genes <- rownames(counts)
gene_map <- orig_genes[gene_idx]
names(gene_map) <- rownames(big_counts)
big_gene_sets <- lapply(gene_sets, function(gs) {
    mapped <- names(gene_map)[gene_map %in% gs]
    if (length(mapped) >= 5) mapped else character(0)
})
big_gene_sets <- big_gene_sets[vapply(big_gene_sets, length, integer(1)) >= 5]

# Build Annot for big data
annot_df_mapped <- do.call(rbind, lapply(seq_along(gene_map), function(i) {
    matches <- annot_df[annot_df$GeneID == gene_map[i], ]
    if (nrow(matches) > 0) { matches$GeneID <- names(gene_map)[i]; matches } else NULL
}))
big_annot <- new("Annot", species = "human", anntype = "KEGG", keytype = "SYMBOL", annot = annot_df_mapped)

# Run big data benchmarks
cat("Running big data scGSVA ssGSEA...\n")
t_big_scgsva <- system.time(scgsva(big_seu, big_annot, method = "ssgsea", cores = 1, use.original = FALSE, batch = 5000))

cat("Running big data GSVA ssGSEA...\n")
t_big_gsva <- system.time(scgsva(big_seu, big_annot, method = "ssgsea", cores = 1, use.original = TRUE, batch = 5000))

cat("Running big data AUCell...\n")
big_expr <- as.matrix(Seurat::GetAssayData(big_seu, assay = "RNA", layer = "counts"))
t_big_aucell <- system.time({
    rnk <- AUCell_buildRankings(big_expr, plotStats = FALSE, verbose = FALSE)
    AUCell_calcAUC(big_gene_sets, rnk, verbose = FALSE)
})

cat("Running big data UCell...\n")
t_big_ucell <- system.time({
    ScoreSignatures_UCell(big_expr, features = big_gene_sets, ncores = 1,
                           maxRank = nrow(big_expr))
})

big_bench <- data.frame(
    Tool = c("scGSVA ssGSEA", "GSVA ssGSEA", "AUCell", "UCell"),
    Time_sec = round(c(t_big_scgsva["elapsed"], t_big_gsva["elapsed"],
                       t_big_aucell["elapsed"], t_big_ucell["elapsed"]), 3),
    stringsAsFactors = FALSE
)
big_bench$Speedup_vs_GSVA <- round(t_big_gsva["elapsed"] / pmax(big_bench$Time_sec, 0.001), 1)
cat("\nBig Data Benchmark:\n")
print(big_bench)
write.csv(big_bench, file.path(results_dir, "benchmark_big_all_tools.csv"), row.names = FALSE)

# --- Figure E: Big data runtime comparison ---
big_bench$Tool <- factor(big_bench$Tool, levels = big_bench$Tool[order(big_bench$Time_sec)])

p_bigrt <- ggplot(big_bench, aes(x = Tool, y = Time_sec, fill = Tool)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = paste0(Time_sec, "s")), hjust = -0.1, size = 3.5) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Runtime: Big Data (2000 genes x 500 cells)",
         x = "", y = "Time (seconds)") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "tool_runtime_big_data.pdf"), p_bigrt, width = 9, height = 4)
ggsave(file.path(results_dir, "tool_runtime_big_data.png"), p_bigrt, width = 9, height = 4, dpi = 300)

# =============================================================================
# 8. Summary
# =============================================================================
cat("\n===================================================================\n")
cat("  All comparison results saved to:", results_dir, "\n")
cat("===================================================================\n")
cat("\nFiles generated:\n")
cat("  Tables:\n")
cat("    - benchmark_all_tools.csv (small data, all tools)\n")
cat("    - benchmark_big_all_tools.csv (big data)\n")
cat("    - cross_tool_correlation_matrix.csv\n")
cat("    - cross_tool_pairwise_correlations.csv\n")
cat("    - feature_comparison_table.csv\n")
cat("  Figures:\n")
cat("    - tool_runtime_comparison.pdf/png (small data)\n")
cat("    - tool_runtime_big_data.pdf/png (big data)\n")
cat("    - cross_tool_correlation_heatmap.pdf/png\n")
cat("    - scatter_scgsva_vs_tools.pdf/png\n")
cat("    - pathway_cor_scgsva_vs_tools.pdf/png\n")
