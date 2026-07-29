#!/usr/bin/env Rscript
# =============================================================================
# Scaling analysis: runtime vs number of cells and genes
# for scGSVA (ssGSEA, GSVA, z-score) vs GSVA(orig) vs AUCell vs UCell
#
# Tests:
#   A) Fixed genes (2000), vary cells: 100, 200, 500, 1000, 2000
#   B) Fixed cells (500), vary genes: 500, 1000, 2000, 5000
#
# Output: results/ folder with scaling tables and log-log curve figures
# =============================================================================

library(devtools)
load_all("/Users/bioguo/Library/Mobile Documents/com~apple~CloudDocs/Work/Packages/scGSVA")
library(ggplot2)
library(dplyr)
library(AUCell)
library(UCell)

results_dir <- "/Users/bioguo/Library/Mobile Documents/com~apple~CloudDocs/Work/Packages/scGSVA/results"
dir.create(results_dir, showWarnings = FALSE)

# =============================================================================
# Helper: generate synthetic data and gene sets
# =============================================================================
data(pbmcs)
hsko <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")
base_counts <- as.matrix(Seurat::GetAssayData(pbmcs, assay = "RNA", layer = "counts"))
annot_df <- as.data.frame(hsko)

make_data <- function(n_genes, n_cells) {
    set.seed(42)
    gene_idx <- sample(nrow(base_counts), n_genes, replace = TRUE)
    cell_idx <- sample(ncol(base_counts), n_cells, replace = TRUE)
    mat <- base_counts[gene_idx, cell_idx]
    rownames(mat) <- paste0("Gene", seq_len(n_genes))
    colnames(mat) <- paste0("Cell", seq_len(n_cells))
    mat <- mat + matrix(rpois(n_genes * n_cells, lambda = 0.5), nrow = n_genes)

    orig_genes <- rownames(base_counts)
    gene_map <- orig_genes[gene_idx]
    names(gene_map) <- rownames(mat)

    ad <- do.call(rbind, lapply(seq_along(gene_map), function(i) {
        matches <- annot_df[annot_df$GeneID == gene_map[i], ]
        if (nrow(matches) > 0) { matches$GeneID <- names(gene_map)[i]; matches } else NULL
    }))
    annot_obj <- new("Annot", species = "human", anntype = "KEGG", keytype = "SYMBOL", annot = ad)

    seu <- Seurat::CreateSeuratObject(counts = mat)

    # Gene set list for AUCell/UCell
    gs_full <- split(ad$GeneID, ad$Annot)
    gs_full <- lapply(gs_full, unique)
    gs_int <- lapply(gs_full, function(gs) intersect(gs, rownames(mat)))
    gs_int <- gs_int[vapply(gs_int, length, integer(1)) >= 5]

    list(seu = seu, annot = annot_obj, mat = mat, gene_sets = gs_int)
}

# =============================================================================
# Benchmark function: time a single tool
# =============================================================================
time_tool <- function(tool_name, data_list) {
    seu <- data_list$seu
    annot <- data_list$annot
    mat <- data_list$mat
    gs <- data_list$gene_sets
    ng <- nrow(mat)

    t <- tryCatch(
        system.time(switch(tool_name,
            "scGSVA_ssGSEA" = scgsva(seu, annot, method = "ssgsea", cores = 1, use.original = FALSE, batch = 5000),
            "scGSVA_GSVA"   = scgsva(seu, annot, method = "gsva",   cores = 1, use.original = FALSE, batch = 5000),
            "scGSVA_zscore" = scgsva(seu, annot, method = "zscore", cores = 1, use.original = FALSE, batch = 5000),
            "GSVA_ssGSEA"   = scgsva(seu, annot, method = "ssgsea", cores = 1, use.original = TRUE,  batch = 5000),
            "GSVA_GSVA"     = scgsva(seu, annot, method = "gsva",   cores = 1, use.original = TRUE,  batch = 5000),
            "GSVA_zscore"   = scgsva(seu, annot, method = "zscore", cores = 1, use.original = TRUE,  batch = 5000),
            "AUCell" = {
                rnk <- AUCell_buildRankings(mat, plotStats = FALSE, verbose = FALSE)
                AUCell_calcAUC(gs, rnk, verbose = FALSE)
            },
            "UCell" = ScoreSignatures_UCell(mat, features = gs, ncores = 1, maxRank = ng)
        )),
        error = function(e) { cat("  ERROR:", e$message, "\n"); c(elapsed = NA) }
    )
    t["elapsed"]
}

# =============================================================================
# A) Scaling with number of cells (fixed 2000 genes)
# =============================================================================
cat("=== A) Scaling with N_cells (2000 genes) ===\n")
n_cells_vec <- c(100, 200, 500, 1000, 2000)
n_genes_fixed <- 2000
tools <- c("scGSVA_ssGSEA", "scGSVA_GSVA", "scGSVA_zscore",
           "GSVA_ssGSEA", "GSVA_GSVA", "GSVA_zscore",
           "AUCell", "UCell")

scale_cells_df <- data.frame()
for (nc in n_cells_vec) {
    cat(sprintf("  N_cells = %d ... ", nc))
    d <- make_data(n_genes_fixed, nc)
    for (tool in tools) {
        cat(tool, "")
        elapsed <- time_tool(tool, d)
        scale_cells_df <- rbind(scale_cells_df, data.frame(
            N_Genes = n_genes_fixed, N_Cells = nc,
            Tool = tool, Time_sec = round(elapsed, 4),
            stringsAsFactors = FALSE
        ))
    }
    cat("\n")
}
write.csv(scale_cells_df, file.path(results_dir, "scaling_by_cells.csv"), row.names = FALSE)
cat("\nScaling by cells:\n")
print(tidyr::pivot_wider(scale_cells_df, names_from = Tool, values_from = Time_sec))

# =============================================================================
# B) Scaling with number of genes (fixed 500 cells)
# =============================================================================
cat("\n=== B) Scaling with N_genes (500 cells) ===\n")
n_genes_vec <- c(500, 1000, 2000, 5000)
n_cells_fixed <- 500

scale_genes_df <- data.frame()
for (ng in n_genes_vec) {
    cat(sprintf("  N_genes = %d ... ", ng))
    d <- make_data(ng, n_cells_fixed)
    for (tool in tools) {
        cat(tool, "")
        elapsed <- time_tool(tool, d)
        scale_genes_df <- rbind(scale_genes_df, data.frame(
            N_Genes = ng, N_Cells = n_cells_fixed,
            Tool = tool, Time_sec = round(elapsed, 4),
            stringsAsFactors = FALSE
        ))
    }
    cat("\n")
}
write.csv(scale_genes_df, file.path(results_dir, "scaling_by_genes.csv"), row.names = FALSE)
cat("\nScaling by genes:\n")
print(tidyr::pivot_wider(scale_genes_df, names_from = Tool, values_from = Time_sec))

# =============================================================================
# C) Figures
# =============================================================================
cat("\n=== Generating scaling figures ===\n")

# Tool categories for coloring
tool_cat <- c(
    "scGSVA_ssGSEA" = "scGSVA", "scGSVA_GSVA" = "scGSVA", "scGSVA_zscore" = "scGSVA",
    "GSVA_ssGSEA" = "GSVA(orig)", "GSVA_GSVA" = "GSVA(orig)", "GSVA_zscore" = "GSVA(orig)",
    "AUCell" = "AUCell", "UCell" = "UCell"
)
tool_lty <- c(
    "scGSVA_ssGSEA" = "solid", "scGSVA_GSVA" = "dashed", "scGSVA_zscore" = "dotted",
    "GSVA_ssGSEA" = "solid", "GSVA_GSVA" = "dashed", "GSVA_zscore" = "dotted",
    "AUCell" = "solid", "UCell" = "solid"
)
tool_colors <- c(
    "scGSVA_ssGSEA" = "#1f77b4", "scGSVA_GSVA" = "#2ca02c", "scGSVA_zscore" = "#17becf",
    "GSVA_ssGSEA" = "#ff7f0e", "GSVA_GSVA" = "#d62728", "GSVA_zscore" = "#bcbd22",
    "AUCell" = "#8c564b", "UCell" = "#7f7f7f"
)

# --- Figure A: Runtime vs N_cells (log-log) ---
scale_cells_df$Tool <- factor(scale_cells_df$Tool, levels = tools)

pA <- ggplot(scale_cells_df, aes(x = N_Cells, y = Time_sec, color = Tool, linetype = Tool)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_x_log10(breaks = n_cells_vec) +
    scale_y_log10() +
    scale_color_manual(values = tool_colors) +
    scale_linetype_manual(values = tool_lty) +
    labs(title = "Scaling: Runtime vs Number of Cells",
         subtitle = sprintf("Fixed %d genes, KEGG gene sets", n_genes_fixed),
         x = "Number of Cells (log scale)", y = "Runtime (seconds, log scale)") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"),
          legend.text = element_text(size = 8))

ggsave(file.path(results_dir, "scaling_by_cells.pdf"), pA, width = 10, height = 6)
ggsave(file.path(results_dir, "scaling_by_cells.png"), pA, width = 10, height = 6, dpi = 300)

# --- Figure B: Runtime vs N_genes (log-log) ---
scale_genes_df$Tool <- factor(scale_genes_df$Tool, levels = tools)

pB <- ggplot(scale_genes_df, aes(x = N_Genes, y = Time_sec, color = Tool, linetype = Tool)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_x_log10(breaks = n_genes_vec) +
    scale_y_log10() +
    scale_color_manual(values = tool_colors) +
    scale_linetype_manual(values = tool_lty) +
    labs(title = "Scaling: Runtime vs Number of Genes",
         subtitle = sprintf("Fixed %d cells, KEGG gene sets", n_cells_fixed),
         x = "Number of Genes (log scale)", y = "Runtime (seconds, log scale)") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"),
          legend.text = element_text(size = 8))

ggsave(file.path(results_dir, "scaling_by_genes.pdf"), pB, width = 10, height = 6)
ggsave(file.path(results_dir, "scaling_by_genes.png"), pB, width = 10, height = 6, dpi = 300)

# --- Figure C: Speedup ratio vs N_cells ---
speedup_cells <- scale_cells_df %>%
    tidyr::pivot_wider(names_from = Tool, values_from = Time_sec) %>%
    mutate(
        ssGSEA_speedup = GSVA_ssGSEA / pmax(scGSVA_ssGSEA, 0.001),
        GSVA_speedup   = GSVA_GSVA / pmax(scGSVA_GSVA, 0.001),
        zscore_speedup = GSVA_zscore / pmax(scGSVA_zscore, 0.001)
    ) %>%
    dplyr::select(N_Cells, ssGSEA_speedup, GSVA_speedup, zscore_speedup) %>%
    tidyr::pivot_longer(-N_Cells, names_to = "Method", values_to = "Speedup") %>%
    mutate(Method = sub("_speedup", "", Method))

pC <- ggplot(speedup_cells, aes(x = N_Cells, y = Speedup, color = Method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    scale_x_log10(breaks = n_cells_vec) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Speedup: scGSVA over GSVA(original) vs Data Size",
         subtitle = sprintf("Fixed %d genes", n_genes_fixed),
         x = "Number of Cells (log scale)", y = "Speedup (x times faster)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "speedup_vs_cells.pdf"), pC, width = 9, height = 5)
ggsave(file.path(results_dir, "speedup_vs_cells.png"), pC, width = 9, height = 5, dpi = 300)

# --- Figure D: Speedup ratio vs N_genes ---
speedup_genes <- scale_genes_df %>%
    tidyr::pivot_wider(names_from = Tool, values_from = Time_sec) %>%
    mutate(
        ssGSEA_speedup = GSVA_ssGSEA / pmax(scGSVA_ssGSEA, 0.001),
        GSVA_speedup   = GSVA_GSVA / pmax(scGSVA_GSVA, 0.001),
        zscore_speedup = GSVA_zscore / pmax(scGSVA_zscore, 0.001)
    ) %>%
    dplyr::select(N_Genes, ssGSEA_speedup, GSVA_speedup, zscore_speedup) %>%
    tidyr::pivot_longer(-N_Genes, names_to = "Method", values_to = "Speedup") %>%
    mutate(Method = sub("_speedup", "", Method))

pD <- ggplot(speedup_genes, aes(x = N_Genes, y = Speedup, color = Method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    scale_x_log10(breaks = n_genes_vec) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Speedup: scGSVA over GSVA(original) vs Gene Count",
         subtitle = sprintf("Fixed %d cells", n_cells_fixed),
         x = "Number of Genes (log scale)", y = "Speedup (x times faster)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "speedup_vs_genes.pdf"), pD, width = 9, height = 5)
ggsave(file.path(results_dir, "speedup_vs_genes.png"), pD, width = 9, height = 5, dpi = 300)

# =============================================================================
# Summary
# =============================================================================
cat("\n===================================================================\n")
cat("  Scaling analysis complete. Results saved to:", results_dir, "\n")
cat("===================================================================\n")
cat("Tables: scaling_by_cells.csv, scaling_by_genes.csv\n")
cat("Figures: scaling_by_cells.png, scaling_by_genes.png,\n")
cat("         speedup_vs_cells.png, speedup_vs_genes.png\n")
