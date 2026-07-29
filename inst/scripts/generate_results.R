#!/usr/bin/env Rscript
# Generate comprehensive comparison results: Custom C++ vs Original GSVA Package
# Output: results/ folder with tables (CSV) and figures (PDF/PNG)

library(devtools)
load_all("/Users/bioguo/Library/Mobile Documents/com~apple~CloudDocs/Work/Packages/scGSVA")
library(ggplot2)
library(tidyr)
library(dplyr)

# Create results directory
results_dir <- "/Users/bioguo/Library/Mobile Documents/com~apple~CloudDocs/Work/Packages/scGSVA/results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
cat("Results directory:", results_dir, "\n\n")

# Load data
data(pbmcs)
hsko <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")

# ============================================================================
# 1. Run all methods and collect results
# ============================================================================
methods <- c("ssgsea", "gsva", "plage", "zscore")
all_results <- list()
benchmark_df <- data.frame()

for (method in methods) {
    cat("Running", toupper(method), "...\n")

    t_custom <- system.time(res_custom <- scgsva(pbmcs, hsko, method = method, cores = 1, use.original = FALSE))
    t_orig   <- system.time(res_orig   <- scgsva(pbmcs, hsko, method = method, cores = 1, use.original = TRUE))

    common <- intersect(colnames(res_custom@gsva), colnames(res_orig@gsva))
    cors <- sapply(common, function(p) {
        cor(res_custom@gsva[[p]], res_orig@gsva[[p]], use = "complete.obs")
    })
    diffs <- sapply(common, function(p) {
        max(abs(res_custom@gsva[[p]] - res_orig@gsva[[p]]), na.rm = TRUE)
    })

    all_results[[method]] <- list(
        custom = res_custom, original = res_orig,
        common = common, cors = cors, diffs = diffs,
        t_custom = t_custom["elapsed"], t_orig = t_orig["elapsed"]
    )

    benchmark_df <- rbind(benchmark_df, data.frame(
        Method = toupper(method),
        Custom_Time = round(t_custom["elapsed"], 3),
        Original_Time = round(t_orig["elapsed"], 3),
        Speedup = round(t_orig["elapsed"] / max(t_custom["elapsed"], 0.001), 1),
        Mean_Cor = round(mean(cors, na.rm = TRUE), 6),
        Min_Cor = round(min(cors, na.rm = TRUE), 6),
        N_Pathways = length(common),
        Cor_gt_099 = sum(cors >= 0.99, na.rm = TRUE),
        Max_Abs_Diff = round(max(diffs, na.rm = TRUE), 6),
        stringsAsFactors = FALSE
    ))
}

# ============================================================================
# 2. Save benchmark summary table
# ============================================================================
write.csv(benchmark_df, file.path(results_dir, "benchmark_summary.csv"), row.names = FALSE)
cat("\nBenchmark Summary:\n")
print(benchmark_df)

# ============================================================================
# 3. Per-pathway correlation tables for each method
# ============================================================================
for (method in methods) {
    r <- all_results[[method]]
    pathway_df <- data.frame(
        Pathway = names(r$cors),
        Correlation = round(r$cors, 6),
        Max_Abs_Diff = round(r$diffs, 6),
        stringsAsFactors = FALSE
    )
    pathway_df <- pathway_df[order(pathway_df$Correlation), ]
    write.csv(pathway_df, file.path(results_dir, paste0("pathway_cor_", method, ".csv")), row.names = FALSE)
}

# ============================================================================
# 4. Score comparison tables (first 10 pathways x first 10 cells)
# ============================================================================
for (method in methods) {
    r <- all_results[[method]]
    common <- r$common
    show_paths <- head(common, 10)
    show_cells <- head(rownames(r$custom@gsva), 10)

    custom_sub <- round(r$custom@gsva[show_cells, show_paths], 4)
    orig_sub   <- round(r$original@gsva[show_cells, show_paths], 4)

    write.csv(custom_sub, file.path(results_dir, paste0("scores_custom_", method, ".csv")))
    write.csv(orig_sub,   file.path(results_dir, paste0("scores_original_", method, ".csv")))
}

# ============================================================================
# 5. Figures
# ============================================================================

# --- Figure 1: Benchmark bar plot (time comparison) ---
bench_long <- benchmark_df[, c("Method", "Custom_Time", "Original_Time")]
bench_long <- tidyr::pivot_longer(bench_long, cols = c("Custom_Time", "Original_Time"),
                 names_to = "Implementation", values_to = "Time")
bench_long$Implementation <- ifelse(bench_long$Implementation == "Custom_Time", "Custom (C++)", "Original GSVA")

p1 <- ggplot(bench_long, aes(x = Method, y = Time, fill = Implementation)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    geom_text(aes(label = paste0(Time, "s")), position = position_dodge(0.6),
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("Custom (C++)" = "#2196F3", "Original GSVA" = "#FF5722")) +
    labs(title = "Performance: Custom C++ vs Original GSVA Package",
         subtitle = paste0("Dataset: pbmcs (", ncol(pbmcs), " cells, ", nrow(pbmcs), " genes)"),
         x = "Method", y = "Time (seconds)", fill = "") +
    theme_minimal() +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "benchmark_barplot.pdf"), p1, width = 8, height = 5)
ggsave(file.path(results_dir, "benchmark_barplot.png"), p1, width = 8, height = 5, dpi = 300)

# --- Figure 2: Speedup bar plot ---
p2 <- ggplot(benchmark_df, aes(x = Method, y = Speedup, fill = Method)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = paste0(Speedup, "x")), vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Speedup: Custom C++ over Original GSVA Package",
         x = "Method", y = "Speedup (x times faster)") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "speedup_barplot.pdf"), p2, width = 7, height = 5)
ggsave(file.path(results_dir, "speedup_barplot.png"), p2, width = 7, height = 5, dpi = 300)

# --- Figure 3: Scatter plots (Custom vs Original) for each method ---
pdf(file.path(results_dir, "scatter_custom_vs_original.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (method in methods) {
    r <- all_results[[method]]
    # Flatten all scores for plotting
    custom_vals <- unlist(r$custom@gsva[, r$common])
    orig_vals   <- unlist(r$original@gsva[, r$common])
    valid <- is.finite(custom_vals) & is.finite(orig_vals)

    plot(orig_vals[valid], custom_vals[valid],
         pch = 16, cex = 0.3, col = rgb(0.2, 0.4, 0.8, 0.3),
         xlab = "Original GSVA Package", ylab = "Custom (C++)",
         main = paste0(toupper(method), " (cor = ",
                       round(cor(custom_vals[valid], orig_vals[valid]), 4), ")"))
    abline(0, 1, col = "red", lwd = 2, lty = 2)
}

dev.off()

# Also as PNG
png(file.path(results_dir, "scatter_custom_vs_original.png"), width = 1200, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (method in methods) {
    r <- all_results[[method]]
    custom_vals <- unlist(r$custom@gsva[, r$common])
    orig_vals   <- unlist(r$original@gsva[, r$common])
    valid <- is.finite(custom_vals) & is.finite(orig_vals)

    plot(orig_vals[valid], custom_vals[valid],
         pch = 16, cex = 0.3, col = rgb(0.2, 0.4, 0.8, 0.3),
         xlab = "Original GSVA Package", ylab = "Custom (C++)",
         main = paste0(toupper(method), " (cor = ",
                       round(cor(custom_vals[valid], orig_vals[valid]), 4), ")"))
    abline(0, 1, col = "red", lwd = 2, lty = 2)
}

dev.off()

# --- Figure 4: Per-pathway correlation distribution ---
cor_df <- data.frame()
for (method in methods) {
    r <- all_results[[method]]
    cor_df <- rbind(cor_df, data.frame(
        Method = toupper(method),
        Correlation = r$cors,
        stringsAsFactors = FALSE
    ))
}

p4 <- ggplot(cor_df, aes(x = Correlation, fill = Method)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    facet_wrap(~Method, scales = "free_y") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Per-Pathway Correlation Distribution (Custom vs Original)",
         x = "Pearson Correlation", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "correlation_distribution.pdf"), p4, width = 10, height = 7)
ggsave(file.path(results_dir, "correlation_distribution.png"), p4, width = 10, height = 7, dpi = 300)

# --- Figure 5: Heatmap comparison for GSVA method ---
r <- all_results[["gsva"]]
show_paths <- head(r$common, 20)
show_cells <- head(rownames(r$custom@gsva), 30)

pdf(file.path(results_dir, "heatmap_gsva_comparison.pdf"), width = 14, height = 8)
par(mfrow = c(1, 2), mar = c(5, 8, 3, 2))

custom_mat <- as.matrix(r$custom@gsva[show_cells, show_paths])
orig_mat   <- as.matrix(r$original@gsva[show_cells, show_paths])

image(t(custom_mat), col = colorRampPalette(c("blue", "white", "red"))(100),
      axes = FALSE, main = "Custom (C++)")
axis(2, at = seq(0, 1, length.out = ncol(custom_mat)),
     labels = colnames(custom_mat), las = 2, cex.axis = 0.6)

image(t(orig_mat), col = colorRampPalette(c("blue", "white", "red"))(100),
      axes = FALSE, main = "Original GSVA Package")
axis(2, at = seq(0, 1, length.out = ncol(orig_mat)),
     labels = colnames(orig_mat), las = 2, cex.axis = 0.6)

dev.off()

# ============================================================================
# 6. Big Data Benchmark (2000 genes x 500 cells)
# ============================================================================
cat("\n=== Big Data Benchmark (2000 genes x 500 cells) ===\n")
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

# Create mapped annotation
orig_genes <- rownames(counts)
gene_map <- orig_genes[gene_idx]
names(gene_map) <- rownames(big_counts)
annot_df <- as.data.frame(hsko)
mapped_annot <- do.call(rbind, lapply(seq_along(gene_map), function(i) {
    matches <- annot_df[annot_df$GeneID == gene_map[i], ]
    if (nrow(matches) > 0) { matches$GeneID <- names(gene_map)[i]; matches } else NULL
}))
big_annot <- new("Annot", species = "human", anntype = "KEGG", keytype = "SYMBOL", annot = mapped_annot)

benchmark_big <- data.frame()
for (method in c("ssgsea", "gsva", "zscore")) {
    cat("Running big data", toupper(method), "...\n")
    t_custom <- system.time(res_custom <- scgsva(big_seu, big_annot, method = method, use.original = FALSE, batch = 5000))
    t_orig   <- system.time(res_orig   <- scgsva(big_seu, big_annot, method = method, use.original = TRUE, batch = 5000))
    common <- intersect(colnames(res_custom@gsva), colnames(res_orig@gsva))
    cors <- sapply(common, function(p) cor(res_custom@gsva[[p]], res_orig@gsva[[p]], use = "complete.obs"))
    spdup <- round(t_orig["elapsed"] / max(t_custom["elapsed"], 0.001), 1)
    benchmark_big <- rbind(benchmark_big, data.frame(
        Method = toupper(method),
        N_Cells = n_cells_target, N_Genes = n_genes_target,
        Custom_Time = round(t_custom["elapsed"], 3),
        Original_Time = round(t_orig["elapsed"], 3),
        Speedup = spdup,
        Mean_Cor = round(mean(cors, na.rm = TRUE), 6),
        Min_Cor = round(min(cors, na.rm = TRUE), 6),
        stringsAsFactors = FALSE
    ))
}
write.csv(benchmark_big, file.path(results_dir, "benchmark_big_data.csv"), row.names = FALSE)
cat("\nBig Data Benchmark:\n")
print(benchmark_big)

# --- Figure 6: Combined benchmark (small + big) ---
benchmark_big$Dataset <- paste0(n_genes_target, " genes x ", n_cells_target, " cells")
benchmark_df$Dataset <- paste0(nrow(pbmcs), " genes x ", ncol(pbmcs), " cells")
combined <- rbind(
    benchmark_df[, c("Method", "Custom_Time", "Original_Time", "Speedup", "Dataset")],
    benchmark_big[, c("Method", "Custom_Time", "Original_Time", "Speedup", "Dataset")]
)

comb_long <- tidyr::pivot_longer(combined[, c("Method", "Custom_Time", "Original_Time", "Dataset")],
                                  cols = c("Custom_Time", "Original_Time"),
                                  names_to = "Implementation", values_to = "Time")
comb_long$Implementation <- ifelse(comb_long$Implementation == "Custom_Time", "Custom (C++)", "Original GSVA")

p6 <- ggplot(comb_long, aes(x = Method, y = Time, fill = Implementation)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    geom_text(aes(label = paste0(round(Time, 2), "s")), position = position_dodge(0.6),
              vjust = -0.5, size = 2.5) +
    facet_wrap(~Dataset, scales = "free_y") +
    scale_fill_manual(values = c("Custom (C++)" = "#2196F3", "Original GSVA" = "#FF5722")) +
    labs(title = "Performance: Custom C++ vs Original GSVA Package",
         subtitle = "Small data (pbmcs) and Big data (simulated)",
         x = "Method", y = "Time (seconds)", fill = "") +
    theme_minimal() +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "benchmark_combined.pdf"), p6, width = 12, height = 5)
ggsave(file.path(results_dir, "benchmark_combined.png"), p6, width = 12, height = 5, dpi = 300)

cat("\n===================================================================\n")
cat("  All results saved to:", results_dir, "\n")
cat("===================================================================\n")
cat("\nFiles generated:\n")
cat("  Tables:\n")
cat("    - benchmark_summary.csv (small data)\n")
cat("    - benchmark_big_data.csv (big data)\n")
cat("    - pathway_cor_*.csv (per method)\n")
cat("    - scores_custom_*.csv / scores_original_*.csv (per method)\n")
cat("  Figures:\n")
cat("    - benchmark_barplot.pdf/png (small data)\n")
cat("    - speedup_barplot.pdf/png\n")
cat("    - scatter_custom_vs_original.pdf/png\n")
cat("    - correlation_distribution.pdf/png\n")
cat("    - heatmap_gsva_comparison.pdf\n")
cat("    - benchmark_combined.pdf/png (small + big data)\n")
