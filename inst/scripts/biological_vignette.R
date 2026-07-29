#!/usr/bin/env Rscript
# Biological application vignette: pathway activity analysis of PBMC subsets
# Demonstrates scGSVA end-to-end workflow for the manuscript

library(scGSVA)
library(ggplot2)
library(dplyr)

# Use serial backend to avoid BiocParallel core limit during R CMD check
if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::register(BiocParallel::SerialParam())
}

set.seed(42)
data(pbmcs)

cat("=== Biological Application Vignette ===\n")
cat("Dataset: PBMC (", ncol(pbmcs), "cells,", nrow(pbmcs), "genes, groups:",
    paste(unique(pbmcs$groups), collapse=" vs "), ")\n\n")

# --- Step 1: Build KEGG annotation ---
hsko <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")
cat("Annotation: ", length(hsko@annot), "KEGG pathways\n")

# --- Step 2: Run all methods via scGSVA ---
cat("\nRunning scGSVA with ssGSEA...\n")
t0 <- proc.time()
res_ssgsea <- scgsva(pbmcs, hsko, method = "ssgsea", use.original = FALSE)
t_ssgsea <- (proc.time() - t0)[3]

cat("Running scGSVA with GSVA...\n")
t0 <- proc.time()
res_gsva <- scgsva(pbmcs, hsko, method = "gsva", use.original = FALSE)
t_gsva <- (proc.time() - t0)[3]

cat("Running scGSVA with z-score...\n")
t0 <- proc.time()
res_zscore <- scgsva(pbmcs, hsko, method = "zscore", use.original = FALSE)
t_zscore <- (proc.time() - t0)[3]

cat(sprintf("  ssGSEA: %.3fs, GSVA: %.3fs, z-score: %.3fs\n",
            t_ssgsea, t_gsva, t_zscore))

# --- Step 3: Differential pathway analysis ---
cat("\nRunning differential pathway analysis (g1 vs g2)...\n")
sig <- sigPathway(res_ssgsea, group = "groups")
sig <- sig[sig$p.adj < 0.05, ]
cat("Significant pathways (ssGSEA, adjusted p < 0.05):", nrow(sig), "\n")

if (nrow(sig) > 0) {
    cat("\nTop 10 differential pathways:\n")
    top10 <- head(sig, 10)
    for (i in seq_len(nrow(top10))) {
        cat(sprintf("  %2d. %-45s p.adj=%.2e\n",
                    i, top10$Path[i], top10$p.adj[i]))
    }
}

# --- Step 4: Generate publication-quality figures ---
outdir <- "results/vignette"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 4a. Violin plot of top pathways
cat("\nGenerating violin plot...\n")
scores <- res_ssgsea@gsva
if (nrow(sig) >= 3) {
    top_paths <- head(sig$Path, 6)
} else {
    # Use pathways with highest variance
    vars <- apply(scores, 2, var)
    top_paths <- names(sort(vars, decreasing = TRUE))[1:6]
}

p_vln <- vlnPlot(res_ssgsea, features = top_paths[1:min(3, length(top_paths))],
                 group_by = "groups")
ggsave(file.path(outdir, "vln_top_pathways.png"), p_vln,
       width = 10, height = 4, dpi = 300)

# 4b. Ridge plot
cat("Generating ridge plot...\n")
p_ridge <- ridgePlot(res_ssgsea, features = top_paths[1:min(3, length(top_paths))],
                     group_by = "groups")
ggsave(file.path(outdir, "ridge_top_pathways.png"), p_ridge,
       width = 10, height = 4, dpi = 300)

# 4c. Dot plot of top 15 pathways
cat("Generating dot plot...\n")
if (nrow(sig) >= 3) {
    p_dot <- dotPlot(res_ssgsea, features = head(sig$Path, 15),
                     group_by = "groups")
    ggsave(file.path(outdir, "dot_top15_pathways.png"), p_dot,
           width = 10, height = 6, dpi = 300)
}

# 4d. Heatmap of pathway scores
cat("Generating heatmap...\n")
if (nrow(sig) >= 5) {
    png(file.path(outdir, "heatmap_top20_pathways.png"), width = 1000, height = 800)
    Heatmap(res_ssgsea, features = head(sig$Path, 20),
            group_by = "groups")
    dev.off()
}

# 4e. Cross-method concordance on top pathways
cat("Generating cross-method concordance...\n")
scores_ss <- res_ssgsea@gsva
scores_gs <- res_gsva@gsva
scores_zs <- res_zscore@gsva
common <- intersect(intersect(colnames(scores_ss), colnames(scores_gs)),
                    colnames(scores_zs))
if (length(common) >= 5) {
    cors_gs <- sapply(common, function(p)
        cor(scores_ss[, p], scores_gs[, p], use = "complete.obs"))
    cors_zs <- sapply(common, function(p)
        cor(scores_ss[, p], scores_zs[, p], use = "complete.obs"))

    cor_df <- data.frame(
        Pathway = rep(common, 2),
        Correlation = c(cors_gs, cors_zs),
        Comparison = rep(c("ssGSEA vs GSVA", "ssGSEA vs z-score"),
                         each = length(common))
    )

    p_cor <- ggplot(cor_df, aes(x = Correlation, fill = Comparison)) +
        geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
        theme_minimal(base_size = 12) +
        labs(title = "Cross-method concordance on PBMC pathways",
             x = "Pearson r", y = "Number of pathways") +
        scale_fill_manual(values = c("#3498db", "#e74c3c")) +
        theme(legend.position = "top")
    ggsave(file.path(outdir, "cross_method_concordance.png"), p_cor,
           width = 7, height = 5, dpi = 300)
}

# --- Step 5: Summary statistics table ---
cat("\nSaving summary statistics...\n")
summary_df <- data.frame(
    Method = c("ssGSEA", "GSVA", "z-score"),
    Runtime_seconds = round(c(t_ssgsea, t_gsva, t_zscore), 3),
    N_pathways = c(ncol(scores_ss), ncol(scores_gs), ncol(scores_zs)),
    stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(outdir, "vignette_summary.csv"), row.names = FALSE)

if (nrow(sig) > 0) {
    write.csv(sig, file.path(outdir, "differential_pathways.csv"), row.names = FALSE)
}

cat("\n=== Vignette complete ===\n")
cat("Output directory:", outdir, "\n")
cat("Files:", paste(list.files(outdir), collapse = ", "), "\n")
