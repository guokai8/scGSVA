#!/usr/bin/env Rscript
# =============================================================================
# Comprehensive all-vs-all comparison:
#   scGSVA (ssGSEA, GSVA, PLAGE, z-score) vs competing tools
#   (GSVA original, AUCell, UCell, singscore)
#
# Generates: results/ folder with tables (CSV) and figures (PDF/PNG)
# =============================================================================

library(devtools)
load_all("/Users/bioguo/Library/Mobile Documents/com~apple~CloudDocs/Work/Packages/scGSVA")
library(ggplot2)
library(tidyr)
library(dplyr)
library(AUCell)
library(UCell)
library(singscore)
library(SummarizedExperiment)

results_dir <- "/Users/bioguo/Library/Mobile Documents/com~apple~CloudDocs/Work/Packages/scGSVA/results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Prepare data and gene sets
# =============================================================================
cat("=== Loading data ===\n")
data(pbmcs)
hsko <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")
expr_mat <- as.matrix(Seurat::GetAssayData(pbmcs, assay = "RNA", layer = "counts"))
cat(sprintf("Expression matrix: %d genes x %d cells\n", nrow(expr_mat), ncol(expr_mat)))

annot_df <- as.data.frame(hsko)
gene_sets_full <- split(annot_df$GeneID, annot_df$Annot)
gene_sets_full <- lapply(gene_sets_full, unique)

# Intersected gene sets (genes present in expr_mat)
gene_sets_int <- lapply(gene_sets_full, function(gs) intersect(gs, rownames(expr_mat)))
gene_sets_int <- gene_sets_int[vapply(gene_sets_int, length, integer(1)) >= 5]
cat(sprintf("Gene sets (intersected, size >= 5): %d\n", length(gene_sets_int)))

cells <- colnames(expr_mat)
n_genes <- nrow(expr_mat)

# =============================================================================
# 2. Run ALL methods
# =============================================================================
cat("\n=== Running all methods ===\n")

# --- scGSVA methods (custom C++ backend) ---
cat("scGSVA ssGSEA... ")
t_sc_ssgsea <- system.time(res_sc_ssgsea <- scgsva(pbmcs, hsko, method = "ssgsea", cores = 1, use.original = FALSE))
cat(round(t_sc_ssgsea["elapsed"], 3), "s\n")

cat("scGSVA GSVA... ")
t_sc_gsva <- system.time(res_sc_gsva <- scgsva(pbmcs, hsko, method = "gsva", cores = 1, use.original = FALSE))
cat(round(t_sc_gsva["elapsed"], 3), "s\n")

cat("scGSVA PLAGE... ")
t_sc_plage <- system.time(res_sc_plage <- scgsva(pbmcs, hsko, method = "plage", cores = 1, use.original = FALSE))
cat(round(t_sc_plage["elapsed"], 3), "s\n")

cat("scGSVA z-score... ")
t_sc_zscore <- system.time(res_sc_zscore <- scgsva(pbmcs, hsko, method = "zscore", cores = 1, use.original = FALSE))
cat(round(t_sc_zscore["elapsed"], 3), "s\n")

# --- Original GSVA package methods ---
cat("GSVA(orig) ssGSEA... ")
t_orig_ssgsea <- system.time(res_orig_ssgsea <- scgsva(pbmcs, hsko, method = "ssgsea", cores = 1, use.original = TRUE))
cat(round(t_orig_ssgsea["elapsed"], 3), "s\n")

cat("GSVA(orig) GSVA... ")
t_orig_gsva <- system.time(res_orig_gsva <- scgsva(pbmcs, hsko, method = "gsva", cores = 1, use.original = TRUE))
cat(round(t_orig_gsva["elapsed"], 3), "s\n")

cat("GSVA(orig) PLAGE... ")
t_orig_plage <- system.time(res_orig_plage <- scgsva(pbmcs, hsko, method = "plage", cores = 1, use.original = TRUE))
cat(round(t_orig_plage["elapsed"], 3), "s\n")

cat("GSVA(orig) z-score... ")
t_orig_zscore <- system.time(res_orig_zscore <- scgsva(pbmcs, hsko, method = "zscore", cores = 1, use.original = TRUE))
cat(round(t_orig_zscore["elapsed"], 3), "s\n")

# --- AUCell ---
cat("AUCell... ")
t_aucell <- system.time({
    rnk_aucell <- AUCell_buildRankings(expr_mat, plotStats = FALSE, verbose = FALSE)
    auc_res <- AUCell_calcAUC(gene_sets_int, rnk_aucell, verbose = FALSE)
    aucell_mat <- getAUC(auc_res)  # pathways x cells
})
cat(round(t_aucell["elapsed"], 3), "s\n")

# --- UCell ---
cat("UCell... ")
t_ucell <- system.time({
    ucell_df <- ScoreSignatures_UCell(expr_mat, features = gene_sets_int,
                                       ncores = 1, maxRank = n_genes)
})
cat(round(t_ucell["elapsed"], 3), "s\n")

# --- singscore ---
cat("singscore... ")
t_singscore <- system.time({
    se <- SummarizedExperiment(assays = list(counts = expr_mat))
    ranked_se <- rankGenes(se)
    ss_list <- lapply(names(gene_sets_int), function(gs_name) {
        gs_genes <- gene_sets_int[[gs_name]]
        if (length(gs_genes) < 2) return(NULL)
        tryCatch({
            res <- simpleScore(ranked_se, upSet = gs_genes)
            score <- res$TotalScore
            names(score) <- rownames(res)
            score
        }, error = function(e) NULL)
    })
    names(ss_list) <- names(gene_sets_int)
    ss_list <- ss_list[!vapply(ss_list, is.null, logical(1))]
    singscore_mat <- do.call(rbind, ss_list)  # pathways x cells
})
cat(round(t_singscore["elapsed"], 3), "s\n")

# =============================================================================
# 3. Build unified score extraction
# =============================================================================

# scGSVA results: cells x pathways data.frame with dot-separated names
# Canonical names use spaces; scGSVA uses dots
scgsva_dotnames <- colnames(res_sc_ssgsea@gsva)
canonical_names <- gsub("\\.", " ", scgsva_dotnames)
dot_lookup <- setNames(scgsva_dotnames, canonical_names)

# Available pathways per tool (canonical space-based names)
avail <- list(
    scGSVA_ssGSEA  = canonical_names,
    scGSVA_GSVA    = gsub("\\.", " ", colnames(res_sc_gsva@gsva)),
    scGSVA_PLAGE   = gsub("\\.", " ", colnames(res_sc_plage@gsva)),
    scGSVA_zscore  = gsub("\\.", " ", colnames(res_sc_zscore@gsva)),
    GSVA_ssGSEA    = gsub("\\.", " ", colnames(res_orig_ssgsea@gsva)),
    GSVA_GSVA      = gsub("\\.", " ", colnames(res_orig_gsva@gsva)),
    GSVA_PLAGE     = gsub("\\.", " ", colnames(res_orig_plage@gsva)),
    GSVA_zscore    = gsub("\\.", " ", colnames(res_orig_zscore@gsva)),
    AUCell         = rownames(aucell_mat),
    UCell          = sub("_UCell$", "", colnames(ucell_df)),
    singscore      = rownames(singscore_mat)
)

# Function to extract a pathways x cells matrix for a given tool + canonical path names
extract_scores <- function(tool, paths) {
    dot_paths <- gsub(" ", ".", paths)
    switch(tool,
        "scGSVA_ssGSEA"  = as.matrix(t(res_sc_ssgsea@gsva[cells, dot_paths, drop = FALSE])),
        "scGSVA_GSVA"    = as.matrix(t(res_sc_gsva@gsva[cells, dot_paths, drop = FALSE])),
        "scGSVA_PLAGE"   = as.matrix(t(res_sc_plage@gsva[cells, dot_paths, drop = FALSE])),
        "scGSVA_zscore"  = as.matrix(t(res_sc_zscore@gsva[cells, dot_paths, drop = FALSE])),
        "GSVA_ssGSEA"    = as.matrix(t(res_orig_ssgsea@gsva[cells, dot_paths, drop = FALSE])),
        "GSVA_GSVA"      = as.matrix(t(res_orig_gsva@gsva[cells, dot_paths, drop = FALSE])),
        "GSVA_PLAGE"     = as.matrix(t(res_orig_plage@gsva[cells, dot_paths, drop = FALSE])),
        "GSVA_zscore"    = as.matrix(t(res_orig_zscore@gsva[cells, dot_paths, drop = FALSE])),
        "AUCell"         = as.matrix(aucell_mat[paths, cells, drop = FALSE]),
        "UCell"          = {
            ucols <- paste0(paths, "_UCell")
            t(as.matrix(ucell_df[cells, ucols, drop = FALSE]))
        },
        "singscore"      = singscore_mat[paths, cells, drop = FALSE]
    )
}

# Compute per-pathway correlations between two tools
compute_cors <- function(tool_a, tool_b, paths) {
    mat_a <- extract_scores(tool_a, paths)
    mat_b <- extract_scores(tool_b, paths)
    sapply(seq_len(nrow(mat_a)), function(k) {
        cor(mat_a[k, ], mat_b[k, ], use = "complete.obs")
    })
}

# =============================================================================
# 4. Benchmark table (all tools)
# =============================================================================
cat("\n=== Benchmark Table ===\n")
bench_df <- data.frame(
    Tool = c("scGSVA ssGSEA", "scGSVA GSVA", "scGSVA PLAGE", "scGSVA z-score",
             "GSVA(orig) ssGSEA", "GSVA(orig) GSVA", "GSVA(orig) PLAGE", "GSVA(orig) z-score",
             "AUCell", "UCell", "singscore"),
    Time_sec = round(c(
        t_sc_ssgsea["elapsed"], t_sc_gsva["elapsed"], t_sc_plage["elapsed"], t_sc_zscore["elapsed"],
        t_orig_ssgsea["elapsed"], t_orig_gsva["elapsed"], t_orig_plage["elapsed"], t_orig_zscore["elapsed"],
        t_aucell["elapsed"], t_ucell["elapsed"], t_singscore["elapsed"]
    ), 3),
    Category = c(rep("scGSVA (C++)", 4), rep("GSVA (original)", 4), "AUCell", "UCell", "singscore"),
    stringsAsFactors = FALSE
)
print(bench_df)
write.csv(bench_df, file.path(results_dir, "benchmark_all_methods.csv"), row.names = FALSE)

# =============================================================================
# 5. Pairwise correlation: each scGSVA method vs everything
# =============================================================================
cat("\n=== Cross-method Correlations ===\n")

all_tools <- names(avail)
external_tools <- c("AUCell", "UCell", "singscore")

# For each scGSVA method, compute correlation with GSVA(orig) counterpart + external tools
scgsva_methods <- c("scGSVA_ssGSEA", "scGSVA_GSVA", "scGSVA_PLAGE", "scGSVA_zscore")
gsva_counterparts <- c("GSVA_ssGSEA", "GSVA_GSVA", "GSVA_PLAGE", "GSVA_zscore")

pairwise_results <- data.frame()

for (i in seq_along(scgsva_methods)) {
    sc_method <- scgsva_methods[i]
    gsva_method <- gsva_counterparts[i]

    # Compare with GSVA counterpart
    common <- intersect(avail[[sc_method]], avail[[gsva_method]])
    if (length(common) > 0) {
        cors <- compute_cors(sc_method, gsva_method, common)
        pairwise_results <- rbind(pairwise_results, data.frame(
            scGSVA_Method = sc_method, Compared_To = gsva_method,
            N_Pathways = length(common),
            Mean_Cor = round(mean(cors, na.rm = TRUE), 6),
            Min_Cor = round(min(cors, na.rm = TRUE), 6),
            Max_Cor = round(max(cors, na.rm = TRUE), 6),
            stringsAsFactors = FALSE
        ))
    }

    # Compare with external tools
    for (ext in external_tools) {
        common <- intersect(avail[[sc_method]], avail[[ext]])
        if (length(common) >= 3) {
            cors <- compute_cors(sc_method, ext, common)
            pairwise_results <- rbind(pairwise_results, data.frame(
                scGSVA_Method = sc_method, Compared_To = ext,
                N_Pathways = length(common),
                Mean_Cor = round(mean(cors, na.rm = TRUE), 6),
                Min_Cor = round(min(cors, na.rm = TRUE), 6),
                Max_Cor = round(max(cors, na.rm = TRUE), 6),
                stringsAsFactors = FALSE
            ))
        }
    }

    # Compare with other GSVA(orig) methods for cross-method correlation
    for (gm in gsva_counterparts) {
        if (gm == gsva_method) next
        common <- intersect(avail[[sc_method]], avail[[gm]])
        if (length(common) >= 3) {
            cors <- compute_cors(sc_method, gm, common)
            pairwise_results <- rbind(pairwise_results, data.frame(
                scGSVA_Method = sc_method, Compared_To = gm,
                N_Pathways = length(common),
                Mean_Cor = round(mean(cors, na.rm = TRUE), 6),
                Min_Cor = round(min(cors, na.rm = TRUE), 6),
                Max_Cor = round(max(cors, na.rm = TRUE), 6),
                stringsAsFactors = FALSE
            ))
        }
    }
}

cat("\nPairwise Correlation Results:\n")
print(pairwise_results)
write.csv(pairwise_results, file.path(results_dir, "pairwise_correlations_all.csv"), row.names = FALSE)

# =============================================================================
# 6. Full cross-tool correlation matrix (all 11 tools)
# =============================================================================
cat("\n=== Full Correlation Matrix ===\n")

# Use pathways common to ALL tools
common_all <- Reduce(intersect, avail)
cat(sprintf("Pathways common to ALL %d tools: %d\n", length(all_tools), length(common_all)))

# If too few common across all, use pairwise approach
if (length(common_all) >= 5) {
    full_cor <- matrix(NA, length(all_tools), length(all_tools),
                       dimnames = list(all_tools, all_tools))
    for (i in seq_along(all_tools)) {
        for (j in seq_along(all_tools)) {
            if (i == j) { full_cor[i, j] <- 1.0; next }
            cors <- compute_cors(all_tools[i], all_tools[j], common_all)
            full_cor[i, j] <- round(mean(cors, na.rm = TRUE), 4)
        }
    }
    cat("\nFull correlation matrix:\n")
    print(full_cor)
    write.csv(full_cor, file.path(results_dir, "full_correlation_matrix_11tools.csv"))
}

# =============================================================================
# 7. FIGURES
# =============================================================================
cat("\n=== Generating Figures ===\n")

# Color scheme
tool_colors <- c(
    "scGSVA ssGSEA" = "#1f77b4", "scGSVA GSVA" = "#2ca02c",
    "scGSVA PLAGE" = "#9467bd", "scGSVA z-score" = "#17becf",
    "GSVA(orig) ssGSEA" = "#ff7f0e", "GSVA(orig) GSVA" = "#d62728",
    "GSVA(orig) PLAGE" = "#e377c2", "GSVA(orig) z-score" = "#bcbd22",
    "AUCell" = "#8c564b", "UCell" = "#7f7f7f", "singscore" = "#393b79"
)

# --- Figure 1: Runtime bar plot (all 11 tools) ---
bench_df$Tool <- factor(bench_df$Tool, levels = bench_df$Tool[order(bench_df$Time_sec)])
p1 <- ggplot(bench_df, aes(x = Tool, y = Time_sec, fill = Category)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = paste0(Time_sec, "s")), hjust = -0.1, size = 2.8) +
    coord_flip() +
    scale_fill_manual(values = c("scGSVA (C++)" = "#2196F3",
                                  "GSVA (original)" = "#FF5722",
                                  "AUCell" = "#8c564b",
                                  "UCell" = "#7f7f7f",
                                  "singscore" = "#393b79")) +
    labs(title = "Runtime: All Gene Set Scoring Methods",
         subtitle = sprintf("%d genes x %d cells, %d KEGG sets",
                            n_genes, length(cells), length(gene_sets_int)),
         x = "", y = "Time (seconds)", fill = "Package") +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "runtime_all_methods.pdf"), p1, width = 10, height = 6)
ggsave(file.path(results_dir, "runtime_all_methods.png"), p1, width = 10, height = 6, dpi = 300)

# --- Figure 2: Speedup chart (scGSVA vs GSVA(orig) for each method) ---
speedup_df <- data.frame(
    Method = c("ssGSEA", "GSVA", "PLAGE", "z-score"),
    scGSVA_Time = c(t_sc_ssgsea["elapsed"], t_sc_gsva["elapsed"],
                    t_sc_plage["elapsed"], t_sc_zscore["elapsed"]),
    GSVA_Time = c(t_orig_ssgsea["elapsed"], t_orig_gsva["elapsed"],
                  t_orig_plage["elapsed"], t_orig_zscore["elapsed"]),
    stringsAsFactors = FALSE
)
speedup_df$Speedup <- round(speedup_df$GSVA_Time / pmax(speedup_df$scGSVA_Time, 0.001), 1)

p2 <- ggplot(speedup_df, aes(x = Method, y = Speedup, fill = Method)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = paste0(Speedup, "x")), vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Speedup: scGSVA (C++) over GSVA (original)",
         x = "Method", y = "Speedup (x times faster)") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "speedup_scgsva_vs_gsva.pdf"), p2, width = 7, height = 5)
ggsave(file.path(results_dir, "speedup_scgsva_vs_gsva.png"), p2, width = 7, height = 5, dpi = 300)
write.csv(speedup_df, file.path(results_dir, "speedup_scgsva_vs_gsva.csv"), row.names = FALSE)

# --- Figure 3: Scatter scGSVA vs GSVA(orig) for each of 4 methods ---
pdf(file.path(results_dir, "scatter_scgsva_vs_gsva_4methods.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (i in seq_along(scgsva_methods)) {
    sc <- scgsva_methods[i]
    gv <- gsva_counterparts[i]
    common <- intersect(avail[[sc]], avail[[gv]])
    mat_sc <- extract_scores(sc, common)
    mat_gv <- extract_scores(gv, common)
    v_sc <- as.vector(mat_sc); v_gv <- as.vector(mat_gv)
    ok <- is.finite(v_sc) & is.finite(v_gv)
    r <- round(cor(v_sc[ok], v_gv[ok]), 6)
    method_name <- sub("scGSVA_", "", sc)
    plot(v_gv[ok], v_sc[ok], pch = 16, cex = 0.3, col = rgb(0.2, 0.4, 0.8, 0.3),
         xlab = paste0("GSVA(orig) ", method_name),
         ylab = paste0("scGSVA ", method_name),
         main = paste0(method_name, " (r = ", r, ")"))
    abline(0, 1, col = "red", lwd = 2, lty = 2)
}
dev.off()

png(file.path(results_dir, "scatter_scgsva_vs_gsva_4methods.png"), width = 1200, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (i in seq_along(scgsva_methods)) {
    sc <- scgsva_methods[i]
    gv <- gsva_counterparts[i]
    common <- intersect(avail[[sc]], avail[[gv]])
    mat_sc <- extract_scores(sc, common)
    mat_gv <- extract_scores(gv, common)
    v_sc <- as.vector(mat_sc); v_gv <- as.vector(mat_gv)
    ok <- is.finite(v_sc) & is.finite(v_gv)
    r <- round(cor(v_sc[ok], v_gv[ok]), 6)
    method_name <- sub("scGSVA_", "", sc)
    plot(v_gv[ok], v_sc[ok], pch = 16, cex = 0.3, col = rgb(0.2, 0.4, 0.8, 0.3),
         xlab = paste0("GSVA(orig) ", method_name),
         ylab = paste0("scGSVA ", method_name),
         main = paste0(method_name, " (r = ", r, ")"))
    abline(0, 1, col = "red", lwd = 2, lty = 2)
}
dev.off()

# --- Figure 4: Each scGSVA method vs AUCell, UCell, singscore (4x3 scatter grid) ---
pdf(file.path(results_dir, "scatter_scgsva_methods_vs_external.pdf"), width = 14, height = 16)
par(mfrow = c(4, 3), mar = c(4, 4, 3, 1))
for (sc in scgsva_methods) {
    method_name <- sub("scGSVA_", "", sc)
    for (ext in external_tools) {
        common <- intersect(avail[[sc]], avail[[ext]])
        if (length(common) < 3) {
            plot.new(); title(main = paste0(method_name, " vs ", ext, " (N/A)")); next
        }
        mat_sc <- extract_scores(sc, common)
        mat_ext <- extract_scores(ext, common)
        v_sc <- as.vector(mat_sc); v_ext <- as.vector(mat_ext)
        ok <- is.finite(v_sc) & is.finite(v_ext)
        r <- round(cor(v_sc[ok], v_ext[ok]), 4)
        plot(v_sc[ok], v_ext[ok], pch = 16, cex = 0.3, col = rgb(0.2, 0.5, 0.7, 0.3),
             xlab = paste0("scGSVA ", method_name),
             ylab = ext,
             main = paste0(method_name, " vs ", ext, " (r=", r, ")"))
        abline(lm(v_ext[ok] ~ v_sc[ok]), col = "red", lwd = 2, lty = 2)
    }
}
dev.off()

png(file.path(results_dir, "scatter_scgsva_methods_vs_external.png"), width = 1400, height = 1600, res = 150)
par(mfrow = c(4, 3), mar = c(4, 4, 3, 1))
for (sc in scgsva_methods) {
    method_name <- sub("scGSVA_", "", sc)
    for (ext in external_tools) {
        common <- intersect(avail[[sc]], avail[[ext]])
        if (length(common) < 3) {
            plot.new(); title(main = paste0(method_name, " vs ", ext, " (N/A)")); next
        }
        mat_sc <- extract_scores(sc, common)
        mat_ext <- extract_scores(ext, common)
        v_sc <- as.vector(mat_sc); v_ext <- as.vector(mat_ext)
        ok <- is.finite(v_sc) & is.finite(v_ext)
        r <- round(cor(v_sc[ok], v_ext[ok]), 4)
        plot(v_sc[ok], v_ext[ok], pch = 16, cex = 0.3, col = rgb(0.2, 0.5, 0.7, 0.3),
             xlab = paste0("scGSVA ", method_name),
             ylab = ext,
             main = paste0(method_name, " vs ", ext, " (r=", r, ")"))
        abline(lm(v_ext[ok] ~ v_sc[ok]), col = "red", lwd = 2, lty = 2)
    }
}
dev.off()

# --- Figure 5: Correlation heatmap (all 11 tools, common pathways) ---
if (length(common_all) >= 5) {
    cor_long <- as.data.frame(as.table(full_cor))
    names(cor_long) <- c("Tool_A", "Tool_B", "Correlation")
    cor_long$Correlation <- as.numeric(cor_long$Correlation)

    # Reorder: scGSVA methods first, then GSVA(orig), then external
    tool_order <- c("scGSVA_ssGSEA", "scGSVA_GSVA", "scGSVA_PLAGE", "scGSVA_zscore",
                    "GSVA_ssGSEA", "GSVA_GSVA", "GSVA_PLAGE", "GSVA_zscore",
                    "AUCell", "UCell", "singscore")
    cor_long$Tool_A <- factor(cor_long$Tool_A, levels = rev(tool_order))
    cor_long$Tool_B <- factor(cor_long$Tool_B, levels = tool_order)

    p5 <- ggplot(cor_long, aes(x = Tool_B, y = Tool_A, fill = Correlation)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = sprintf("%.2f", Correlation)), size = 2.5) +
        scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                             midpoint = 0.5, limits = c(-0.2, 1)) +
        labs(title = "Cross-Tool Score Correlation Matrix (11 Methods)",
             subtitle = sprintf("%d common pathways, mean per-pathway Pearson r", length(common_all)),
             x = "", y = "") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(face = "bold"))

    ggsave(file.path(results_dir, "correlation_heatmap_11tools.pdf"), p5, width = 10, height = 8)
    ggsave(file.path(results_dir, "correlation_heatmap_11tools.png"), p5, width = 10, height = 8, dpi = 300)
}

# --- Figure 6: Per-pathway correlation distribution (scGSVA ssGSEA vs each) ---
pathway_cor_df <- data.frame()
comparisons <- c("GSVA_ssGSEA", "GSVA_GSVA", "AUCell", "UCell", "singscore")
for (comp in comparisons) {
    common <- intersect(avail[["scGSVA_ssGSEA"]], avail[[comp]])
    if (length(common) < 3) next
    cors <- compute_cors("scGSVA_ssGSEA", comp, common)
    pathway_cor_df <- rbind(pathway_cor_df, data.frame(
        Comparison = comp, Correlation = cors, stringsAsFactors = FALSE
    ))
}

p6a <- ggplot(pathway_cor_df, aes(x = Correlation, fill = Comparison)) +
    geom_histogram(bins = 30, alpha = 0.8) +
    facet_wrap(~Comparison, ncol = 1, scales = "free_y") +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Per-Pathway Correlation: scGSVA ssGSEA vs Other Methods",
         x = "Pearson Correlation", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "pathway_cor_dist_ssgsea_vs_all.pdf"), p6a, width = 7, height = 10)
ggsave(file.path(results_dir, "pathway_cor_dist_ssgsea_vs_all.png"), p6a, width = 7, height = 10, dpi = 300)

# --- Figure 6b: Per-pathway correlation distribution (scGSVA GSVA vs each) ---
pathway_cor_df2 <- data.frame()
comparisons2 <- c("GSVA_GSVA", "GSVA_ssGSEA", "AUCell", "UCell", "singscore")
for (comp in comparisons2) {
    common <- intersect(avail[["scGSVA_GSVA"]], avail[[comp]])
    if (length(common) < 3) next
    cors <- compute_cors("scGSVA_GSVA", comp, common)
    pathway_cor_df2 <- rbind(pathway_cor_df2, data.frame(
        Comparison = comp, Correlation = cors, stringsAsFactors = FALSE
    ))
}

p6b <- ggplot(pathway_cor_df2, aes(x = Correlation, fill = Comparison)) +
    geom_histogram(bins = 30, alpha = 0.8) +
    facet_wrap(~Comparison, ncol = 1, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Per-Pathway Correlation: scGSVA GSVA vs Other Methods",
         x = "Pearson Correlation", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "pathway_cor_dist_gsva_vs_all.pdf"), p6b, width = 7, height = 10)
ggsave(file.path(results_dir, "pathway_cor_dist_gsva_vs_all.png"), p6b, width = 7, height = 10, dpi = 300)

# --- Figure 7: Summary bar chart of mean correlations ---
summary_cor <- pairwise_results[, c("scGSVA_Method", "Compared_To", "Mean_Cor")]
summary_cor$scGSVA_Method <- sub("scGSVA_", "scGSVA\n", summary_cor$scGSVA_Method)

p7 <- ggplot(summary_cor, aes(x = Compared_To, y = Mean_Cor, fill = Compared_To)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = round(Mean_Cor, 3)), vjust = -0.3, size = 2.5) +
    facet_wrap(~scGSVA_Method, nrow = 1) +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "Mean Per-Pathway Correlation: Each scGSVA Method vs All Others",
         x = "", y = "Mean Pearson r") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "none",
          plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "summary_correlations_barplot.pdf"), p7, width = 14, height = 5)
ggsave(file.path(results_dir, "summary_correlations_barplot.png"), p7, width = 14, height = 5, dpi = 300)

# --- Figure 8: Feature comparison table (updated with all methods) ---
feature_df <- data.frame(
    Feature = c(
        "Language", "Backend", "Methods available",
        "Seurat integration", "AnnData/Scanpy", "Annotation builder",
        "Visualization", "Statistical testing", "Spatial transcriptomics",
        "Species coverage", "Cross-platform", "GSVA-equivalent (r=1.0)",
        "Parallel execution", "Batch processing"
    ),
    scGSVA = c(
        "R + C++", "Rcpp/C++", "ssGSEA, GSVA, PLAGE, z-score, UCell",
        "Yes (native)", "Yes (pygsva)", "Yes (KEGG/GO/MSigDB, all species)",
        "Yes (6+ types)", "Yes (findPathway, sigPathway)", "Yes",
        "All (AnnotationHub)", "Yes (scGSVA + pygsva)", "Yes (all 4 methods)",
        "Yes (BiocParallel)", "Yes (cell batching)"
    ),
    GSVA = c(
        "R + C", "R + C", "ssGSEA, GSVA, PLAGE, z-score",
        "Manual", "No", "No",
        "No", "No", "No",
        "User provides", "No", "Reference",
        "Yes (BiocParallel)", "No"
    ),
    AUCell = c(
        "R", "R", "AUC-based only",
        "Manual", "pySCENIC (separate)", "No",
        "Threshold plots", "No", "No",
        "User provides", "No", "No (different method)",
        "No", "No"
    ),
    UCell = c(
        "R", "R", "Mann-Whitney U only",
        "Yes (AddModuleScore_UCell)", "No", "No",
        "No", "No", "No",
        "User provides", "No", "No (different method)",
        "Yes (multicore)", "Yes (chunk.size)"
    ),
    singscore = c(
        "R", "R", "Rank-based only",
        "Manual", "No", "No",
        "Yes (landscape)", "Yes (permutation)", "No",
        "User provides", "No", "No (different method)",
        "No", "No"
    ),
    stringsAsFactors = FALSE
)
write.csv(feature_df, file.path(results_dir, "feature_comparison_full.csv"), row.names = FALSE)

# =============================================================================
# 8. Big data comparison (2000 genes x 500 cells)
# =============================================================================
cat("\n=== Big Data Comparison (2000 x 500) ===\n")
set.seed(42)
counts <- as.matrix(Seurat::GetAssayData(pbmcs, assay = "RNA", layer = "counts"))
n_g <- 2000; n_c <- 500
gene_idx <- sample(nrow(counts), n_g, replace = TRUE)
cell_idx <- sample(ncol(counts), n_c, replace = TRUE)
big_counts <- counts[gene_idx, cell_idx]
rownames(big_counts) <- paste0("Gene", seq_len(n_g))
colnames(big_counts) <- paste0("Cell", seq_len(n_c))
big_counts <- big_counts + matrix(rpois(n_g * n_c, lambda = 0.5), nrow = n_g)
big_seu <- Seurat::CreateSeuratObject(counts = big_counts)

orig_genes <- rownames(counts)
gene_map <- orig_genes[gene_idx]
names(gene_map) <- rownames(big_counts)
annot_df_m <- do.call(rbind, lapply(seq_along(gene_map), function(i) {
    matches <- annot_df[annot_df$GeneID == gene_map[i], ]
    if (nrow(matches) > 0) { matches$GeneID <- names(gene_map)[i]; matches } else NULL
}))
big_annot <- new("Annot", species = "human", anntype = "KEGG", keytype = "SYMBOL", annot = annot_df_m)

big_gene_sets <- lapply(gene_sets_full, function(gs) {
    mapped <- names(gene_map)[gene_map %in% gs]
    if (length(mapped) >= 5) mapped else character(0)
})
big_gene_sets <- big_gene_sets[vapply(big_gene_sets, length, integer(1)) >= 5]
big_expr <- as.matrix(Seurat::GetAssayData(big_seu, assay = "RNA", layer = "counts"))

# Run all methods on big data
big_bench <- data.frame()
methods_big <- list(
    list(name = "scGSVA ssGSEA", fn = function() scgsva(big_seu, big_annot, method = "ssgsea", cores = 1, use.original = FALSE, batch = 5000)),
    list(name = "scGSVA GSVA", fn = function() scgsva(big_seu, big_annot, method = "gsva", cores = 1, use.original = FALSE, batch = 5000)),
    list(name = "scGSVA z-score", fn = function() scgsva(big_seu, big_annot, method = "zscore", cores = 1, use.original = FALSE, batch = 5000)),
    list(name = "GSVA(orig) ssGSEA", fn = function() scgsva(big_seu, big_annot, method = "ssgsea", cores = 1, use.original = TRUE, batch = 5000)),
    list(name = "GSVA(orig) GSVA", fn = function() scgsva(big_seu, big_annot, method = "gsva", cores = 1, use.original = TRUE, batch = 5000)),
    list(name = "GSVA(orig) z-score", fn = function() scgsva(big_seu, big_annot, method = "zscore", cores = 1, use.original = TRUE, batch = 5000)),
    list(name = "AUCell", fn = function() {
        rnk <- AUCell_buildRankings(big_expr, plotStats = FALSE, verbose = FALSE)
        AUCell_calcAUC(big_gene_sets, rnk, verbose = FALSE)
    }),
    list(name = "UCell", fn = function() {
        ScoreSignatures_UCell(big_expr, features = big_gene_sets, ncores = 1, maxRank = nrow(big_expr))
    })
)

for (m in methods_big) {
    cat("Big data:", m$name, "... ")
    t <- system.time(m$fn())
    cat(round(t["elapsed"], 3), "s\n")
    big_bench <- rbind(big_bench, data.frame(
        Tool = m$name, Time_sec = round(t["elapsed"], 3),
        stringsAsFactors = FALSE
    ))
}
cat("\nBig Data Benchmark:\n")
print(big_bench)
write.csv(big_bench, file.path(results_dir, "benchmark_big_all_methods.csv"), row.names = FALSE)

# --- Figure 9: Big data runtime ---
big_bench$Tool <- factor(big_bench$Tool, levels = big_bench$Tool[order(big_bench$Time_sec)])
p9 <- ggplot(big_bench, aes(x = Tool, y = Time_sec, fill = Tool)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = paste0(Time_sec, "s")), hjust = -0.1, size = 3) +
    coord_flip() +
    scale_fill_brewer(palette = "Paired") +
    labs(title = sprintf("Runtime: Big Data (%d genes x %d cells)", n_g, n_c),
         x = "", y = "Time (seconds)") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(file.path(results_dir, "runtime_big_all_methods.pdf"), p9, width = 10, height = 5)
ggsave(file.path(results_dir, "runtime_big_all_methods.png"), p9, width = 10, height = 5, dpi = 300)

# =============================================================================
# 9. Summary
# =============================================================================
cat("\n===================================================================\n")
cat("  All results saved to:", results_dir, "\n")
cat("===================================================================\n\n")
cat("Tables:\n")
cat("  benchmark_all_methods.csv        - Small data, all 11 methods\n")
cat("  benchmark_big_all_methods.csv    - Big data, 8 methods\n")
cat("  pairwise_correlations_all.csv    - Each scGSVA method vs all others\n")
cat("  full_correlation_matrix_11tools.csv - 11x11 correlation matrix\n")
cat("  speedup_scgsva_vs_gsva.csv       - Speedup per method\n")
cat("  feature_comparison_full.csv      - Feature comparison table\n")
cat("\nFigures:\n")
cat("  runtime_all_methods.pdf/png              - All 11 tools runtime\n")
cat("  runtime_big_all_methods.pdf/png          - Big data runtime\n")
cat("  speedup_scgsva_vs_gsva.pdf/png           - Speedup bars\n")
cat("  scatter_scgsva_vs_gsva_4methods.pdf/png  - 4 scatter: scGSVA vs GSVA(orig)\n")
cat("  scatter_scgsva_methods_vs_external.pdf/png - 4x3 scatter grid\n")
cat("  correlation_heatmap_11tools.pdf/png      - 11x11 heatmap\n")
cat("  pathway_cor_dist_ssgsea_vs_all.pdf/png   - scGSVA ssGSEA vs others histogram\n")
cat("  pathway_cor_dist_gsva_vs_all.pdf/png     - scGSVA GSVA vs others histogram\n")
cat("  summary_correlations_barplot.pdf/png     - Summary bar chart\n")
