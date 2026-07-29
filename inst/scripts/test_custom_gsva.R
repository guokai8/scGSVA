# Test script: Compare custom GSVA implementation vs original GSVA-package-based
# Run this after devtools::load_all() to verify the custom code produces
# comparable results to the original GSVA package.
#
# Usage:
#   devtools::load_all()
#   source("tests/test_custom_gsva.R")

library(scGSVA)

cat("=== Loading test data ===\n")
data(pbmcs)

cat("=== Building KEGG annotation ===\n")
hsko <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")

cat("\n========================================\n")
cat("  Test 1: ssGSEA (default method)\n")
cat("========================================\n")

cat("Running CUSTOM ssGSEA...\n")
t1 <- system.time(res_custom <- scgsva(pbmcs, hsko, method = "ssgsea", cores = 1, use.original = FALSE))
cat("Custom time:", t1["elapsed"], "sec\n")

cat("Running ORIGINAL ssGSEA...\n")
t2 <- system.time(res_original <- scgsva(pbmcs, hsko, method = "ssgsea", cores = 1, use.original = TRUE))
cat("Original time:", t2["elapsed"], "sec\n")

# Compare dimensions
cat("\nDimensions - Custom:", dim(res_custom@gsva), "Original:", dim(res_original@gsva), "\n")

# Compare column names (pathway names)
common_paths <- intersect(colnames(res_custom@gsva), colnames(res_original@gsva))
cat("Common pathways:", length(common_paths), "\n")

# Correlation between scores
if (length(common_paths) > 0) {
    cors <- sapply(common_paths, function(p) {
        cor(res_custom@gsva[[p]], res_original@gsva[[p]], use = "complete.obs")
    })
    cat("Mean Pearson correlation across pathways:", round(mean(cors, na.rm = TRUE), 4), "\n")
    cat("Min correlation:", round(min(cors, na.rm = TRUE), 4), "\n")
    cat("Pathways with cor < 0.9:", sum(cors < 0.9, na.rm = TRUE), "/", length(cors), "\n")
}

cat("\n========================================\n")
cat("  Test 2: GSVA method\n")
cat("========================================\n")

cat("Running CUSTOM GSVA...\n")
t3 <- system.time(res_gsva_custom <- scgsva(pbmcs, hsko, method = "gsva", cores = 1, use.original = FALSE))
cat("Custom time:", t3["elapsed"], "sec\n")

cat("Running ORIGINAL GSVA...\n")
t4 <- system.time(res_gsva_original <- scgsva(pbmcs, hsko, method = "gsva", cores = 1, use.original = TRUE))
cat("Original time:", t4["elapsed"], "sec\n")

common_paths2 <- intersect(colnames(res_gsva_custom@gsva), colnames(res_gsva_original@gsva))
cat("Common pathways:", length(common_paths2), "\n")

if (length(common_paths2) > 0) {
    cors2 <- sapply(common_paths2, function(p) {
        cor(res_gsva_custom@gsva[[p]], res_gsva_original@gsva[[p]], use = "complete.obs")
    })
    cat("Mean Pearson correlation across pathways:", round(mean(cors2, na.rm = TRUE), 4), "\n")
    cat("Min correlation:", round(min(cors2, na.rm = TRUE), 4), "\n")
    cat("Pathways with cor < 0.9:", sum(cors2 < 0.9, na.rm = TRUE), "/", length(cors2), "\n")
}

cat("\n========================================\n")
cat("  Test 3: PLAGE method\n")
cat("========================================\n")

cat("Running CUSTOM PLAGE...\n")
t5 <- system.time(res_plage_custom <- scgsva(pbmcs, hsko, method = "plage", cores = 1, use.original = FALSE))
cat("Custom time:", t5["elapsed"], "sec\n")

cat("Running ORIGINAL PLAGE...\n")
t6 <- system.time(res_plage_original <- scgsva(pbmcs, hsko, method = "plage", cores = 1, use.original = TRUE))
cat("Original time:", t6["elapsed"], "sec\n")

common_paths3 <- intersect(colnames(res_plage_custom@gsva), colnames(res_plage_original@gsva))
cat("Common pathways:", length(common_paths3), "\n")

if (length(common_paths3) > 0) {
    cors3 <- sapply(common_paths3, function(p) {
        cor(res_plage_custom@gsva[[p]], res_plage_original@gsva[[p]], use = "complete.obs")
    })
    cat("Mean Pearson correlation across pathways:", round(mean(cors3, na.rm = TRUE), 4), "\n")
    cat("Min correlation:", round(min(cors3, na.rm = TRUE), 4), "\n")
}

cat("\n========================================\n")
cat("  Test 4: z-score method\n")
cat("========================================\n")

cat("Running CUSTOM z-score...\n")
t7 <- system.time(res_zs_custom <- scgsva(pbmcs, hsko, method = "zscore", cores = 1, use.original = FALSE))
cat("Custom time:", t7["elapsed"], "sec\n")

cat("Running ORIGINAL z-score...\n")
t8 <- system.time(res_zs_original <- scgsva(pbmcs, hsko, method = "zscore", cores = 1, use.original = TRUE))
cat("Original time:", t8["elapsed"], "sec\n")

common_paths4 <- intersect(colnames(res_zs_custom@gsva), colnames(res_zs_original@gsva))
cat("Common pathways:", length(common_paths4), "\n")

if (length(common_paths4) > 0) {
    cors4 <- sapply(common_paths4, function(p) {
        cor(res_zs_custom@gsva[[p]], res_zs_original@gsva[[p]], use = "complete.obs")
    })
    cat("Mean Pearson correlation across pathways:", round(mean(cors4, na.rm = TRUE), 4), "\n")
    cat("Min correlation:", round(min(cors4, na.rm = TRUE), 4), "\n")
}

cat("\n========================================\n")
cat("  Summary\n")
cat("========================================\n")
cat("All tests completed. Review correlations above.\n")
cat("High correlations (>0.9) indicate custom implementation is consistent with original.\n")
cat("Note: exact numeric match is NOT expected due to different KDE implementations.\n")
