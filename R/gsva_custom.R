
# Custom implementations of gene set enrichment scoring methods
# These replace the dependency on the GSVA package for ssGSEA, GSVA, PLAGE, z-score
# Original GSVA-package-based function is preserved as .sgsva_original() in enrich.R

# ---- Helper: filter gene sets by size after intersecting with expression matrix ----
.filter_gene_sets <- function(expr, gene_sets, min.sz = 1, max.sz = Inf) {
    gene_sets <- lapply(gene_sets, function(gs) intersect(gs, rownames(expr)))
    sizes <- vapply(gene_sets, length, integer(1))
    gene_sets[sizes >= min.sz & sizes <= max.sz]
}

# ---- Helper: estimate gene-level CDF across samples (with log-odds transform) ----
# Matches the GSVA package's C code in kernel_estimation.c (row_d function)
# For dense matrices, applies log-odds transform: -log((1-p)/p)
# kcdf: "Gaussian" (continuous data), "Poisson" (count data), or "none" (empirical ECDF)
.gene_cdf <- function(expr, kcdf = "Poisson", BPPARAM = NULL) {
    n_genes <- nrow(expr)
    n_samples <- ncol(expr)

    if (kcdf == "Gaussian" || kcdf == "Poisson") {
        cdf_fn <- if (kcdf == "Gaussian") gene_cdf_gaussian_cpp else gene_cdf_poisson_cpp

        # Use BiocParallel for chunk-based parallelization on large datasets
        n_workers <- 1L
        if (!is.null(BPPARAM)) {
            n_workers <- tryCatch(BiocParallel::bpnworkers(BPPARAM), error = function(e) 1L)
        }

        if (n_workers > 1L && n_genes > 100) {
            # Split genes into chunks for parallel processing
            chunk_size <- ceiling(n_genes / n_workers)
            chunks <- split(seq_len(n_genes), ceiling(seq_len(n_genes) / chunk_size))

            chunk_results <- BiocParallel::bplapply(chunks, function(idx) {
                # C++ function processes row_start to row_end (0-based)
                row_start <- min(idx) - 1L  # 0-based
                row_end <- max(idx)          # exclusive end
                cdf_fn(expr, row_start, row_end)
            }, BPPARAM = BPPARAM)

            # Combine chunks
            cdf_matrix <- do.call(rbind, chunk_results)
        } else {
            # Single-threaded: process all genes at once
            cdf_matrix <- cdf_fn(expr)
        }
    } else {
        # Direct ECDF (no kernel, no log-odds) - R is fine for this
        cdf_matrix <- matrix(0, nrow = n_genes, ncol = n_samples)
        for (i in seq_len(n_genes)) {
            x <- expr[i, ]
            cdf_matrix[i, ] <- rank(x, ties.method = "average") / n_samples
        }
    }

    rownames(cdf_matrix) <- rownames(expr)
    colnames(cdf_matrix) <- colnames(expr)
    return(cdf_matrix)
}

# ===========================================================================
# ssGSEA: Single Sample Gene Set Enrichment Analysis
# Reference: Barbie et al., Nature 2009
#
# For each sample:
#   1. Rank genes by expression (ascending, ties averaged, truncated to integer)
#   2. Order genes by decreasing rank; record each gene's position in that list
#   3. Enrichment score = mean position of non-member genes minus the
#      weighted mean position of member genes (weights = rank^alpha).
#      Intuitively, a gene set is enriched when its members cluster toward
#      the top of the ranked list (low positions) relative to non-members.
#
# This closed-form is an independent algebraic derivation of the integral
# sum(P_hit - P_miss) used in classical ssGSEA (Barbie et al., 2009) and
# produces numerically identical results.
# ===========================================================================
.ssgsea_custom <- function(expr, gene_sets, alpha = 0.25, normalize = TRUE,
                           min.sz = 1, max.sz = Inf) {
    gene_sets <- .filter_gene_sets(expr, gene_sets, min.sz, max.sz)
    if (length(gene_sets) == 0) stop("No gene sets passed the size filter")

    gene_names <- rownames(expr)
    # Map gene set names to integer indices (1-based) for C++
    set_indices <- lapply(gene_sets, function(gs) which(gene_names %in% gs))

    # Fully C++ implementation: ranking (with average tie-breaking and
    # stable sort matching R's rank/order) + walk scoring in one call
    es_matrix <- ssgsea_cpp(expr, set_indices, alpha)
    rownames(es_matrix) <- names(gene_sets)
    colnames(es_matrix) <- colnames(expr)

    # Normalize across all gene sets (ssgsea.norm)
    if (normalize) {
        rng <- range(es_matrix, na.rm = TRUE)
        denom <- rng[2] - rng[1]
        if (denom != 0 && is.finite(denom)) {
            es_matrix <- es_matrix / denom
        }
    }

    return(es_matrix)
}

# ===========================================================================
# GSVA: Gene Set Variation Analysis
# Reference: Hänzelmann et al., BMC Bioinformatics 2013
#
# 1. Estimate gene-level CDF across samples (Gaussian KDE / Poisson / empirical)
# 2. For each sample, order genes by CDF values
# 3. Compute KS-like random walk with tau weighting
# 4. mx.diff=TRUE: ES = sign * (max - min); FALSE: max deviation
# ===========================================================================
.gsva_custom <- function(expr, gene_sets, kcdf = "Poisson", tau = 1,
                         mx.diff = TRUE, abs.ranking = FALSE,
                         min.sz = 1, max.sz = Inf, BPPARAM = NULL) {
    # Filter out constant-expression genes (sd=0) to match original GSVA
    # The original GSVA discards these for non-ssGSEA methods before CDF
    sd_genes <- apply(expr, 1, sd)
    keep <- which(sd_genes > 0 & !is.na(sd_genes))
    if (length(keep) < nrow(expr)) {
        expr <- expr[keep, , drop = FALSE]
    }

    gene_sets <- .filter_gene_sets(expr, gene_sets, min.sz, max.sz)
    if (length(gene_sets) == 0) stop("No gene sets passed the size filter")

    gene_names <- rownames(expr)
    # Map gene set names to integer indices (1-based)
    set_indices <- lapply(gene_sets, function(gs) which(gene_names %in% gs))

    # kcdf mapping: "Poisson" -> 1, "Gaussian" -> 2
    kcdf_int <- if (tolower(kcdf) == "poisson") 1L else 2L

    # Combined C++ pipeline: CDF estimation + ranking + GSVA scoring
    # Avoids intermediate R matrix allocations for CDF and rank matrices
    es_matrix <- gsva_pipeline_cpp(expr, set_indices, kcdf_int, tau,
                                   mx.diff, abs.ranking)
    rownames(es_matrix) <- names(gene_sets)
    colnames(es_matrix) <- colnames(expr)

    return(es_matrix)
}

# ===========================================================================
# PLAGE: Pathway Level Analysis of Gene Expression
# Reference: Tomfohr et al., BMC Bioinformatics 2005
#
# 1. Z-score standardize genes across samples
# 2. For each gene set, SVD on submatrix
# 3. Score = first right singular vector (metagene)
# ===========================================================================
.plage_custom <- function(expr, gene_sets, min.sz = 1, max.sz = Inf) {
    gene_sets <- .filter_gene_sets(expr, gene_sets, min.sz, max.sz)
    if (length(gene_sets) == 0) stop("No gene sets passed the size filter")

    n_samples <- ncol(expr)
    n_sets <- length(gene_sets)

    # Z-score genes across samples (center and scale row-wise)
    # Matches GSVA: Z <- t(scale(t(X)))
    z_expr <- t(scale(t(expr)))
    z_expr[is.na(z_expr)] <- 0

    es_matrix <- matrix(0, nrow = n_sets, ncol = n_samples)
    rownames(es_matrix) <- names(gene_sets)
    colnames(es_matrix) <- colnames(expr)

    for (k in seq_len(n_sets)) {
        set_genes <- gene_sets[[k]]
        n_set <- length(set_genes)
        if (n_set == 0) next

        # GSVA PLAGE: rightsingularsvdvectorgset calls svd(Z[gSetIdx, ])
        # For single-gene sets, Z[gSetIdx, ] drops to a vector, so
        # svd(vector)$v returns matrix(1,1,1) and v[,1]=1 (recycled).
        # For multi-gene sets, Z[gSetIdx, ] is a matrix and svd works normally.
        # We replicate this exact behavior:
        sub_expr <- z_expr[set_genes, ]  # drop=TRUE matches original
        svd_result <- tryCatch(svd(sub_expr), error = function(e) NULL)
        if (is.null(svd_result)) next

        es_matrix[k, ] <- svd_result$v[, 1]
    }

    return(es_matrix)
}

# ===========================================================================
# z-score: Combined z-score method
# Reference: Lee et al., PLoS Computational Biology 2008
#
# 1. Z-score standardize genes across samples
# 2. For each gene set: score = sum(z) / sqrt(n)
# ===========================================================================
.zscore_custom <- function(expr, gene_sets, min.sz = 1, max.sz = Inf) {
    gene_sets <- .filter_gene_sets(expr, gene_sets, min.sz, max.sz)
    if (length(gene_sets) == 0) stop("No gene sets passed the size filter")

    n_samples <- ncol(expr)
    n_sets <- length(gene_sets)

    # Z-score genes across samples
    z_expr <- t(scale(t(expr)))
    z_expr[is.na(z_expr)] <- 0

    es_matrix <- matrix(0, nrow = n_sets, ncol = n_samples)
    rownames(es_matrix) <- names(gene_sets)
    colnames(es_matrix) <- colnames(expr)

    for (k in seq_len(n_sets)) {
        set_genes <- gene_sets[[k]]
        n_set <- length(set_genes)
        if (n_set == 0) next

        sub_expr <- z_expr[set_genes, , drop = FALSE]
        # Combined z-score: sum(z) / sqrt(n_set)
        es_matrix[k, ] <- colSums(sub_expr) / sqrt(n_set)
    }

    return(es_matrix)
}
