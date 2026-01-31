
#' genes
#'
#' @name genes
#' @rdname genes-methods
#' @title gene method
#' @method genes GSVA
#' @param object GSVA object
#' @param features A vector of features to extract
#' @param useTerm use Term or use id (default: TRUE)
#' @param with.expr extract the expression value or not (default: TRUE)
#' @param assay Assay to use in GSVA analysis ('RNA','SCT' or 'Spatial' if spatial transcriptomics)
#' @param slot Specific assay data to get or set
#' @return data.frame
#' @export
#' @author Kai Guo
setGeneric("genes", function(object, features,  useTerm = TRUE,
                            with.expr = TRUE,assay = NULL, slot = "counts")
    standardGeneric("genes"))

#' @title gene method
#' @method genes GSVA
#' @param object GSVA object
#' @param features A vector of features to extract
#' @param useTerm use Term or use id (default: TRUE)
#' @param with.expr extract the expression value or not (default: TRUE)
#' @param assay Assay to use in GSVA analysis ('RNA','SCT' or 'Spatial' if spatial transcriptomics)
#' @param slot Specific assay data to get or set
#' @examples
#' \dontrun{
#' data(pbmc_small)
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' res<-scgsva(pbmc_small, hsako)
#' head(genes(res,"Acute myeloid leukemia"))
#' }
#' @return data.frame
#' @export
#' @author Kai Guo
setMethod("genes", signature(object = "GSVA"),
          definition = function(object, features, useTerm = TRUE,
                                with.expr = TRUE,assay = NULL, slot = "counts") {
            annot <- object@annot
            exp <- object@obj
            #input <- exp@assays[["RNA"]]@counts
            if(is.null(assay)) assay <- "RNA"
            input <- GetAssayData(exp,assay = assay, layer = slot)
            if(isTRUE(useTerm)){
                if(grepl('\\.',features)){
                    features <- gsub('\\.',' ', features)
                }
                res <- annot[annot[,3]%in%features, ]
            }else{
                res <- annot[annot[,2]%in%features, ]
            }
            if(isTRUE(with.expr)){
                res <- cbind(res,input[res[,1], ,drop = FALSE])
            }
            return(res)
})

#' @export
setMethod(
    f = "colnames",
    signature = "GSVA",
    function(x){
        gsva <- x@gsva
        return( colnames(gsva) );
    });
#' @export
setMethod(
    f = "rownames",
    signature = "GSVA",
    function(x){
        gsva <- x@gsva
        return( rownames(gsva) );
    })

#' @export
setMethod(
    f = "colnames",
    signature = "Annot",
    function(x){
        annot <- x@annot
        return( colnames(annot) );
    });

#' @export
setMethod(
    f = "rownames",
    signature = "Annot",
    function(x){
        annot <- x@annot
        return( rownames(annot) );
    })

#' topGenes
#'
#' @name topGenes
#' @title Extract top influential genes for a pathway
#' @description Returns genes in a pathway ranked by their contribution to the
#' enrichment score. Genes are sorted by mean expression across cells, which
#' indicates their influence on GSVA/UCell scores.
#' @param object GSVA object
#' @param features Pathway name to extract genes from
#' @param n Number of top genes to return (default: 10, use Inf for all)
#' @param useTerm use Term or use id (default: TRUE)
#' @param group Optional grouping variable from metadata to calculate group-wise statistics
#' @param assay Assay to use ('RNA','SCT' or 'Spatial')
#' @param slot Specific assay data to get or set
#' @return data.frame with gene names, mean expression, percent expressed, and optionally group-wise stats
#' @examples
#' \dontrun{
#' data(pbmc_small)
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' res<-scgsva(pbmc_small, hsako)
#' topGenes(res, "Wnt signaling pathway", n = 10)
#' topGenes(res, "Wnt signaling pathway", group = "groups")
#' }
#' @export
#' @author Kai Guo
topGenes <- function(object, features, n = 10, useTerm = TRUE,
                     group = NULL, assay = NULL, slot = "counts") {
    if (!inherits(object, "GSVA")) {
        stop("object must be a GSVA object")
    }
    annot <- object@annot
    exp <- object@obj
    if(is.null(assay)) assay <- "RNA"
    input <- GetAssayData(exp, assay = assay, layer = slot)

    # Get genes for the pathway
    if(isTRUE(useTerm)){
        if(grepl('\\.', features)){
            features <- gsub('\\.', ' ', features)
        }
        pathway_genes <- annot[annot[,3] %in% features, 1]
    } else {
        pathway_genes <- annot[annot[,2] %in% features, 1]
    }

    # Filter to genes present in expression matrix
    pathway_genes <- pathway_genes[pathway_genes %in% rownames(input)]

    if(length(pathway_genes) == 0) {
        warning("No genes found for the specified pathway")
        return(NULL)
    }

    # Extract expression for pathway genes
    expr_mat <- as.matrix(input[pathway_genes, , drop = FALSE])

    # Calculate statistics
    mean_expr <- rowMeans(expr_mat)
    pct_expressed <- rowSums(expr_mat > 0) / ncol(expr_mat) * 100

    res <- data.frame(
        gene = pathway_genes,
        mean_expr = mean_expr[pathway_genes],
        pct_expressed = pct_expressed[pathway_genes],
        row.names = NULL
    )

    # Add group-wise statistics if requested
    if(!is.null(group)) {
        meta <- exp@meta.data
        groups <- meta[colnames(expr_mat), group]
        unique_groups <- unique(groups)

        for(g in unique_groups) {
            cells_in_group <- which(groups == g)
            group_mean <- rowMeans(expr_mat[, cells_in_group, drop = FALSE])
            res[[paste0("mean_", g)]] <- group_mean[pathway_genes]
        }
    }

    # Sort by mean expression (descending)
    res <- res[order(res$mean_expr, decreasing = TRUE), ]

    # Return top n genes
    if(!is.infinite(n) && n < nrow(res)) {
        res <- res[1:n, ]
    }

    rownames(res) <- NULL
    return(res)
}



