#' @title Fit linear model for each pathway
#' @importFrom limma lmFit eBayes topTable
#' @param object A GSVA objectect or data.frame
#' @param group Name of one or more metadata columns to group cells to compare
#' or a vector if the obj is a data.frame
#' @param ref reference group
#' @param pvalue pvalue cut off
#' @param padj padj cut off
#' @export
#' @author Kai Guo
#'
findPathway <- function(object, group = NULL, ref = NULL){
    if (inherits(x = object, what = "GSVA")) {
        seu <- object@obj
        gsva <- object@gsva
        meta <- seu@meta.data
        gsva <- gsva[rownames(gsva)%in%rownames(meta), ]
        if(is.null(group)){
            group = "seurat_clusters"
        }
        group <- meta[rownames(gsva), group]
    }else{
        gsva <- object
    }
    gsva <- t(gsva)
    level <- unique(group)
    if(!is.null(ref)){
        level <- c(ref,setdiff(level, ref))
    }
    group <- factor(as.vector(group), levels = level, ordered = F)
    design <- model.matrix(~group)
    colnames(design) <- levels(group)
    fit <- lmFit(gsva, design)
    fit1 <- eBayes(fit)
    res <- topTable(fit1, adjust = 'BH', coef = 2, number = Inf)
    return(res)
}

