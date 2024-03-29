#' @title Fit linear model for each pathway
#' @importFrom limma lmFit eBayes topTable
#' @param object A GSVA objectect or data.frame
#' @param group Name of one or more metadata columns to group cells to compare
#' or a vector if the obj is a data.frame
#' @param ref reference group
#' @examples
#' set.seed(123)
#' library(scGSVA)
#' data(pbmc_small)
#' hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
#' sc<-scgsva(pbmc_small,hsko)
#' res <- findPathway(sc, group = "groups", ref = "g1")
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
    }else{
        ref <- level[1]
    }
    lev <- expand.grid(setdiff(level,ref),ref)
    levn <- apply(lev, 1, function(x)paste0(x[1],"_vs_",x[2],"@@"))
    group <- factor(as.vector(group), levels = level, ordered = F)
    design <- model.matrix(~group)
    colnames(design) <- c(ref,levn)
    fit <- lmFit(gsva, design)
    fit1 <- eBayes(fit)
    if(length(level)>2){
        res <- do.call(rbind,sapply(levn,function(x)topTable(fit1,
                        adjust = 'BH', coef=x, number=Inf),simplify=F))
        res$term <- sub('^\\.','',sub('.*@@','',rownames(res)))
        res$comparision <- sub('@@.*','',rownames(res))
        rownames(res)<-NULL
    }else{
        res <- topTable(fit1, adjust = 'BH', coef = 2, number = Inf)
        res$comparision <- sub('@@.*','',levn)
    }
    return(res)
}
#' Significance testing between groups.
#' @param object A GSVA objectect or data.frame
#' @param group Name of one or more metadata columns to group cells to compare
#' or a vector if the obj is a data.frame
#' @param ref reference group
#' @param method correction method, a character string
#' @examples
#' set.seed(123)
#' library(scGSVA)
#' data(pbmc_small)
#' hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
#' sc<-scgsva(pbmc_small,hsko)
#' res <- sigPathway(sc, group = "groups")
#' @export
#' @author Kai Guo
sigPathway<-function(object, group = NULL, test.use = "wilcox", ref = NULL,
                method = "BH"){
    if (inherits(x = object, what = "GSVA")) {
        seu <- object@obj
        gsva <- object@gsva
        meta <- seu@meta.data
        if(is.null(group)){
            group = "seurat_clusters"
        }
        gsva$group <- meta[rownames(gsva), group]
    }else{
        gsva <- object
        gsva <- cbind(gsva, group = group)
        }
    if(test.use == "t.test" | test.use == "t"){
        res <- .do_ttest(gsva, group = "group", ref = ref, method = method)
    }else if(test.use == "wilcox" | test.use == "w"){
        res <- .do_wilcox(gsva, group = "group", ref = ref, method = method)
    }else if(test.use == "anova" | test.use == "aov" | test.use == "a"){
        res <- .do_aov(gsva, group = "group")
    }else{
        cat("Please specify the test method: t.test, wilcox or anova!\n")
    }
    return(res)
}
