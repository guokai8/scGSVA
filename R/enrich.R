
#' @title GSVA function for single cell data or data.frame with expression value
#' @param obj The count matrix, Seurat, or SingleCellExperiment object.
#' @param annotation annotation object
#' @param cores The number of cores to use for parallelization.
#'
#' @importFrom GSVA gsva
#' @importFrom SingleCellExperiment counts
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assays
#' @importFrom Matrix colSums
#' @importFrom Seurat as.Seurat
#' @importFrom BiocParallel SerialParam
#' @importFrom Matrix summary
#' @author Kai Guo
#' @export
scgsva <- function(obj, annot = NULL, cores = 4,
                   method="ssgsea",kcdf="Poisson",
                   abs.ranking=FALSE,min.sz=1,
                   max.sz=Inf,
                   mx.diff=TRUE,
                   ssgsea.norm=TRUE,
                   useTerm=TRUE,
                   verbose=TRUE) {
    tau=switch(method, gsva=1, ssgsea=0.25, NA)
    if(is.null(annot)) {
        stop("Please provide anotation object or data.frame")
    } else {
        if(isTRUE(useTerm)){
            annotation <- split(annot[,1],annot[,3])
        }else{
            annotation <- split(annot[,1],annot[,2])
        }
    }
    if (inherits(x = obj, what = "Seurat")) {
        input <- obj@assays[["RNA"]]@counts
        input<- input[tabulate(summary(input)$i) != 0, , drop = FALSE]
        input <- as.matrix(input)
    } else if (inherits(x = obj, what = "SingleCellExperiment")) {
        input <- counts(obj)
        if(!"logcounts"%in%names(assays(obj))){
            libsizes <- colSums(assay(obj,"counts"))
            size.factors <- libsizes/mean(libsizes)
            logcounts(obj) <- as.matrix(log(t(t(input)/size.factors) + 1))
        }else{
            logcounts(obj) <- as.matrix(logcounts(obj))
        }
        input<- input[tabulate(summary(input)$i) != 0, , drop = FALSE]
        input <- as.matrix(input)
        obj<-as.Seurat(obj)

    } else {
        input <- obj
    }
    out<- .sgsva(input=input,annotation = annotation,method=method,kcdf=kcdf,
                 abs.ranking=abs.ranking,
                 min.sz=min.sz,
                 max.sz=max.sz,cores=cores,
                 tau=tau,ssgsea.norm=ssgsea.norm,
                 verbose=verbose
                 )
    res<-new("GSVA",
             obj=obj,
             gsva=out,
             annot=annot)
    return(res)
}

.sgsva <- function(input,annotation,method="ssgsea",kcdf="Poisson",
                   abs.ranking=FALSE,min.sz=1,
                   max.sz=Inf,
                   cores=1L,
                   mx.diff=TRUE,
                   tau=switch(method, gsva=1, ssgsea=0.25, NA),
                   ssgsea.norm=TRUE,
                   verbose=TRUE){
    input <- input[rowSums(input > 0) != 0, ]
    out<- suppressWarnings(gsva(input, annotation, method = method,kcdf = kcdf,tau=tau,
                                   ssgsea.norm = ssgsea.norm,  parallel.sz = cores,
                                   BPPARAM = SerialParam(progressbar=verbose)))
    output <- data.frame(t(out),check.names = F)
    return(output)
}

