
#' @title GSVA function for single cell data or data.frame with expression value
#' @param obj The count matrix, Seurat, or SingleCellExperiment object.
#' @param annot annotation object
#' @param assay Assay to use in GSVA analysis ('RNA','SCT' or 'Spatial' if spatial transcriptomics)
#' @param method to employ in the estimation of gene-set enrichment scores per
#' sample. By default this is set to ssgsea, you can also set it as UCell if you would like use the UCell method
#' @param kcdf Character string denoting the kernel to use during the
#' non-parametric estimation of the cumulative distribution function of
#' expression levels across samples when method="ssgsea".
#' By default, kcdf="Poisson"
#' @param abs.ranking Flag used only when mx.diff=TRUE.
#' @param min.sz Minimum size of the resulting gene sets
#' @param max.sz Maximum size of the resulting gene sets.
#' @param mx.diff Offers two approaches to calculate the enrichment
#' statistic (ES) from the KS random walk statistic.
#' @param ssgsea.norm Logical, set to TRUE (default) with method="ssgsea"
#' runs the SSGSEA method
#' @param useTerm use Term or use id (default: TRUE)
#' @param cores The number of cores to use for parallelization.
#' @param verbose Gives information about each calculation step. Default: FALSE.
#' @importFrom GSVA gsva
#' @importFrom SingleCellExperiment counts
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assays
#' @importFrom Matrix colSums
#' @importFrom Seurat as.Seurat
#' @importFrom BiocParallel SerialParam
#' @importFrom Matrix summary
#' @importFrom UCell ScoreSignatures_UCell
#' @examples
#' set.seed(123)
#' library(scGSVA)
#' data(pbmc_small)
#' hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
#' res<-scgsva(pbmc_small,hsko)
#' @author Kai Guo
#' @export
scgsva <- function(obj, annot = NULL, assay = NULL,
                   batch = 1000,
                   method="ssgsea",kcdf="Poisson",
                   abs.ranking=FALSE,min.sz=1,
                   max.sz=Inf,
                   mx.diff=TRUE,
                   ssgsea.norm=TRUE,
                   useTerm=TRUE,
                   BPPARAM = SnowParam(),
                   cores = 4,
                   verbose=TRUE,...) {
    tau=switch(method, gsva=1, ssgsea=0.25, NA)
    if(is.null(annot)) {
        stop("Please provide annotation object or data.frame")
    } else {
        if(isTRUE(useTerm)){
            annotation <- split(annot[,1],annot[,3])
        }else{
            annotation <- split(annot[,1],annot[,2])
        }
    }
    if (inherits(x = obj, what = "Seurat")) {
        if(is.null(assay)) assay <- "RNA"
        input <- obj@assays[[assay]]@counts
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
    if(method == "UCell"){
        out<- suppressWarnings(ScoreSignatures_UCell(input, features=annotation,
                                                                chunk.size = batch, ncores = cores,
                                                     BPPARAM=SerialParam(progressbar=verbose),...))
        colnames(out)<-gsub(' ','\\.',sub('_UCell','',colnames(out)))
        out<-as.data.frame(out)
    }else{
        if(ncol(obj) > batch){
            split.data <- split.data.matrix(matrix=input, chunk.size=batch)
            out<- lapply(split.data,function(x).sgsva(input=x,annotation = annotation,
                                                      method=method,kcdf=kcdf,
                     abs.ranking=abs.ranking,
                     min.sz=min.sz,
                     max.sz=max.sz,cores=cores,
                     tau=tau,ssgsea.norm=ssgsea.norm,
                     verbose=verbose
        ))
            out <- do.call(rbind,out)
    }else{
        out<- .sgsva(input=input,annotation = annotation, method=method,kcdf=kcdf,
                     abs.ranking=abs.ranking,
                     min.sz=min.sz,
                     max.sz=max.sz,cores=cores,
                     tau=tau,ssgsea.norm=ssgsea.norm,
                     verbose=verbose
        )
    }
    }
    annot <- annot[annot[,1]%in%rownames(input),]
    if(isTRUE(useTerm)){
      annot <- annot[order(annot[,3]),]
      }else{
      annot <- annot[order(annot[,2]),]
      }
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
                   replace_empty = TRUE,
                   verbose=TRUE){
    input <- input[rowSums(input > 0) != 0, ]
    out<- suppressWarnings(gsva(input, annotation, method = method,kcdf = kcdf,tau=tau,
                                   ssgsea.norm = ssgsea.norm,  parallel.sz = cores,
                                   BPPARAM = SerialParam(progressbar=verbose)))
    output <- data.frame(t(out))
    #output <- data.frame(t(out), check.names=FALSE)
    return(output)
}

