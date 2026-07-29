
#' @title GSVA function for single cell data or data.frame with expression value
#' @param obj The count matrix, Seurat, or SingleCellExperiment object.
#' @param annot Annotation object built by \code{buildAnnot()} or \code{buildMSIGDB()}.
#' @param assay Assay to use in GSVA analysis ('RNA','SCT' or 'Spatial' if spatial transcriptomics).
#' @param slot Specific assay data to get or set.
#' @param batch Chunk size for batch processing of large datasets. Default: 1000
#' @param method Enrichment method: \code{"ssgsea"} (default), \code{"gsva"},
#' \code{"plage"}, \code{"zscore"}, \code{"UCell"}, or \code{"plaid"}.
#' @param kcdf Character string denoting the kernel to use during the
#' non-parametric estimation of the cumulative distribution function of
#' expression levels across samples when method="gsva".
#' \code{"Poisson"} (default, for count data) or \code{"Gaussian"} (continuous data).
#' @param abs.ranking Logical flag used only when mx.diff=TRUE. Default: FALSE
#' @param min.sz Minimum size of the resulting gene sets. Default: 1
#' @param max.sz Maximum size of the resulting gene sets. Default: Inf
#' @param mx.diff Logical, two approaches to calculate the enrichment
#' statistic (ES) from the KS random walk statistic. \code{TRUE} (default)
#' uses max difference; \code{FALSE} uses maximum deviation.
#' @param ssgsea.norm Logical, whether to normalize ssGSEA scores. Default: TRUE
#' @param useTerm Logical, use pathway Term names (\code{TRUE}) or IDs (\code{FALSE}). Default: TRUE
#' @param maxRank Maximum number of genes to rank per cell (UCell method only).
#' Gene signatures are automatically filtered to genes present in the expression matrix,
#' so this is auto-adjusted when needed. Default: 1500
#' @param cores Number of parallel workers for computation. Default: 4
#' @param verbose Logical, print progress information. Default: TRUE
#' @param sc.keep Logical, keep the full single cell data in the result object. Default: TRUE
#' @param use.original Logical, if \code{TRUE} use the original GSVA Bioconductor package
#' implementation; if \code{FALSE} (default) use the fast custom C++ implementation.
#' @param ... Additional arguments passed to internal scoring functions.
#' @importFrom SingleCellExperiment counts
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assays
#' @importFrom Matrix colSums Matrix
#' @importFrom Seurat as.Seurat GetAssayData
#' @importFrom BiocParallel SerialParam MulticoreParam
#' @importFrom Matrix summary
#' @examples
#' \dontrun{
#' set.seed(123)
#' library(scGSVA)
#' data(pbmc_small)
#' hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
#' res<-scgsva(pbmc_small,hsko)
#' }
#' @author Kai Guo
#' @export
scgsva <- function(obj, annot = NULL, assay = NULL, slot = "counts",
                   batch = 1000,
                   method="ssgsea",kcdf="Poisson",
                   abs.ranking=FALSE,min.sz=1,
                   max.sz=Inf,
                   mx.diff=TRUE,
                   ssgsea.norm=TRUE,
                   useTerm=TRUE,
                   maxRank=1500,
                   cores = 4,
                   verbose=TRUE,
                   sc.keep=TRUE,
                   use.original=FALSE,...) {
    tau=switch(method, gsva=1, ssgsea=0.25, plaid=NA, NA)
    chk_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk_limit)) {
        cores <- min(cores, 2L)
    }
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
       # input <- obj@assays[[assay]]@counts
        input <- GetAssayData(obj,assay = assay, layer = slot)
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
    input <- input[rowSums(input > 0) != 0, ]
    # Select the enrichment engine
    .engine <- if (isTRUE(use.original)) .sgsva_original else .sgsva
    if(method == "UCell"){
        expressed_genes <- rownames(input)
        annotation <- lapply(annotation, function(g) intersect(g, expressed_genes))
        annotation <- annotation[lengths(annotation) > 0]
        out<- suppressWarnings(UCell::ScoreSignatures_UCell(input, features=annotation,
                                                                chunk.size = batch, ncores = cores,
                                                                maxRank = maxRank,
                                                     BPPARAM=SerialParam(progressbar=verbose),...))
        colnames(out)<-gsub(' ','\\.',sub('_UCell','',colnames(out)))
        out<-as.data.frame(out)
    }else{
        if(ncol(obj) > batch){
            split.data <- split_data_matrix(matrix=input, chunk.size=batch)
            out<- lapply(split.data,function(x).engine(input=x,annotation = annotation,
                                                      method=method,kcdf=kcdf,
                     abs.ranking=abs.ranking,
                     min.sz=min.sz,
                     max.sz=max.sz,cores=cores,
                     mx.diff=mx.diff,
                     tau=tau,ssgsea.norm=FALSE,
                     verbose=verbose
        ))
            out <- do.call(rbind,out)
            if(isTRUE(ssgsea.norm)){
                rng <- range(out) ## na.rm increases execution time and memory consumption
            if (any(is.na(rng) | !is.finite(rng))) rng <- range(out, na.rm=TRUE) ## discard always NA values to calculate the
            out <- out[1:nrow(out),, drop=FALSE] / (rng[2] - rng[1])
        }
    }else{
        out<- .engine(input=input,annotation = annotation, method=method,kcdf=kcdf,
                     abs.ranking=abs.ranking,
                     min.sz=min.sz,
                     max.sz=max.sz,cores=cores,
                     mx.diff=mx.diff,
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
    if (!isTRUE(sc.keep)) {
        if (is.null(assay)) assay <- DefaultAssay(obj)
        empty_counts <- Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE)
        empty_counts <- as(as(empty_counts, "generalMatrix"), "CsparseMatrix")

        # Assign the empty dgCMatrix to the counts slot of the specified assay
        obj <- SetAssayData(object = obj, assay = assay, layer = "counts", new.data = empty_counts)
    }
    res<-new("GSVA",
             obj=obj,
             gsva=out,
             annot=annot)
    return(res)
}

# New .sgsva() dispatcher: routes to custom implementations in gsva_custom.R
.sgsva <- function(input,annotation,method="ssgsea",kcdf="Poisson",
                   abs.ranking=FALSE,min.sz=1,
                   max.sz=Inf,
                   cores=1L,
                   mx.diff=TRUE,
                   tau=switch(method, gsva=1, ssgsea=0.25, NA),
                   ssgsea.norm=TRUE,
                   replace_empty = TRUE,
                   verbose=TRUE){

    if (method == "ssgsea") {
        out <- .ssgsea_custom(input, annotation, alpha = tau,
                              normalize = ssgsea.norm,
                              min.sz = min.sz, max.sz = max.sz)
    } else if (method == "gsva") {
        bpparam <- if (cores > 1L) {
            if (.Platform$OS.type == "windows") {
                BiocParallel::SnowParam(cores)
            } else {
                BiocParallel::MulticoreParam(cores)
            }
        } else NULL
        out <- .gsva_custom(input, annotation, kcdf = kcdf, tau = tau,
                            mx.diff = mx.diff, abs.ranking = abs.ranking,
                            min.sz = min.sz, max.sz = max.sz,
                            BPPARAM = bpparam)
    } else if (method == "plage") {
        out <- .plage_custom(input, annotation,
                             min.sz = min.sz, max.sz = max.sz)
    } else if (method == "plaid") {
        out <- .plaid_custom(input, annotation,
                             normalize = TRUE,
                             min.sz = min.sz, max.sz = max.sz)
    } else {
        out <- .zscore_custom(input, annotation,
                              min.sz = min.sz, max.sz = max.sz)
    }

    output <- data.frame(t(out))
    return(output)
}

# ---- Original GSVA-package-based function (backup for comparison) ----
.sgsva_original <- function(input,annotation,method="ssgsea",kcdf="Poisson",
                   abs.ranking=FALSE,min.sz=1,
                   max.sz=Inf,
                   cores=1L,
                   mx.diff=TRUE,
                   tau=switch(method, gsva=1, ssgsea=0.25, NA),
                   ssgsea.norm=TRUE,
                   replace_empty = TRUE,
                   verbose=TRUE){
    if (!requireNamespace("GSVA", quietly = TRUE)) {
        stop("Package 'GSVA' is required when use.original=TRUE. ",
             "Install it with: BiocManager::install('GSVA')")
    }

    if (method == "gsva") {
        param <- GSVA::gsvaParam(
            input,
            annotation,
            minSize = min.sz,
            maxSize = max.sz,
            kcdf = kcdf,
            tau = tau,
            maxDiff = mx.diff,
            absRanking = abs.ranking
        )
    } else if (method == "ssgsea") {
        param <- GSVA::ssgseaParam(
            input,
            annotation,
            normalize = ssgsea.norm,
            alpha = tau,
            minSize = min.sz,
            maxSize = max.sz
        )
    }else if(method == "plage"){
        param <- GSVA::plageParam(
            input,
            annotation,
            minSize = min.sz,
            maxSize = max.sz
        )
    }else{
        param <- GSVA::zscoreParam(
            input,
            annotation,
            minSize = min.sz,
            maxSize = max.sz
        )
    }
    chk_limit_orig <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if(cores>1 && !nzchar(chk_limit_orig)){
        if (.Platform$OS.type=="windows") {
            BPPARAM <- BiocParallel::SnowParam(cores,type="SOCK",progressbar=verbose)
        } else {
            BPPARAM <- BiocParallel::MulticoreParam(cores,progressbar=verbose)
        }
       out <- suppressWarnings(GSVA::gsva(param, BPPARAM = BPPARAM))
    }else{
       out <- suppressWarnings(GSVA::gsva(param, BPPARAM = BiocParallel::SerialParam(progressbar=verbose)))
    }


    # out<- suppressWarnings(gsva(input, annotation, method = method,kcdf = kcdf,tau=tau,
    #                                ssgsea.norm = ssgsea.norm,  parallel.sz = cores,
    #                                BPPARAM = SerialParam(progressbar=verbose)))
    output <- data.frame(t(out))
    #output <- data.frame(t(out), check.names=FALSE)
    return(output)
}

