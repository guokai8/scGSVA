
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
#' @return data.frame
#' @export
#' @author Kai Guo
setGeneric("genes", function(object, features,  useTerm = TRUE,
                            with.expr = TRUE)
    standardGeneric("genes"))

#' @title gene method
#' @method genes GSVA
#' @param object GSVA object
#' @param features A vector of features to extract
#' @param useTerm use Term or use id (default: TRUE)
#' @param with.expr extract the expression value or not (default: TRUE)
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
                                with.expr = TRUE) {
    annot <- object@annot
    exp <- object@obj
    input <- exp@assays[["RNA"]]@counts
    if(isTRUE(useTerm)){
        features <- gsub(' ','\\.', features)
        res <- annot[annot[,3]%in%features, ]
    }else{
        res <- annot[annot[,2]%in%features, ]
    }
    if(isTRUE(with.expr)){
        res <- cbind(res,input[res[,1], ])
    }
    return(res)
})
