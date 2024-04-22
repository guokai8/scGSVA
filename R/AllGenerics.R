
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



