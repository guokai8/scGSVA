
##' result generic
##' @param x richResult object
##' @return result return dataframe and print summary
##' @export
result<-function(x){
  UseMethod("result",x)
}
##' detail generic
##' @param x richResult object
##' @return detail return detial for these significant genes
##' @export
detail<-function(x){
  UseMethod("detail",x)
}
#' kappa cluster analysis
##' @param x richResult object or dataframe
##' @param gene (Optional).a vector of gene list
##' @param useTerm to use the term or not (TRUE/FALSE)
##' @param cutoff kappa score threshold for significant dispersion results
##' @param overlap cutoff value of the overlap between two Terms
##' @param minSize minimal number of terms in the cluster
##' @param escore kappa enrichment score cutoff value (default: 3)
#' @examples
#' \dontrun{
#'   hsago<-buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#'   gene=sample(unique(hsago$GeneID),1000)
#'   res<-richGO(gene,godata = hsago,ontology ="BP")
#'   resc<-richCluster(res)
#' }
#' @export
#' @author Kai Guo
setGeneric("richCluster", function(x,gene=NULL,useTerm=FALSE,cutoff=0.5,overlap=0.5,minSize=5,escore=3)
  standardGeneric("richCluster"))
