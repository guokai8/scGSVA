##' richGO
##'
##' @name richGO
##' @rdname richGO-methods
##' @title richGO method
#' @param x vector contains gene names or dataframe with DEGs information
#' @param godata GO annotation data
#' @param ontology BP,MF or CC
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
##' @return richResult
##' @examples
#' \dontrun{
#' hsago<-buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#' hsago<-as.data.frame(hsago)
#' gene=sample(unique(hsago$GeneID),1000)
#' res<-richGO(gene,godata = hsago,ontology ="BP")
#' }
##' @export
##' @author Kai Guo
setGeneric("richGO", function(x,godata,ontology="BP",pvalue=0.05,padj=NULL,organism=NULL,keytype=NULL,minSize=2,maxSize=500,
                              keepRich=TRUE,filename=NULL,padj.method="BH",sep=",",...)
  standardGeneric("richGO"))

##' richKEGG
##'
##' @name richKEGG
##' @rdname richKEGG-methods
##' @title richKEGG method
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata GO annotation data
#' @param ontology KEGG
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param bulitin use KEGG bulit in KEGG annotation or not(set FALSE if you want use newest KEGG data)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
##' @return richResult
#' @examples
#' \dontrun{
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' gene=sample(unique(hsako$GeneID),1000)
#' res<-richKEGG(gene,kodata = hsako)
#' }
##' @export
##' @author Kai Guo
setGeneric("richKEGG", function(x,kodata,pvalue=0.05,padj=NULL,organism=NULL,ontology="KEGG",
                                keytype=NULL,minSize=2,maxSize=500,
                                keepRich=TRUE,filename=NULL,padj.method="BH",builtin=TRUE,sep=",",...)
  standardGeneric("richKEGG"))

##' richGSEA
##'
##' @name richGSEA
##' @rdname richGSEA-methods
##' @title richGSEA method
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off method
#' @param padj.method p value adjust method
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes when concatenating
##' @return GSEAResult
#' @examples
#' \dontrun{
#' set.seed(123)
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' }
##' @export
##' @author Kai Guo
setGeneric("richGSEA", function(x,object,keytype="",pvalue=0.05,padj=NULL,minSize=15,ontology="",
                                maxSize=500,padj.method="BH",organism=NULL,table=TRUE,sep=",",...)
  standardGeneric("richGSEA"))
##' enrich
##'
##' @name enrich
##' @rdname enrich-methods
##' @title enrich method
#' @param x vector contains gene names or dataframe with DEGs information
#' @param annot ontology type
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
##' @return richResult
##' @examples
##' \dontrun{
##' library(bioAnno)
##' fromKEGG(species="ath")
##' athgo<-buildOwn(dbname="org.ath.eg.db",anntype="GO")
##' athgo<-as.data.frame(athgo)
##' athgo<-na.omit(athgo)
##' gene=sample(unique(athgo$GeneID),1000)
##' res<-enrich(gene,athgo)
##' }
##' @export
##' @author Kai Guo
setGeneric("enrich", function(x,annot,ontology="",pvalue=0.05,padj=NULL,organism=NULL,
                              keytype="",filename=NULL,minSize=2,maxSize=500,
                              keepRich=TRUE,padj.method="BH",sep=",",...)
  standardGeneric("enrich"))


##' ggdot
##'
##' @name ggdot
##' @rdname ggdot-methods
##' @title ggdot method
##' @param object enrichment result or data.frame
##' @param top number of terms you want to display,
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param low low color
##' @param high high color
##' @param alpha transparency alpha
##' @param font.x font of x axis
##' @param font y font of y axis
##' @param fontsize.x fontsize of x axis
##' @param fontsize.y fontsize of y axis
##' @param short automatic short name or not
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param font.size font size for xlim or ylim
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @return ggplot2 object
##' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#'   ggdot(res)
#' }
##' @export
##' @author Kai Guo
setGeneric("ggdot", function(object,top=50,pvalue=0.05,order=FALSE,
                             low="lightpink",high="red",alpha=0.7,
                             font.x="bold",font.y="bold",fontsize.x=10,fontsize.y=10,
                             short=FALSE,
                             padj=NULL,usePadj=TRUE,filename=NULL,width=10,height=8,...)
    standardGeneric("ggdot"))

##' ggbar
##'
##' @name ggbar
##' @rdname ggbar-methods
##' @title ggbar method
##' @param object enrichment result or data.frame
##' @param top number of terms you want to display,
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param low low color
##' @param high high color
##' @param alpha transparency alpha
##' @param font.x font of x axis
##' @param font y font of y axis
##' @param fontsize.x fontsize of x axis
##' @param fontsize.y fontsize of y axis
##' @param short automatic short name or not
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param font.size font size for xlim or ylim
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @param horiz horiz or not
##' @return ggplot2 object
##' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#'   ggbar(res)
#' }
##' @export
##' @author Kai Guo
setGeneric("ggbar", function(object,top=50,pvalue=0.05,padj=NULL,order=FALSE,
                             usePadj=TRUE,fontsize.x=10,fontsize.y=10,short=FALSE,fontsize.text=3,angle=0,filename=NULL,
                             width=10,height=8,horiz=TRUE,...)
    standardGeneric("ggbar"))
##' ggnetplot
##' @rdname ggnetplot-method
##' @title ggnetplot method
##' @param object richResult or dataframe
##' @param top number of terms to show (default: 50)
##' @param pvalue cutoff p value for enrichment result
##' @param padj cutoff p adjust value for enrichment result
##' @param usePadj use adjust p value for the color or not
##' @param low color used for small value
##' @param high color used for large value
##' @param useTerm use terms for nodes (default: TRUE)
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param label.color label color
##' @param label.size label size
##' @param node.shape vector of shape and names of the vector should be the terms (default: 20)
##' @param layout layout method
##' @param savefig save figures or not
##' @param filename output figure name
##' @param width width for output figure
##' @param height height for output figure
##' @param node.alpha alpha-transparency scales
##' @param node.shape shape of the node
##' @param repel use ggrepel text function or not
##' @param segment.size segment size for ggrepel text
##' @param sep character string used to separate the genes when concatenating
##' @return ggplot2 object
##' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#'   ggnetplot(res)
#' }
##' @export
##' @author Kai Guo
setGeneric("ggnetplot",function(object,top=50, pvalue=0.05, padj=NULL,
                                usePadj =TRUE, useTerm=TRUE,low="orange",high="red",
                                writeCyt=FALSE, cytoscapeFile=NULL,
                                label.color = "black", label.size = 2, node.shape=NULL,
                                layout = layout.fruchterman.reingold,savefig=FALSE,filename=NULL,
                                width=7,height=7,node.alpha=0.7,repel=TRUE,segment.size=0.2,sep=",",...)
  standardGeneric("ggnetplot"))


##' ggnetwork
##'
##' @name ggnetwork
##' @rdname ggnetwork-methods
##' @title ggnetwork method
##' @param object richResult,GSEAResult object or dataframe
##' @param gene vector contains gene names or dataframe with DEGs information
##' @param top number of terms to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param weightcut cutoff valule for edge
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param layout layout for the network (layout.fruchterman.reingold)
##' @param low color used for small value
##' @param high color used for large value
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param segment.size size for label segment
##' @param node.alpha alpha-transparency scales
##' @param label.font label font
##' @param label.color label color
##' @param label.size label size
##' @param filename figure output name
##' @param savefig save the figure or not
##' @param width figure width
##' @param height figure height
##' @return ggplot2 object
##' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#'   ggnetwork(res)
#' }
##' @export
setGeneric("ggnetwork", function(object,gene,top = 50, pvalue = 0.05, padj = NULL,usePadj=TRUE,low = "orange",high = "red",
                                 weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = NULL,
                                 label.color = "black", label.size = 2,node.shape=NULL, layout = layout.fruchterman.reingold,savefig=FALSE,
                                 visNet=FALSE,smooth=TRUE,nodeselect=FALSE,edit=FALSE,savehtml=FALSE,filename=NULL,
                                 width=7,height=7,segment.size=0.2,node.alpha=0.7,sep=",",...)
  standardGeneric("ggnetwork"))


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
