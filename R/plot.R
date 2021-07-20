#' @title Visualize 'Pathways or Terms' on a dimensional reduction plot
#' @description Colors single cells on a dimensional reduction plot according to a
#' 'pathways' (i.e. KEGG, GO terms)
#'
#' @author Kai Guo
#' @export
featurePlot<-function(obj,features,dims,cells,cols,pt.size=0.5,order,repel=TRUE,label = FALSE, label.size = 4){
    features <- features
    seu <- obj@obj

}

#' @title VlnPlot for the pathways
#' @description Draws a violin plot of single cell data (KEGG, GO)
#'
#' @author Kai Guo
#' @export
vlnPlot<-function(obj,features){
    features <- features
    seu <- obj@obj
}
#' @title DotPlot for the specific pathway
#' @description Intuitive way of visualizing how pathway changes across
#' different identity classes (clusters).
#'
#' @author Kai Guo
#' @export
dotPlot<-function(obj,features){
    features <- features
    seu <- obj@obj
}


#' @author Kai Guo
#' @export
ridgePlot<-function(obj,features){
    features <- features
    seu <- obj@obj
}


#' @author Kai Guo
#' @export
Heatmap<-function(obj){
    features <- features
    mat <- as.data.framet(obj@gsva)
}


#' @author Kai Guo
#' @export
ggnet<-function(obj){

}
