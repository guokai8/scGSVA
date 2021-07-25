#' @title Visualize 'Pathways or Terms' on a dimensional reduction plot
#' @description Colors single cells on a dimensional reduction plot according to a
#' 'pathways' (i.e. KEGG, GO terms)
#' @importFrom viridis scale_color_viridis
#' @importFrom dplyr group_by summarize
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot geom_point facet_wrap aes_string
#' @importFrom ggplot2 scale_color_gradient aes theme_classic
#' @importFrom ggrepel geom_text_repel
#' @importFrom Seurat Embeddings
#' @importFrom tidyr gather
#' @importFrom viridis scale_color_viridis
#' @param object A GSVA objectect or data.frame
#' @param features A vector of features to plot,
#' @param reduction Which dimensionality reduction to use. default("umap")
#' @param color Colors to use for identity class plotting
#' @param group_by Name of one or more metadata columns to group (color) cells by
#' @param label Name of one or more metadata columns to label the cells by
#' @param dims Dimensions to plot, must be a two-length numeric vector
#' specifying x- and y-dimensions
#' @param pt.size Size of the points on the plot
#' @param pt.shape If NULL, all points are circles (default)
#' @param nrow Number of rows
#' @param ncol Number of columns
#' @param basesize base font size, given in pts.
#' @param label.size Sets the size of the labels
#' @param label.color Sets the color of the label text
#' @examples
#' set.seed(123)
#' library(scGSVA)
#' data(pbmc_small)
#' hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
#' res<-scgsva(pbmc_small,hsko)
#' featurePlot(res,features="Wnt.signaling.pathway")
#' @author Kai Guo
#' @export
featurePlot<-function(object, features, reduction = "umap", color = NULL,
                      group_by = NULL, label = NULL,
                      dims = c(1, 2), pt.size = 1, pt.shape = 19,
                      nrow=NULL, ncol = NULL,
                      basesize = 12,
                      label.size = 4, label.color="black"){
    seu <- object@obj
    meta <- seu@meta.data
    red <- Embeddings(seu,reduction)
    gsva <- object@gsva[,features,drop=F]
    gsva <-cbind(red[rownames(gsva),dims],gsva)
    if(!is.null(label)){
        df <- data.frame(text = meta[rownames(gsva),label],
            x = gsva[,1],y=gsva[,2])
        df <- df %>%
            group_by(text) %>%
            summarize(x = median(x), y = median(y))

    }
    if(!is.null(group_by)){
        gsva <- cbind(gsva, meta[rownames(gsva),group_by])
        colnames(gsva)[ncol(gsva)] <- group_by
    }
    if(length(features)>1){
        if(!is.null(group_by)){
            gsva <- gather(gsva, path, val, -1, -2, -ncol(gsva))
        }else{
            gsva <- gather(gsva, path, val, -1, -2)
        }
    }
    xlabel <- colnames(red)[dims[1]]
    ylabel <- colnames(red)[dims[2]]
    p<-ggplot(gsva, aes_string(x = xlabel,
                              y = ylabel
                              ))
    if(length(features)>1){
        p <- p + geom_point(aes(color = val), size = pt.size,
                shape = pt.shape)
    }else{
        p <- p + geom_point(aes_string(color = features),
                size = pt.size, shape = pt.shape)
    }

    if(is.null(color)){
        p <- p + scale_color_viridis()
    }else{
        if(length(color)==2){
            p <- p + scale_color_gradient(low = color[1],high = color[2])
        }else{
            p <- p + scale_color_gradient(low = "white",high = color[1])
        }
    }
    if(!is.null(label)){
        p <- p + geom_text_repel(data = df,aes(x = x, y = y, label = text),
                                 size = label.size, color = label.color)
    }
    if(length(features) > 1){
        if(!is.null(group_by)){
            p <- p + facet_wrap(as.formula(paste0("path",'~',group_by)),ncol=ncol,nrow=nrow)
        }else{
            p <- p + facet_wrap(as.formula(paste0( '.~',"path")),ncol=ncol,nrow=nrow)
        }
    }else{
        if(!is.null(group_by)){
            p <- p + facet_wrap(as.formula(paste('.~',group_by)))
        }
    }
    p <- p + theme_classic(base_size=basesize) + labs(color="NES")
    p
}

#' @title VlnPlot for the pathways
#' @description Draws a violin plot of single cell data (KEGG, GO)
#' @importFrom viridis scale_color_viridis
#' @importFrom dplyr group_by summarize
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot geom_point facet_wrap aes_string
#' @importFrom ggplot2 scale_fill_manual aes theme_classic
#' @importFrom ggplot2 geom_violin geom_jitter
#' @importFrom tidyr gather
#' @importFrom viridis scale_color_viridis
#' @param object A GSVA objectect or data.frame
#' @param features A vector of features to plot
#' @param group_by Name of one or more metadata columns to group (color) cells by
#' @param color Colors to use for identity class plotting
#' @param split.by Factor to split the groups by
#' @param pt.size Size of the points on the plot
#' @param pt.shape If NULL, all points are circles (default)
#' @param nrow Number of rows
#' @param ncol Number of columns
#' @param basesize base font size, given in pts.
#' @examples
#' set.seed(123)
#' library(scGSVA)
#' data(pbmc_small)
#' hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
#' res<-scgsva(pbmc_small,hsko)
#' vlnPlot(res,features="Wnt.signaling.pathway")
#' @author Kai Guo
#' @export
vlnPlot<-function(object, features, group_by = NULL,color = NULL,
                  split.by = NULL,
                  pt.size = 0, pt.shape = 19,
                  nrow=NULL, ncol = NULL,
                  basesize = 12){

    if (inherits(x = object, what = "GSVA")) {
        seu <- object@obj
        gsva <- object@gsva[, features, drop=F]
        meta <- seu@meta.data
        if(is.null(group_by)){
            group_by = "seurat_clusters"
        }
        gsva$group <- meta[rownames(gsva), group_by]
        if(!is.null(split.by)){
            gsva$facet <- meta[rownames(gsva), split.by]
        }
    }else{
        gsva <- object[, features, drop=F]
        gsva <- cbind(gsva,group=group_by)
    }
    if(is.null(color)){
        color <- distcolor[seq_len(length(unique(gsva$group)))]
    }
    if(length(features)>1){
        if(!is.null(split.by)){
            gsva <- gather(gsva, path, val, -group, -facet)
        }else{
            gsva <- gather(gsva, path, val, -group)
        }
    }
    if(length(features) > 1){
        p <- ggplot(gsva,aes_string(x="group",y="val",fill="group"))
    }else{
        p <- ggplot(gsva,aes_string(x="group",y=features,fill="group"))
    }
    p <- p + geom_violin(trim = F)
    if(pt.size>0){
        p <- p + geom_jitter(size = pt.size, shape = pt.shape)
    }
    p <- p + theme_classic(base_size=basesize)+ scale_fill_manual(values=color)+
        theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
        guides(fill = FALSE)+xlab("")
    if(length(features) > 1){
        if(!is.null(split.by)){
            p <- p + facet_wrap(as.formula(paste("facet", '~', "path")),
                            ncol=ncol,nrow=nrow)
        }else{
            p <- p + facet_wrap(as.formula(paste( '.~',"path")),
                            ncol=ncol,nrow=nrow)
        }
        p <- p + ylab("Normalized Enrichment Score")
    }else{
        if(!is.null(split.by)){
            p <- p + facet_wrap(as.formula(paste('.~','facet')),
                            ncol=ncol,nrow=nrow)
        }
        p <- p + ylab(paste0(features,"(NES)"))
    }
    p
}
#' @title DotPlot for the specific pathway
#' @description Intuitive way of visualizing how pathway changes across
#' different identity classes (clusters).
#' @importFrom viridis scale_color_viridis
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot geom_point facet_wrap aes_string theme_classic
#' @importFrom ggplot2 scale_color_gradient aes theme element_text
#' @importFrom tidyr gather
#' @importFrom viridis scale_color_viridis
#' @param object A GSVA objectect or data.frame
#' @param features A vector of features to plot
#' @param color Colors to use for identity class plotting
#' @param split.by Factor to split the groups by
#' @param group_by Name of one or more metadata columns to group (color) cells by
#' @param pt.size Size of the points on the plot
#' @param pt.shape If NULL, all points are circles (default)
#' @param nrow Number of rows
#' @param ncol Number of columns
#' @param basesize base font size, given in pt
#' @examples
#' set.seed(123)
#' library(scGSVA)
#' data(pbmc_small)
#' hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
#' res<-scgsva(pbmc_small,hsko)
#' dotPlot(res,features="Wnt.signaling.pathway")
#' @author Kai Guo
#' @export
dotPlot<-function(object,features,group_by=NULL,split.by=NULL,color=NULL,
                  pt.size = 3, pt.shape = 19,
                  nrow=NULL, ncol = NULL,
                  basesize = 12){
    if (inherits(x = object, what = "GSVA")) {
        seu <- object@obj
        gsva <- object@gsva[, features, drop=F]
        meta <- seu@meta.data
        if(is.null(group_by)){
            group_by = "seurat_clusters"
        }
        gsva$group <- meta[rownames(gsva), group_by]
        if(!is.null(split.by)){
            gsva$facet <- meta[rownames(gsva), split.by]
        }
    }else{
        gsva <- object[, features, drop = F]
        gsva <- cbind(gsva,group = group_by)
    }
    if(length(features)>1){
        if(!is.null(split.by)){
            gsva <- gather(gsva, path, val, -group, -facet)
            gsva <- gsva%>%group_by(path,group,facet)%>%
                summarise(val = mean(val))
        }else{
            gsva <- gather(gsva, path, val, -group)
            gsva <- gsva%>%group_by(path,group)%>%
                summarise(val = mean(val))
        }
    }else{
        colnames(gsva)[1] <- "path"
        gsva <- gsva%>%group_by(group)%>%summarise(val = mean(path))
        gsva$path <- features
    }
    p <- ggplot(gsva,aes_string(x = "group", y = "path"))
    p <- p + geom_point(aes(color = val), size = pt.size,
                            shape = pt.shape)
    if(is.null(color)){
        p <- p + scale_color_viridis()
    }else{
        if(length(color)==2){
            p <- p + scale_color_gradient(low = color[1],high = color[2])
        }else{
            p <- p + scale_color_gradient(low = "white",high = color[1])
        }
    }
    if(!is.null(split.by)){
        p <- p + facet_wrap(as.formula(paste0("facet", '~.')),
                        ncol=ncol,nrow=nrow)
    }
    p <- p + xlab("") + ylab("") + theme_classic(base_size=basesize)+
        labs(color = "Average NES")+
        theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
    p
}

#' @title Generate a ridge plot to examine enrichment distributions
#' @description This function allows to the user to generate the distribution of
#' enrichment across groups with a ridge plot.
#' @importFrom ggridges geom_density_ridges geom_density_ridges2 position_points_jitter
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes_string scale_fill_manual facet_wrap theme_classic
#' @importFrom ggplot2 ylab xlab guides labs
#' @param object A GSVA objectect or data.frame
#' @param features A vector of features to plot
#' @param group_by Name of one or more metadata columns to group (color) cells by
#' or a vector to show the group
#' @param color Colors to use for identity class plotting
#' @param split.by Factor to split the groups by
#' @param rug Adds a rug representation or not
#' @param nrow Number of rows
#' @param ncol Number of columns
#' @param basesize base font size, given in pt
#' @examples
#' set.seed(123)
#' library(scGSVA)
#' data(pbmc_small)
#' hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
#' res<-scgsva(pbmc_small,hsko)
#' ridgePlot(res,features="Wnt.signaling.pathway")
#' @author Kai Guo
#' @export
ridgePlot<-function(object, features, group_by = NULL, color = NULL,
                    split.by = NULL,
                    rug = TRUE,
                    nrow=NULL, ncol = NULL,
                    basesize = 12){
    if (inherits(x = object, what = "GSVA")) {
        seu <- object@obj
        gsva <- object@gsva[, features, drop=F]
        meta <- seu@meta.data
        if(is.null(group_by)){
            group_by = "seurat_clusters"
        }
        gsva$group <- meta[rownames(gsva), group_by]
        if(!is.null(split.by)){
            gsva$facet <- meta[rownames(gsva), split.by]
        }
    }else{
        gsva <- object[, features, drop=F]
        gsva <- cbind(gsva,group=group_by)
        group_by <- ""
    }
    if(is.null(color)){
        color <- distcolor[seq_len(length(unique(gsva$group)))]
    }
    if(length(features)>1){
        if(!is.null(split.by)){
            gsva <- gather(gsva, path, val, -group, -facet)
        }else{
            gsva <- gather(gsva, path, val, -group)
        }
        p<-ggplot(gsva, aes_string(x = "val", y = "group", fill = "group"))
    }else{
        colnames(gsva)[2] <- "group"
        p<-ggplot(gsva, aes_string(x = features, y = "group", fill = "group"))
    }

    if(isTRUE(rug)){
        p <- p + geom_density_ridges(jittered_points = TRUE,
            position = position_points_jitter(width = 0.02, height = 0),
            point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7)

    }else{
        p <- p + geom_density_ridges2(alpha = 0.7)
    }
    p <- p + theme_classic(base_size=basesize) + scale_fill_manual(values=color)+
        ylab(group_by) +
        guides(fill = FALSE)
    if(length(features) > 1){
        if(!is.null(split.by)){
            p <- p + facet_wrap(as.formula(paste0("facet", '~', "path")),
                            ncol=ncol,nrow=nrow)
        }else{
            p <- p + facet_wrap(as.formula(paste0( '.~',"path")),
                            ncol=ncol,nrow=nrow)
        }
        p <- p + xlab("")
    }else{
        if(!is.null(split.by)){
            p <- p + facet_wrap(as.formula(paste0('.~','facet')),
                            ncol=ncol,nrow=nrow)
        }
        p <- p + xlab(paste0(features,"(NES)"))
    }
    p
}

#' Generate heatmap for the enrichment scores
#' @importFrom pheatmap pheatmap
#' @importFrom viridis viridis
#' @param object A GSVA objectect or data.frame
#' @param features A vector of features to plot
#' @param group_by Name of one or more metadata columns to group (color) cells by
#' @param color vector of colors used in heatmap.
#' @param scale character indicating if the values should be centered and
#' scaled in either the row direction or the column direction, or none.
#' Corresponding values are "row", "column" and "none"
#' @param fontsize_row fontsize for rownames (Default: fontsize)
#' @param fontsize_col fontsize for colnames (Default: fontsize)
#' @param border_color color of cell borders on heatmap, use NA if
#' no border should be drawn.
#' @param annotation_col data frame that specifies the annotations shown on
#' top side of the heatmap.
#' @param annotation_colors ist for specifying annotation_row and annotation_col
#' track colors manually. It is possible to define the colors for only some
#' of the features.
#' @param cluster_rows boolean values determining if rows should be clustered
#' or hclust object,
#' @param cluster_cols boolean values determining if columns should be clustered
#' or hclust object.
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @param ... parameters used in pheatmap
#' set.seed(123)
#' library(scGSVA)
#' data(pbmc_small)
#' hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
#' res<-scgsva(pbmc_small,hsko)
#' Heatmap(res,group="groups")
#' @author Kai Guo
#' @export
Heatmap<-function(object, features=NULL, group_by = NULL,
                  color = NULL,scale = "row",
                  fontsize_row = 5,
                  fontsize_col = 5,
                  border_color = "grey60",
                  annotation_col = NA,
                  annotation_colors = NA,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = TRUE, ...){
    if (inherits(x = object, what = "GSVA")) {
        seu <- object@obj
        if(is.null(features)){
            gsva <- object@gsva
            features <- colnames(gsva)
        }else{
            gsva <- object@gsva[, features, drop = FALSE]
        }
        meta <- seu@meta.data
        if(is.null(group_by)){
            group_by = "seurat_clusters"
        }
        gsva <- cbind(gsva,meta[rownames(gsva), group_by, drop = FALSE])
        colnames(gsva)[(length(features)+1):ncol(gsva)]<-group_by

    }else{
        if(is.null(features)){
            gsva <- object
            features <- colnames(gsva)
        }else{
            gsva <- object[, features, drop = FALSE]
        }
        gsva <- cbind(gsva,group_by)
        colnames(gsva)[(length(features)+1):ncol(gsva)]<-group_by
    }
    if(is.null(color)){
        color = viridis(100)
    }
    dat <- t(gsva[,seq_len(length(features))])
    if(is.na(annotation_col)){
        anncol <- gsva[,group_by, drop = FALSE]
    }else{
        anncol <- annotation_col
    }
    if(is.na(annotation_colors)){
        groupl <- apply(gsva[,group_by, drop = FALSE], 2,
                            function(x)length(unique(x)))
        len <- sum(groupl)
        acolors <- c(distcolor,lightcolor)[1:len]
        names(acolors) <- as.vector(apply(gsva[,group_by,drop = FALSE], 2,
                                          function(x)unique(x)))
        anncolors <- split(acolors,rep(names(groupl),times=groupl))

    }else{
        anncolors <- annotation_colors
    }
    pheatmap(dat,scale=scale,color = color,annotation_col = anncol,
             annotation_colors = anncolors,cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,
             fontsize_row = fontsize_row,
             fontsize_col = fontsize_col,
             border_color = border_color,
             show_rownames = show_rownames,show_colnames = show_colnames,...)
}

