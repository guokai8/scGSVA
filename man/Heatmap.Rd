% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{Heatmap}
\alias{Heatmap}
\title{Generate heatmap for the enrichment scores}
\usage{
Heatmap(
  object,
  features = NULL,
  group_by = NULL,
  color = NULL,
  scale = "row",
  average = FALSE,
  order_by = NULL,
  decreasing = FALSE,
  fontsize_row = 5,
  fontsize_col = 5,
  border_color = "grey60",
  annotation_col = NULL,
  annotation_colors = NULL,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  ...
)
}
\arguments{
\item{object}{A GSVA objectect or data.frame}

\item{features}{A vector of features to plot}

\item{group_by}{Name of one or more metadata columns to group (color) cells by}

\item{color}{vector of colors used in heatmap.}

\item{scale}{character indicating if the values should be centered and
scaled in either the row direction or the column direction, or none.
Corresponding values are "row", "column" and "none"}

\item{average}{use the average enrichment score or not (default: FALSE)}

\item{order_by}{a charcter indicate the top of heatmap should order, only
work when cluster_col = FALSE}

\item{decreasing}{logical.  Should the sort order be increasing or decreasing?}

\item{fontsize_row}{fontsize for rownames (Default: fontsize)}

\item{fontsize_col}{fontsize for colnames (Default: fontsize)}

\item{border_color}{color of cell borders on heatmap, use NA if
no border should be drawn.}

\item{annotation_col}{data frame that specifies the annotations shown on
top side of the heatmap.}

\item{annotation_colors}{ist for specifying annotation_row and annotation_col
track colors manually. It is possible to define the colors for only some
of the features.}

\item{cluster_rows}{boolean values determining if rows should be clustered
or hclust object,}

\item{cluster_cols}{boolean values determining if columns should be clustered
or hclust object.}

\item{show_rownames}{boolean specifying if column names are be shown.}

\item{show_colnames}{boolean specifying if column names are be shown.}

\item{...}{parameters used in pheatmap
set.seed(123)
library(scGSVA)
data(pbmc_small)
hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
res<-scgsva(pbmc_small,hsko)
Heatmap(res,group="groups")}
}
\description{
Generate heatmap for the enrichment scores
}
\author{
Kai Guo
}
