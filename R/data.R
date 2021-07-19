#' PBMC 3k

#' @name pbmc
#' @docType data
#' @format A \link[Seurat]{Seurat} object.
#' @keywords datasets
#' @examples
#' data("pbmc")
"pbmc"

#' Example Single Cell RNA-Seq data in SingleCellExperiment Object,
#' subset of 10x public dataset
#' https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k
#' A subset of 390 barcodes and top 200 genes were included in this example.
#' Within 390 barcodes, 195 barcodes are empty droplet, 150 barcodes are cell
#' barcode and 45 barcodes are doublets predicted by scrublet and doubletFinder
#' package. This example only serves as a proof of concept and a tutoriol on how
#' to run the functions in this package. The results should not be used for
#' drawing scientific conclusions.

#' @name sce
#' @docType data
#' @format A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @keywords datasets
#' @examples
#' data("scExample")
"sce"

#' MSigDB gene get Cctegory table
#'
#' A table of gene set categories that can be download from MSigDB. The
#' categories and descriptions can be found here:
#' https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp. The IDs in the
#' first column can be used to retrieve the gene sets for these categories
#' using the \link{importGeneSetsFromMSigDB} function.

#' @name msigb
#' @docType data
#' @format A data.frame.
#' @keywords datasets
#' @examples
#' data("msigdb_table")
"msigdb_table"
