#' A small example version of the PBMC dataset
#'
#' 10X Genomics' 3k PBMC dataset
#'
#' @format A Seurat object with the following slots filled
#' \describe{
#'   \item{assays}{
#'   \itemize{Currently only contains one assay ("RNA" - scRNA-seq expression data)
#'   \item{counts - Raw expression data}
#'   \item{data - Normalized expression data}
#'   \item{scale.data - Scaled expression data}
#'   \item{var.features - names of the current features selected as variable}
#'   \item{meta.features - Assay level metadata such as mean and variance}
#'    }}
#'   \item{meta.data}{Cell level metadata}
#'   \item{active.assay}{Current default assay}
#'   \item{active.ident}{Current default idents}
#'   \item{graphs}{Neighbor graphs computed, currently stores the SNN}
#'   \item{reductions}{Dimensional reductions: currently PCA and tSNE}
#'   \item{version}{Seurat version used to create the object}
#'   \item{commands}{Command history}
#' }
#' @source `[pbmc_matrix]` processed by Seurat
#' @examples
#' data("pbmc")
"pbmc"




#' A small example version of the PBMC dataset
#'
#' A subsetted version of 10X Genomics' 3k PBMC dataset
#'
#' @format A Seurat object with the following slots filled
#' \describe{
#'   \item{assays}{
#'   \itemize{Currently only contains one assay ("RNA" - scRNA-seq expression data)
#'   \item{counts - Raw expression data}
#'   \item{data - Normalized expression data}
#'   \item{scale.data - Scaled expression data}
#'   \item{var.features - names of the current features selected as variable}
#'   \item{meta.features - Assay level metadata such as mean and variance}
#'    }}
#'   \item{meta.data}{Cell level metadata}
#'   \item{active.assay}{Current default assay}
#'   \item{active.ident}{Current default idents}
#'   \item{graphs}{Neighbor graphs computed, currently stores the SNN}
#'   \item{reductions}{Dimensional reductions: currently PCA and tSNE}
#'   \item{version}{Seurat version used to create the object}
#'   \item{commands}{Command history}
#' }
#' @source `[pbmc_matrix]` processed by Seurat
#'
"pbmc_small"

#' Small SingleCellExperiment object
#' @name sce_small
#' @docType data
#' @format A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @family data
#' @source \url{"https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/segerstolpe.rds"}
#' @examples
#' data("sce_small")
"sce_small"



