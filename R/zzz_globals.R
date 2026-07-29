#' @import methods
#' @importFrom stats as.formula model.matrix na.omit p.adjust
#' @importFrom utils data install.packages
#' @importFrom graphics text
NULL

# Suppress R CMD check NOTEs for NSE (non-standard evaluation) variables
# used in dplyr/tidyr/ggplot2 pipelines
utils::globalVariables(c(
  # ggplot2 aes variables
  ".data", "val", "path", "group", "facet", "text", "x", "y",
  "pathway", "value", "mean", "se",
  # tidyr/dplyr pipeline variables
  "type", "GeneID", "GOALL", "ONTOLOGYALL",
  # MSigDB variables
  "gs_collection", "gs_subcollection",
  # annotation variables
  "GO.db", "ONTOLOGY", "kegg", "kegg.db", "module",
  "reactomePATHID2EXTID", "reactomePATHNAME2ID",
  # Seurat/SCE variables
  "Images", "SetAssayData", "logcounts<-",
  # internal functions referenced across files
  ".get_kgm_dat",
  # data objects
  "data", "godata"
))
