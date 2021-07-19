##' Class "Annot"
##' This class represents the Annotation information
##' @name Annot-class
##' @aliases Annot-class
##'   summary, Annot-method
##' @docType class
##' @slot species the species of the annotation file
##' @slot anntype the type of the annotation file
##' @slot keytype Gene ID type
##' @slot annot Annotation information data.frame
##' @exportClass Annot
##' @author Kai Guo
##' @keywords classes
setClass("Annot",
         representation = representation(
           species="character",
           anntype="character",
           keytype="character",
           annot="data.frame"
         ))

##' Class "GSVA"
##' This class represents the Annotation information
##' @name GSVA-class
##' @aliases GSVA-class
##'   summary, GSVA-method
##' @docType class
##' @slot species the species of the annotation file
##' @slot anntype the type of the annotation file
##' @slot keytype Gene ID type
##' @slot annot Annotation information data.frame
##' @exportClass Annot
##' @author Kai Guo
##' @keywords classes
setClass("GSVA",
         representation = representation(
             obj="ANY",
             gsva="data.frame",
             annot = "ANY")
         )
