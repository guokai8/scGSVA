##' @method as.data.frame Annot
##' @export
as.data.frame.Annot<-function(x,...){
  as.data.frame(x@annot)
}
##' @method row.names Annot
##' @export
row.names.Annot<-function(x,...){
  row.names(x@annot)
}

##' @method names Annot
##' @export
names.Annot<-function(x,...){
  names(x@annot)
}
##' @importFrom utils head
##' @method head Annot
##' @export
head.Annot<-function(x,n=6L,...){
  cat("=== species is:",x@species,"and Annotation is",x@anntype," keytype is",x@keytype,"===\n")
  head(x@annot,n,...)
}

##' @importFrom utils tail
##' @method tail Annot
##' @export
tail.Annot<-function(x,n=6L,...){
  cat("=== species is:",x@species,"and Annotation is",x@anntype," keytype is",x@keytype,"===\n")
  tail(x@annot,n,...)
}

##' @method dim Annot
##' @export
dim.Annot <- function(x) {
  dim(x@annot)
}
##' @method as.data.frame GSVA
##' @export
as.data.frame.GSVA<-function(x,...){
  as.data.frame(x@gsva)
}
##' @method row.names GSVA
##' @export
row.names.GSVA<-function(x,...){
  row.names(x@gsva)
}

##' @method names GSVA
##' @export
names.GSVA<-function(x,...){
  names(x@gsva)
}
##' @importFrom utils head
##' @method head GSVA
##' @export
head.GSVA<-function(x,n=6L,...){
  head(x@gsva,n,...)
}

##' @importFrom utils tail
##' @method tail GSVA
##' @export
tail.GSVA<-function(x,n=6L,...){
  tail(x@gsva,n,...)
}

##' @method dim GSVA
##' @export
dim.GSVA <- function(x) {
  dim(x@gsva)
}
##' @method [ Annot
##' @export
`[.Annot` <- function(x, i, j) {
  x@annot[i,j]
}

##' @method $ Annot
##' @export
`$.Annot` <-  function(x, name) {
  x@annot[, name]
}

##' @method $ Annot
##' @export
`$.GSVA` <-  function(x, name) {
  x@gsva[, name]
}

##' @method [ GSVA
##' @export
`[.GSVA` <- function(x, i, j) {
  x@gsva[i,j]
}

##' @method [ GSVA
##' @export
`[<-.GSVA` <- function(x, i, j,..., value) {
    x@gsva[i,j] <- value
    return(x)
}

##' @method $ GSVA
##' @export
`$<-.GSVA` <-  function(x, name,..., value) {
  x@gsva[,name] <- value
  return(x)
}
#' @title subset the GSVA object
##' @method subset GSVA
##' @param x a GSVA object
##' @param ... arguments used in subset.Seurat
##' @export
subset.GSVA <-  function(x, ...) {
 seu <- x@obj
 slot(object = x, name = "obj") <- subset(seu,...)
 gsva <- x@gsva[rownames(x@obj@meta.data),]
 return(new("GSVA",obj = x@obj, gsva = gsva, annot = x@annot))
 #return(x)
}


##' @importFrom AnnotationDbi keys
.get_go_dat<-function(ont="BP"){
  require(GO.db)
  key <- keys(GO.db)
  suppressMessages(go_dat <- AnnotationDbi::select(GO.db, keys = key,
                          columns = c("TERM","ONTOLOGY"),keytype = "GOID"))
  if(ont=="BP") res<-as.data.frame(subset(go_dat, ONTOLOGY == "BP"))
  if(ont=="CC") res<-as.data.frame(subset(go_dat, ONTOLOGY == "CC"))
  if(ont=="MF") res<-as.data.frame(subset(go_dat, ONTOLOGY == "MF"))
  rownames(res) <- res[,1]
  res<-res[, 2, drop = FALSE]
  colnames(res) <- "annotation"
  return(res)
}
##' @importFrom KEGGREST keggList
.get_kg_dat<-function(builtin=TRUE){
  if(isTRUE(builtin)){
    data(kegg)
    return(kegg.db)
  }else{
    pathway<-cbind(keggList('pathway'))
    rownames(pathway)<-sub('path:map','',rownames(pathway))
    colnames(pathway)<-"annotation"
    pathway<-as.data.frame(pathway)
    pathway$annotation<-as.vector(pathway$annotation)
    return(pathway)
  }
}
##' @importFrom KEGGREST keggList
##'
.get_kgm.data <- function(){
  module <-  cbind(keggList('module'))
  rownames(module)<-sub('md:','',rownames(module))
  colnames(module)<-"annotation"
  module<-as.data.frame(module)
  module$annotation<-as.vector(module$annotation)
  return(module)
}

##' build annotaion for kegg
##' @param ontype GO or KEGG
##' @author Kai Guo
.getann<-function(ontype="GO"){
  if(ontype == "GO"){
    res <- rbind(.get_go_dat("BP"),.get_go_dat("MF"),.get_go_dat("CC"))
  }
  if(ontype == "KEGG"){
    res<-.get_kg_dat(builtin = FALSE)
  }
  if(ontype == "Module"){
    res <-.get_kgm_dat()
  }
  return(res)
}


#' Convert ID between ENTREZID to SYMBOL or other type ID based on bioconductor
#' annotation package
#' @title Convert id between different type
#' @importFrom AnnotationDbi mapIds
#' @param species species of gene
#' @param keys input gene id
#' @param fkeytype input gene id type
#' @param tkeytype output gene id type
#' @examples
#' \dontrun{
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' hsako<-as.data.frame(hsako)
#' gene=sample(unique(hsako$GeneID),1000)
#' id<-idconvert(species="human",fkeytype="SYMBOL",tkeytype="ENTREZID")
#' }
#' @return vector
#' @export
#' @author Kai Guo
idconvert<-function(species,keys,fkeytype,tkeytype){
  dbname <- .getdbname(species);
  suppressMessages(require(dbname, character.only = TRUE))
  dbname <- eval(parse(text = dbname))
  unlist(mapIds(dbname, keys = as.vector(keys),
         column = tkeytype,
         keytype = fkeytype,
         multiVals = "first"))
}
.getdbname<-function(species = "human"){
  dbname <- .getdb(species = species);
  if(is.null(dbname)){
    cat("You must check if your request database is avaliable by using showData,
        If not you could make your database by using makeOwnGO and makeOwnKO
        and give a user defined database\n")
    stop("databse error!")
  }
  return(dbname)
}
.getdb <- function(species = "human"){
  species <- tryCatch(match.arg(species,c("anopheles", "arabidopsis",
            "bovine", "celegans", "canine", "fly", "zebrafish",
            "ecoli", "ecsakai", "chicken", "human", "mouse",
            "rhesus","malaria","chipm","rat",
            "toxoplasma", "streptomyces", "pig", "yeast", "xenopus", "warm")),
            error=function(cond){return("unsupported")})
  if (species == "anopheles") {
    dbname <- "org.Ag.eg.db"
  } else if (species == "arabidopsis") {
    dbname <- "org.At.tair.db"
  } else if (species == "bovine") {
    dbname <- "org.Bt.eg.db"
  } else if (species == "canine") {
    dbname <- "org.Cf.eg.db"
  } else if (species == "worm" || species == "celegans") {
    dbname <- "org.Ce.eg.db"
  } else if (species == "chicken") {
    dbname <- "org.Gg.eg.db"
  } else if (species == "ecolik12") {
    dbname <- "org.EcK12.eg.db"
  } else if (species == "ecsakai") {
    dbname <- "org.EcSakai.eg.db"
  } else if (species == "fly") {
    dbname <- "org.Dm.eg.db"
  } else if (species == "human") {
    dbname <- "org.Hs.eg.db"
  } else if (species == "malaria") {
    dbname <- "org.Pf.plasmo.db"
  } else if (species == "chipm") {
    dbname <- "org.Pt.eg.db"
  }else if (species == "mouse") {
    dbname <- "org.Mm.eg.db"
  } else if (species == "pig") {
    dbname <- "org.Ss.eg.db"
  } else if (species == "rat") {
    dbname <- "org.Rn.eg.db"
  } else if (species == "rhesus") {
    dbname <- "org.Mmu.eg.db"
  } else if (species == "xenopus") {
    dbname <- "org.Xl.eg.db"
  } else if (species == "yeast") {
    dbname <- "org.Sc.sgd.db"
  } else if (species == "streptomyces") {
    dbname <- "org.Sco.eg.db"
  } else if (species == "zebrafish") {
    dbname <- "org.Dr.eg.db"
  } else if (species == "toxoplasma"){
    dbname<- "org.Tgondii.eg.db"
  } else {
    dbname <- NULL
  }
  return(dbname)
}
.getspeices<-function(species="human"){
  species=tryCatch(match.arg(species,c("anopheles","arabidopsis","bovine","celegans","canine","fly","zebrafish",
                                       "ecoli","ecsakai","chicken","human","mouse","rhesus","malaria","chipm","rat",
                                       "toxoplasma","sco","pig","yeast","xenopus","warm")),
                   error=function(cond){return("unsupported")})
  if (species == "anopheles") {
    species <- "aga"
  } else if (species == "arabidopsis") {
    species <- "ath"
  } else if (species == "bovine") {
    species <- "bta"
  } else if (species == "canine") {
    species <- "cfa"
  } else if (species == "chicken") {
    species <- "gga"
  } else if (species == "chipm") {
    species <- "ptr"
  } else if (species == "ecolik12") {
    species <- "eco"
  } else if (species == "ecsakai") {
    species <- "ecs"
  } else if (species == "fly") {
    species <- "dme"
  } else if (species == "human") {
    species <- "hsa"
  } else if (species == "malaria") {
    species <- "pfa"
  } else if (species == "mouse") {
    species <- "mmu"
  } else if (species == "pig") {
    species <- "ssc"
  } else if (species == "rat") {
    species <- "rno"
  } else if (species == "rhesus") {
    species <- "mcc"
  } else if (species == "worm" || species == "celegans") {
    species <- "cel"
  } else if (species == "xenopus") {
    species <- "xla"
  } else if (species == "yeast") {
    species <- "sce"
  } else if (species =="streptomyces"){
    species <- "sco"
  } else if (species == "zebrafish") {
    species <- "dre"
  } else {
    species <- NULL
  }
  return(species)
}
#' show avaliable data based on bioconductor annotation package
#' @export
#' @author Kai Guo
showData<-function(){
  species=c("anopheles","arabidopsis","bovine","celegans","canine","fly","zebrafish",
            "ecoli","ecsakai","chicken","human","mouse","rhesus","malaria","chipm","rat",
            "toxoplasma","sco","pig","yeast","xenopus")
  dbname=c("org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db","org.Cf.eg.db","org.Dm.eg.db",
           "org.Dr.eg.db","org.EcK12.eg.db","org.EcSakai.eg.db","org.Gg.eg.db","org.Hs.eg.db","org.Mm.eg.db",
           "org.Mmu.eg.db","org.Pf.plasmo.db","org.Pt.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Sco.eg.db",
           "org.Ss.eg.db","org.Tgondii.eg.db","org.Xl.eg.db")
  dbdata<-data.frame(dbname=dbname,species=species)
  dbdata
}

.vec_to_df<-function(x, name){
  dd <- data.frame(names(x), x)
  colnames(dd) <- name
  return(dd)
}

##' msigdb support species
##' @param species with common name
.getmsig<-function(species="human"){
  out<-NULL
  if(species=="human"){
    out<-"Homo sapiens"
  }else if(species=="mouse"){
    out<-"Mus musculus"
  }else if(species=="rat"){
    out<-"Rattus norvegicus"
  }else if(species=="celegans"){
    out<-"Caenorhabditis elegans"
  }else if(species=="fly"){
    out<-"rosophila melanogaster"
  }else if(species=="yeast"){
    out<-"Saccharomyces cerevisiae"
  }else if(species=="bovine"){
    out<-"Bos taurus"
  }else if(species=="canine"){
    out<-"Canis lupus familiaris"
  }else if(species=="pig"){
    out<-"Sus scrofa"
  }else if(species=="chicken"){
    out<-"Gallus gallus"
  }else if(species=="zebrafish"){
    out<-"Danio rerio"
  }else{
    out<-NULL
  }
}
##' Print MSIGDB infomation
##' @export
msigdbinfo <- function() {
  cat("#--------------------------------------------------------------#\n")
  cat("# Molecular Signatures Database                        v6.2.1  #\n")
  cat("#--------------------------------------------------------------#\n")
  cat("# Category | Subcategory # Details ----------------------------#\n")
  cat("# C1               # Positional (326)                          #\n")
  cat("# C2 | CGP         # Chemical and Genetic Perturbations (3433) #\n")
  cat("# C2 | CP          # Canonical Pathways (252)                  #\n")
  cat("# C2 | BIOCARTA # Canonical BIOCARTA (217)                     #\n")
  cat("# C2 | KEGG     # Canonical KEGG (186)                         #\n")
  cat("# C2 | CPREACTOME  # Canonical REACTOME (674)                  #\n")
  cat("# C3 | MIR         # Motif miRNA Targets (221)                 #\n")
  cat("# C3 | TFT         # Motif Transcription Factor Targets (615)  #\n")
  cat("# C4 | CGN         # Cancer Gene Neighborhoods (427)           #\n")
  cat("# C4 | CM          # Cancer Modules (431)                      #\n")
  cat("# C5 | BP          # GO Biological Process (4436)              #\n")
  cat("# C5 | CC          # GO Cellular Component (580)               #\n")
  cat("# C5 | MF          # GO Molecular Function (901)               #\n")
  cat("# C6               # Oncogenic Signatures (189)                #\n")
  cat("# C7               # Immunologic Signatures (4872)             #\n")
  cat("# H                # Hallmark (50)                             #\n")
  cat("#--------------------------------------------------------------#\n")
  cat("# Source: http://software.broadinstitute.org/gsea/msigdb       #\n")
  cat("#--------------------------------------------------------------#\n")
  listspe<-c("human","mouse","rat","celegans","fly","yeast","bovine","canine",
             "pig","chicken","zebrafish")
  cat("# Support species:                                             #\n")
  cat(sort(listspe),"\n")
}
.getrodbname <- function(species){
 # "Schizosaccharomyces pombe","Taeniopygia guttata",,"Mycobacterium tuberculosis"
  spe<-c("Homo sapiens","Dictyostelium discoideum","Plasmodium falciparum",
         "Saccharomyces cerevisiae","Caenorhabditis elegans",
         "Sus scrofa","Bos taurus","Canis familiaris","Mus musculus","Rattus norvegicus",
         "Xenopus tropicalis","Danio rerio",
         "Drosophila melanogaster","Arabidopsis thaliana","Oryza sativa","Gallus gallus")
  names(spe)<-c("human","dictyostelium","malaria","yeast","celegans","pig",
                    "bovine","dog","mouse","rat","xenopus","zebrafish","fly","arabidopsis","rice","chicken")
  return(spe[species])
}
##'
setAs(from = "data.frame", to = "Annot", def = function(from){
  keytype <- character()
  species <- character()
  anntype <- character()
  GeneID <- as.vector(from[,1])
  Term<-as.vector(from[,2])
  Annot=from$Annot
  annot <- data.frame(GeneID, Term, Annot)
  new("Annot",
      species = species,
      anntype = anntype,
      keytype = keytype,
      annot = annot
  )
})

.paste.char<-function(x){
  return(gsub("([^ ]+ [^ ]+ [^ ]+ [^ ]+) ", "\\1\n", x))
}

.clean.char<-function(x){
  return(gsub('\\\n',' ',x))
}

##'
setAs(from = "GSVA", to = "data.frame", def = function(from){
  result <- as.data.frame(from@gsva)
  result
})
##'
setAs(from = "Annot", to = "data.frame", def = function(from){
  result <- as.data.frame(from@annot)
  result
})


.is_inst <- function(pkg) {
  nzchar(system.file(package = pkg))
}
.load_pkg<-function(dbname){
  if (!.is_inst(dbname)){
    if(!.is_inst("BiocManager")){
      install.packages("BiocManager")
    }else{
      BiocManager::install(dbname,update = FALSE)
    }
    suppressMessages(require(dbname, character.only = TRUE,quietly = TRUE))
  }else{
    suppressMessages(require(dbname, character.only = TRUE,quietly = TRUE))
  }
  dbname<-eval(parse(text=dbname))
  return(dbname)
}
#' distinguish colors for making figures
#' @author Kai Guo
#' @export
distcolor<-c("#66A61E","#1B9E77","#E7298A","#7570B3","#E6AB02","#A6761D",
             "#D95F02",
             "#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20",
             "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
             "#A6761D","#D95F02","#66A61E","#1B9E77","#E7298A","#7570B3",
             "#E6AB02",'#e6194b', '#3cb44b', '#4363d8',
             '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c',
             '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8',
             '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075',
             '#808080', '#ffffff', '#000000')

#' light colors for making figures
#' @author Kai Guo
#' @export
lightcolor<-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
              '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
              '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
              '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
              '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
              '#968175','#e6194b', '#3cb44b', '#ffe119', '#4363d8','#f58231', '#911eb4',
              '#46f0f0', '#f032e6', '#bcf60c',
              '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8',
              '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075',
              '#808080'
)


#' do anova test and return results as data.frame
#' @importFrom rstatix anova_test
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @param x data.frame with sample id as the column name, genes or otu as rownames
#' @param group group factor used for comparison
#' @param ... parameters to anova_test
#' @examples
#' {
#' data("ToothGrowth")
#' do_aov(ToothGrowth,group="supp")
#' }
#' @author Kai Guo
.do_aov<-function(x,group,...){
  d<-x[,setdiff(colnames(x),group)]
  d$group<-x[,group]
  d<-d%>%gather(type,val,-group)
  res<-d%>%group_by(type)%>%anova_test(val~group,...)
  return(res)
}

#' do t.test
#' @importFrom rstatix t_test
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @param x data.frame with sample id as the column name, genes or otu as rownames
#' @param group group factor used for comparison
#' @param ref reference group
#' @param method correction method, a character string
#' @param ... parameters to t_test
#' @examples
#' {
#' data("mtcars")
#' do_ttest(mtcars,group="vs")
#' do_ttest(mtcars,group="cyl",ref="4")
#' }
#' @author Kai Guo
.do_ttest<-function(x,group,ref=NULL,method = "BH",...){
  d<-x[,setdiff(colnames(x),group)]
  d$group<-x[,group]
  d<-d%>%gather(type,val,-group)
  res<-d%>%group_by(type)%>%t_test(val~group,ref.group = ref,...)
  res$p.adj<-p.adjust(res$p,method = method)
  return(res)
}

#' do wilcox test
#' @importFrom rstatix wilcox_test
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @param x data.frame with sample id as the column name, genes or otu as rownames
#' @param group group factor used for comparison
#' @param ref reference group
#' @param method correction method, a character string
#' @param ... parameters to wilcox_test
#' @examples
#' {
#' data("mtcars")
#' do_wilcox(mtcars,group="vs")
#' do_wilcox(mtcars,group="cyl",ref="4")
#' }
#' @author Kai Guo
.do_wilcox<-function(x,group,ref=NULL, method = "BH",...){
  d<-x[,setdiff(colnames(x),group)]
  d$group<-x[,group]
  d<-d%>%gather(type,val,-group)
  res<-d%>%group_by(type)%>%wilcox_test(val~group,ref.group = ref,...)
  res$p.adj<-p.adjust(res$p,method=method)
  return(res)
}
