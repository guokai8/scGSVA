.onLoad <- function(libname, pkgname) {
    options(stringsAsFactors = FALSE)
    invisible()
}
##code from Seurat
AttachDeps <- function(deps) {
    for (d in deps) {
        if (!paste0('package:', d) %in% search()) {
            packageStartupMessage("Attaching ", d)
            attachNamespace(ns = d)
        }
    }
}

.onAttach <- function(libname, pkgname) {
    AttachDeps(deps = 'SeuratObject')
}


.onAttach <- function(libname, pkgname) {
    AttachDeps(deps = 'SeuratObject')
}
