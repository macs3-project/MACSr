#' @importFrom reticulate py_module_available conda_list use_condaenv
.onLoad <- function(libname, pkgname){
    if(!reticulate::py_module_available("MACS3")){
        if("MACS" %in% reticulate::conda_list()$name){
            warning("Please load MACS env first!")
        }else{
            warning("Please install MACS first!")
        }
    }
}
