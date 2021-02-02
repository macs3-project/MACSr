#' macsList
#'
#' @param arguments The arguments used in the function.
#' @param outputs The outputs from the function.
#' @param log The run logs.
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#' @export
setClass("macsList",
         contains = "SimpleList",
         prototype = prototype(
             listData = list(arguments = list(),
                             outputs = list(),
                             log = character())))

macsList <- function(arguments = list(),
                     outputs = list(), log = character()){
    if(length(log) > 0){
        log <- unlist(strsplit(log, split = "\n"))
    }
    new("macsList",
        SimpleList(arguments = arguments,
                   outputs = outputs,
                   log = log))
}

setMethod("show", "macsList", function(object){
    cat("macsList class\n")
    cat("$outputs:\n", paste(object$outputs, collapse = "\n "), "\n")
    cat("$arguments:", paste0(names(object$arguments[-1]),
                              collapse = ", "), "\n")
    cat("$log:\n ")
    cat(paste
    (head(object$log, 5), collapse = "\n "))
    if(length(object$log) > 5){
        cat("\n...\n")
    }
})
