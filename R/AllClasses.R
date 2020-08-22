#' macsList
#'
#' @param fun The function executed.
#' @param arguments The arguments used in the function.
#' @param outputs The outputs from the function.
#' @param log The run logs.
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#' @export
setClass("macsList",
         contains = "SimpleList",
         prototype = prototype(
             listData = list(fun = character(),
                             arguments = list(),
                             outputs = list(),
                             log = character())))

macsList <- function(fun = character(), arguments = list(),
                     outputs = list(), log = character()){
    new("macsList",
        SimpleList(fun = fun,
                   arguments = arguments,
                   outputs = outputs,
                   log = log))
}
