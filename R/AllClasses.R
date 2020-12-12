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
    new("macsList",
        SimpleList(arguments = arguments,
                   outputs = outputs,
                   log = log))
}
