#' @importFrom reticulate import py_run_string
.MACS <- local({
    .MACS <- NULL
    function() {
        if (is.null(.MACS))
            .MACS <<- import("MACS2")
        .MACS
    }
})

.namespace <- local({
    .namespace <- NULL
    function() {
        if (is.null(.namespace)){
            .namespace <<- py_run_string("
class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
")
        }
        .namespace
    }
})

## .namespace <<- py_run_string("
## class Namespace:
##     def __init__(self, **kwargs):
##         self.__dict__.update(kwargs)
## ")

.filterdup <- local({
    .filterdup <- NULL
    function() {
        if (is.null(.filterdup))
            .filterdup <<- import("MACS2.filterdup_cmd")
        .filterdup
    }
})
