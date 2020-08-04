#' @importFrom reticulate import py_eval
#' @export
.MACS <- local({
    .MACS <- NULL
    function() {
        if (is.null(.MACS))
            .MACS <<- import("MACS2")
        .MACS
    }
})
