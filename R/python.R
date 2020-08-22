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

.logging <- local({
    .logging <- NULL
    function() {
        if (is.null(.logging)){
            .logging <<- py_run_string("
import logging
def run():
    for h in logging.root.handlers[:]:
        logging.root.removeHandler(h)
        h.close()
")
        }
        .logging
    }
})

.filterdup <- local({
    .filterdup <- NULL
    function() {
        if (is.null(.filterdup))
            .filterdup <<- import("MACS2.filterdup_cmd")
        .filterdup
    }
})

.predictd <- local({
    .predictd <- NULL
    function() {
        if (is.null(.predictd))
            .predictd <<- import("MACS2.predictd_cmd")
        .predictd
    }
})

.pileup <- local({
    .pileup <- NULL
    function() {
        if (is.null(.pileup))
            .pileup <<- import("MACS2.pileup_cmd")
        .pileup
    }
})

.callpeak <- local({
    .callpeak <- NULL
    function() {
        if (is.null(.callpeak))
            .callpeak <<- import("MACS2.callpeak_cmd")
        .callpeak
    }
})

