#' @importFrom reticulate import py_run_string
.MACS <- local({
    .MACS <- NULL
    function() {
        if (is.null(.MACS))
            .MACS <<- import("MACS3")
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
            .filterdup <<- import("MACS3.Commands.filterdup_cmd")
        .filterdup
    }
})

.predictd <- local({
    .predictd <- NULL
    function() {
        if (is.null(.predictd))
            .predictd <<- import("MACS3.Commands.predictd_cmd")
        .predictd
    }
})

.pileup <- local({
    .pileup <- NULL
    function() {
        if (is.null(.pileup))
            .pileup <<- import("MACS3.Commands.pileup_cmd")
        .pileup
    }
})

.callpeak <- local({
    .callpeak <- NULL
    function() {
        if (is.null(.callpeak))
            .callpeak <<- import("MACS3.Commands.callpeak_cmd")
        .callpeak
    }
})

.bdgpeakcall <- local({
    .bdgpeakcall <- NULL
    function() {
        if (is.null(.bdgpeakcall))
            .bdgpeakcall <<- import("MACS3.Commands.bdgpeakcall_cmd")
        .bdgpeakcall
    }
})

.bdgbroadcall <- local({
    .bdgbroadcall <- NULL
    function() {
        if (is.null(.bdgbroadcall))
            .bdgbroadcall <<- import("MACS3.Commands.bdgbroadcall_cmd")
        .bdgbroadcall
    }
})

.bdgcmp <- local({
    .bdgcmp <- NULL
    function() {
        if (is.null(.bdgcmp))
            .bdgcmp <<- import("MACS3.Commands.bdgcmp_cmd")
        .bdgcmp
    }
})

.bdgopt <- local({
    .bdgopt <- NULL
    function() {
        if (is.null(.bdgopt))
            .bdgopt <<- import("MACS3.Commands.bdgopt_cmd")
        .bdgopt
    }
})

.cmbreps <- local({
    .cmbreps <- NULL
    function() {
        if (is.null(.cmbreps))
            .cmbreps <<- import("MACS3.Commands.cmbreps_cmd")
        .cmbreps
    }
})

.bdgdiff <- local({
    .bdgdiff <- NULL
    function() {
        if (is.null(.bdgdiff))
            .bdgdiff <<- import("MACS3.Commands.bdgdiff_cmd")
        .bdgdiff
    }
})

.randsample <- local({
    .randsample <- NULL
    function() {
        if (is.null(.randsample))
            .randsample <<- import("MACS3.Commands.randsample_cmd")
        .randsample
    }
})

.refinepeak <- local({
    .refinepeak <- NULL
    function() {
        if (is.null(.refinepeak))
            .refinepeak <<- import("MACS3.Commands.refinepeak_cmd")
        .refinepeak
    }
})

.callvar <- local({
    .callvar <- NULL
    function() {
        if (is.null(.callvar))
            .callvar <<- import("MACS3.Commands.callvar_cmd")
        .callvar
    }
})


