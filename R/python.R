#' @importFrom reticulate import py_run_string
.namespace <- local({
    .namespace <- NULL
    function() {
        if (is.null(.namespace)){
            .namespace <<- reticulate::py_run_string("
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
            .logging <<- reticulate::py_run_string("
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
