#' cmbreps
#'
#' Combine BEDGraphs of scores from replicates. Note: All regions on
#' the same chromosome in the bedGraph file should be continuous so
#' only bedGraph files from MACS3 are accpetable.
#'
#' @param ifiles MACS score in bedGraph for each replicate. Require at
#'     least 2 files such as '-i A B C D'. REQUIRED
#' @param method to use while combining scores from replicates. 1)
#'     fisher: Fisher's combined probability test. It requires scores
#'     in ppois form (-log10 pvalues) from bdgcmp. Other types of
#'     scores for this method may cause cmbreps unexpected errors. 2)
#'     max: take the maximum value from replicates for each genomic
#'     position. 3) mean: take the average value. Note, except for
#'     Fisher's method, max or mean will take scores AS IS which means
#'     they won't convert scores from log scale to linear scale or
#'     vice versa.", default="fisher"
#' @param outputfile Output BEDGraph filename for combined scores.
#' @param outdir The output directory.
#' @param log Whether to capture logs.
#' @param verbose Set verbose level of runtime message. 0: only show
#'     critical message, 1: show additional warning message, 2: show
#'     process information, 3: show debug messages. DEFAULT:2
#' @return `macsList` object.
#' @export
#' @examples
#' eh <- ExperimentHub::ExperimentHub()
#' CHIP <- eh[["EH4558"]]
#' CTRL <- eh[["EH4563"]]
#' c1 <- callpeak(CHIP, CTRL, gsize = 5.2e7, cutoff_analysis = TRUE,
#'                outdir = tempdir(), name = "callpeak_narrow0",
#'                store_bdg = TRUE)
#' cmbreps(ifiles = list(c1$outputs[1], c1$outputs[7]),
#'         method = "max", outdir = tempdir(), outputfile = "cmbreps")
cmbreps <- function(ifiles = list(),
                    method = c("fisher", "max", "mean"),
                    outputfile = character(),
                    outdir = ".", log = TRUE, verbose = 2L){
    method <- match.arg(method)
    names(ifiles) <- NULL
    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.namespace, outdir){
        opts <- .namespace()$Namespace(ifile = ifiles,
                                       method = method,
                                       ofile = outputfile,
                                       outdir = outdir,
                                       verbose = verbose)
        .cmbreps <- reticulate::import("MACS3.Commands.cmbreps_cmd")
        if(log){
            reticulate::py_capture_output(.cmbreps$run(opts))
        }else{
            .cmbreps$run(opts)
        }
    }, .namespace = .namespace, outdir = outdir)
    if(log){
        message(res)
    }

    ofile <- file.path(outdir, outputfile)
    args <- as.list(match.call())
    macsList(arguments = args, outputs = ofile, log = res)
}
