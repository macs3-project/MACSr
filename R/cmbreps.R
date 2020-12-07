#' cmbreps
#'
#' Combine BEDGraphs of scores from replicates. Note: All regions on
#' the same chromosome in the bedGraph file should be continuous so
#' only bedGraph files from MACS3 are accpetable.
#'
#' @param ifiles MACS score in bedGraph for each replicate. Require at
#'     least 2 files such as '-i A B C D'. REQUIRED
#' @param weights Weight for each replicate. Default is 1.0 for
#'     each. When given, require same number of parameters as IFILE.
#' @param method to use while combining scores from replicates. 1)
#'     fisher: Fisher's combined probability test. It requires scores
#'     in ppois form (-log10 pvalues) from bdgcmp. Other types of
#'     scores for this method may cause cmbreps unexpected errors. 2)
#'     max: take the maximum value from replicates for each genomic
#'     position. 3) mean: take the average value. Note, except for
#'     Fisher's method, max or mean will take scores AS IS which means
#'     they won't convert scores from log scale to linear scale or
#'     vice versa.", default="fisher"
#' @param outputfile Output filename. Mutually exclusive with
#'     --o-prefix. The number and the order of arguments for --ofile
#'     must be the same as for -m.
#' @param outdir The output directory.
#' @param log Whether to capture logs.
#' @examples
#' \dontrun{
#' bdgopt("run_callpeak_narrow0_treat_pileup.bdg",
#' method = "min", extraparam = 10, outdir = "/tmp", outputfile = "bdgopt_min.bdg")
#' }

cmbreps <- function(ifiles = list(), weights = 1.0,
                    method = c("fisher", "max", "mean"),
                    outputfile = character(),
                    outdir = ".", log = TRUE){
    method <- match.arg(method)
    opts <- .namespace()$Namespace(ifile = ifiles,
                                   weights = weights,
                                   method = method,
                                   ofile = outputfile,
                                   outdir = outdir)
    if(log){
        .logging()$run()
        res <- py_capture_output(.cmbreps()$run(opts))
        message(res)
    }else{
        res <- .cmbreps()$run(opts)
    }

    ofile <- file.path(outdir, outputfile)
    args <- as.list(match.call())
    macsList(fun = args[[1]], arguments = args[-1], outputs = ofile, log = res)
}
