#' bdgopt
#'
#' Operations on score column of bedGraph file. Note: All regions on
#' the same chromosome in the bedGraph file should be continuous so
#' only bedGraph files from MACS3 are accpetable.
#'
#' @param ifile MACS score in bedGraph. Note: this must be a bedGraph
#'     file covering the ENTIRE genome. REQUIRED
#' @param method Method to modify the score column of bedGraph file. Available choices are: multiply, add, max, min, or p2q. 1) multiply, the EXTRAPARAM is required and will be multiplied to the score column. If you intend to divide the score column by X, use value of 1/X as EXTRAPARAM. 2) add, the EXTRAPARAM is required and will be added to the score column. If you intend to subtract the score column by X, use value of -X as EXTRAPARAM. 3) max, the EXTRAPARAM is required and will take the maximum value between score and the EXTRAPARAM. 4) min, the EXTRAPARAM is required and will take the minimum value between score and the EXTRAPARAM. 5) p2q, this will convert p-value scores to q-value scores using Benjamini-Hochberg process. The EXTRAPARAM is not required. This method assumes the scores are -log10 p-value from MACS3. Any other types of score will cause unexpected errors. Default="p2q"
#' 
#' @param extraparam The extra parameter for METHOD. Check the detail
#'     of -m option.
#' @param outputfile Output BEDGraph filename. Required.
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
#' cfile <- grep("treat_pileup.bdg", c1$outputs, value = TRUE)
#' bdgopt(cfile, method = "min", extraparam = 10,
#'        outdir = tempdir(), outputfile = "bdgopt_min.bdg")
bdgopt <- function(ifile,
                   method = c("multiply", "add", "p2q", "max", "min"),
                   extraparam = numeric(),
                   outputfile = character(),
                   outdir = ".", log = TRUE, verbose = 2L){
    method <- match.arg(method)
    ifile <- normalizePath(ifile)
    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.namespace, outdir){
        opts <- .namespace()$Namespace(ifile = ifile,
                                       method = method,
                                       extraparam = list(extraparam),
                                       ofile = outputfile,
                                       outdir = outdir,
                                       verbose = verbose)
        .bdgopt <- reticulate::import("MACS3.Commands.bdgopt_cmd")
        if(log){
            reticulate::py_capture_output(.bdgopt$run(opts))
        }else{
            .bdgopt$run(opts)
        }
    }, .namespace = .namespace, outdir = outdir)
    if(log){
        message(res)
    }

    ofile <- file.path(outdir, outputfile)
    args <- as.list(match.call())
    macsList(arguments = args, outputs = ofile, log = res)
}

