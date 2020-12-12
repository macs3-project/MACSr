#' bdgbroadcall
#'
#' Call broad peaks from bedGraph output. Note: All regions on the
#' same chromosome in the bedGraph file should be continuous so only
#' bedGraph files from MACS3 are accpetable.
#'
#' @param ifile MACS score in bedGraph. REQUIRED.
#' @param cutoffpeak Cutoff for peaks depending on which method you
#'     used for score track. If the file contains qvalue scores from
#'     MACS3, score 2 means qvalue 0.01. DEFAULT: 2
#' @param cutofflink Cutoff for linking regions/low abundance regions
#'     depending on which method you used for score track. If the file
#'     contains qvalue scores from MACS3, score 1 means qvalue 0.1,
#'     and score 0.3 means qvalue 0.5. DEFAULT: 1", default = 1
#' @param minlen minimum length of peak, better to set it as d value. DEFAULT: 200",
#'     default = 200
#' @param lvl1maxgap maximum gap between significant peaks, better to
#'     set it as tag size. DEFAULT: 30
#' @param lvl2maxgap maximum linking between significant peaks, better
#'     to set it as 4 times of d value. DEFAULT: 800
#' @param trackline Tells MACS not to include trackline with bedGraph
#'     files. The trackline is required by UCSC.
#' @param outputfile The output file.
#' @param outdir The output directory.
#' @param log Whether to capture logs.
#' @examples
#' \dontrun{
#' bdgbroadcall("run_bdgcmp_FE.bdg", cutoffpeak = 2, cutofflink = 1.5)
#' }
bdgbroadcall <- function(ifile, cutoffpeak = 2, cutofflink = 1,
                         minlen = 200L, lvl1maxgap = 30L,
                         lvl2maxgap = 800L, trackline = TRUE,
                         outdir = ".", outputfile = character(),
                         log = TRUE){
    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.logging, .namespace, outdir){
        opts <- .namespace()$Namespace(ifile = file.path(ifile),
                                       cutoffpeak = cutoffpeak,
                                       cutofflink = cutofflink,
                                       minlen = minlen,
                                       lvl1maxgap = lvl1maxgap,
                                       lvl2maxgap = lvl2maxgap,
                                       trackline = trackline,
                                       outdir = outdir,
                                       ofile = outputfile)
        .bdgbroadcall <- reticulate::import("MACS3.Commands.bdgbroadcall_cmd")
        if(log){
            .logging()$run()
            reticulate::py_capture_output(.bdgbroadcall$run(opts))
        }else{
            .bdgbroadcall$run(opts)
        }
    }, .logging = .logging, .namespace = .namespace, outdir = outdir)
    if(log){
        message(res)
    }
    ofile <- file.path(outdir, outputfile)
    args <- as.list(match.call())
    macsList(arguments = args, outputs = ofile, log = res)
}
