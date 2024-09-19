#' bdgbroadcall
#'
#' Call broad peaks from bedGraph output. Note: All regions on the
#' same chromosome in the bedGraph file should be continuous so only
#' bedGraph files from MACS3 are accpetable.
#'
#' @param ifile MACS score in bedGraph. REQUIRED.
#' @param cutoffpeak Cutoff for peaks depending on which method you used for score track. If the file contains qvalue scores from MACS3, score 2 means qvalue 0.01. Regions with signals lower than cutoff will not be considerred as enriched regions. DEFAULT: 2
#' @param cutofflink Cutoff for linking regions/low abundance regions depending on which method you used for score track. If the file contains qvalue scores from MACS3, score 1 means qvalue 0.1, and score 0.3 means qvalue 0.5. DEFAULT: 1
#' @param minlen minimum length of peak, better to set it as d value. DEFAULT: 200",
#'     default = 200
#' @param lvl1maxgap maximum gap between significant peaks, better to
#'     set it as tag size. DEFAULT: 30
#' @param lvl2maxgap maximum linking between significant peaks, better
#'     to set it as 4 times of d value. DEFAULT: 800
#' @param trackline Tells MACS not to include trackline with bedGraph
#'     files. The trackline is required by UCSC.
#' @param outputfile Output file name. Mutually exclusive with --o-prefix
#' @param oprefix The PREFIX of output bedGraph file to write
#'     scores. If it is given as A, and method is 'ppois', output file
#'     will be A_ppois.bdg. Mutually exclusive with -o/--ofile.
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
#' p1 <- pileup(CHIP, outdir = tempdir(),
#'              outputfile = "pileup_ChIP_bed.bdg", format = "BED")
#' p2 <- pileup(CTRL, outdir = tempdir(),
#'              outputfile = "pileup_CTRL_bed.bdg", format = "BED")
#' c1 <- bdgcmp(p1$outputs, p2$outputs, outdir = tempdir(),
#'              oprefix = "bdgcmp", pseudocount = 1, method = "FE")
#' bdgbroadcall(c1$outputs, cutoffpeak = 2, cutofflink = 1.5,
#'              outdir = tempdir(), outputfile = "bdgbroadcall")
bdgbroadcall <- function(ifile, cutoffpeak = 2, cutofflink = 1,
                         minlen = 200L, lvl1maxgap = 30L,
                         lvl2maxgap = 800L, trackline = TRUE,
                         outdir = ".", outputfile = character(),
                         oprefix = character(),
                         log = TRUE, verbose = 2L){
    ifile <- normalizePath(ifile)
    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.namespace, outdir){
        opts <- .namespace()$Namespace(ifile = ifile,
                                       cutoffpeak = cutoffpeak,
                                       cutofflink = cutofflink,
                                       minlen = minlen,
                                       lvl1maxgap = lvl1maxgap,
                                       lvl2maxgap = lvl2maxgap,
                                       trackline = trackline,
                                       outdir = outdir,
                                       ofile = outputfile,
                                       oprefix = oprefix,
                                       verbose = verbose)
        .bdgbroadcall <- reticulate::import("MACS3.Commands.bdgbroadcall_cmd")
        if(log){
            reticulate::py_capture_output(.bdgbroadcall$run(opts))
        }else{
            .bdgbroadcall$run(opts)
        }
    }, .namespace = .namespace, outdir = outdir)
    if(log){
        message(res)
    }
    ofile <- file.path(outdir, outputfile)
    args <- as.list(match.call())
    macsList(arguments = args, outputs = ofile, log = res)
}
