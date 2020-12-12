#' bdgpeakcall
#'
#' Call peaks from bedGraph output. Note: All regions on the same
#' chromosome in the bedGraph file should be continuous so only
#' bedGraph files from MACS3 are accpetable.
#'
#' @param ifile MACS score in bedGraph. REQUIRED.
#' @param cutoff Cutoff depends on which method you used for score
#'     track. If the file contains pvalue scores from MACS3, score 5
#'     means pvalue 1e-5. DEFAULT: 5", default = 5.
#' @param minlen minimum length of peak, better to set it as d
#'     value. DEFAULT: 200", default = 200.
#' @param maxgap maximum gap between significant points in a peak,
#'     better to set it as tag size. DEFAULT: 30", default = 30.
#' @param call_summits If set, MACS will use a more sophisticated
#'     approach to find all summits in each enriched peak region
#'     DEFAULT: False",default=False.
#' @param cutoff_analysis While set, bdgpeakcall will analyze number
#'     or total length of peaks that can be called by different cutoff
#'     then output a summary table to help user decide a better
#'     cutoff. Note, minlen and maxgap may affect the
#'     results. DEFAULT: False", default = False.
#' @param trackline Tells MACS not to include trackline with bedGraph
#'     files. The trackline is required by UCSC.
#' @param outputfile The output file.
#' @param outdir The output directory.
#' @param log Whether to capture logs.
#' @examples
#' \dontrun{
#' bdgpeakcall("run_bdgcmp_FE.bdg", cutoff = 2)
#' }
bdgpeakcall <- function(ifile, cutoff = 5, minlen = 200L, maxgap = 30L,
                        call_summits = FALSE, cutoff_analysis = FALSE,
                        trackline = TRUE, outdir = ".",
                        outputfile = character(), log = TRUE){
    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.logging, .namespace, outdir){
        opts <- .namespace()$Namespace(ifile = file.path(ifile),
                                       cutoff = cutoff,
                                       minlen = minlen,
                                       maxgap = maxgap,
                                       call_summits = call_summits,
                                       cutoff_analysis = cutoff_analysis,
                                       trackline = trackline,
                                       outdir = outdir,
                                       ofile = outputfile)
        .bdgpeakcall <- reticulate::import("MACS3.Commands.bdgpeakcall_cmd")
        if(log){
            .logging()$run()
            reticulate::py_capture_output(.bdgpeakcall$run(opts))
        }else{
            .bdgpeakcall$run(opts)
        }
    }, .logging = .logging, .namespace = .namespace, outdir = outdir)
    if(log){
        message(res)
    }
    ofile <- file.path(outdir, outputfile)
    args <- as.list(match.call())
    macsList(arguments = args, outputs = ofile, log = res)
}
