#' predictd
#' 
#' @param ifile Input file(s).
#' @param gsize Effective genome size. It can be 1.0e+9 or 1000000000,
#'     or shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9),
#'     'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8),
#'     Default:hs.
#' @param format Input file format.
#' @param plot PDF path of peak model and correlation plots.
#' @param tsize Tag size. This will override the auto detected tag
#'     size.
#' @param bw Band width for picking regions to compute fragment
#'     size. This value is only used while building the shifting
#'     model. DEFAULT: 300
#' @param d_min Minimum fragment size in basepair. Any predicted
#'     fragment size less than this will be excluded. DEFAULT: 20
#' @param mfold Select the regions within MFOLD range of
#'     high-confidence enrichment ratio against background to build
#'     model. Fold-enrichment in regions must be lower than upper
#'     limit, and higher than the lower limit. Use as
#'     "-m 10 30". DEFAULT:5 50
#' @param buffer_size Buffer size for incrementally increasing
#'     internal array size to store reads alignment
#'     information. DEFAULT: 100000.
#' @param verbose Set verbose level of runtime message. 0: only show
#'     critical message, 1: show additional warning message, 2: show
#'     process information, 3: show debug messages. DEFAULT:2
#' @param log Whether to capture log.
#' @examples
#' \dontrun{
#' predictd("extdata/CTCF_SE_ChIP_chr22_50k.bed.gz", d_min = 10)
#' }

predictd <- function(ifile, gsize = "hs", format = "AUTO",
                     plot = file.path(tempdir(), "predictd_mode.pdf"),
                     tsize = NULL, bw = 300, d_min = 20, mfold = c(5, 50),
                     buffer_size = 100000, verbose = 2L, log = TRUE){
    if(is.character(ifile)){
        ifile <- as.list(file.path(ifile))
    }
    rfile <- tempfile()
    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.logging, .namespace, rfile){
        opts <- .namespace()$Namespace(ifile = ifile,
                                       gsize = gsize,
                                       format = format,
                                       tsize = tsize,
                                       bw = bw,
                                       d_min = d_min,
                                       mfold = mfold,
                                       outdir = dirname(rfile),
                                       rfile = basename(rfile),
                                       verbose = verbose,
                                       buffer_size = buffer_size)
        .predictd <- reticulate::import("MACS3.Commands.predictd_cmd")
        if(log){
            .logging()$run()
            reticulate::py_capture_output(.predictd$run(opts))
        }else{
            .predictd$run(opts)
        }
    }, .logging = .logging, .namespace = .namespace, rfile = rfile)
    if(log){
        message(res)
    }

    env <- new.env()
    rs <- readLines(rfile)
    rs[grep("^pdf", rs)] <- paste0("pdf('", rfile, "_model.pdf',height=6,width=6)")
    source(textConnection(rs), local = env)

    if(!is.null(plot)){
        fc <- file.copy(paste0(rfile, "_model.pdf"), plot)
        message("model plot:", plot)
    }
    altd <- get("altd", envir = env)
    return(altd)
}
