#' predictd
#' 
#' Predict d or fragment size from alignment results. In case of PE
#' data, report the average insertion/fragment size from all
#' pairs. *Will NOT filter duplicates*
#' 
#' @param ifile ChIP-seq alignment file. If multiple files are given as '-t A B C', then they will all be read and combined. REQUIRED.
#' @param gsize Effective genome size. It can be 1.0e+9 or 1000000000, or shortcuts:'hs' for human (2,913,022,398), 'mm' for mouse (2,652,783,500), 'ce' for C. elegans (100,286,401) and 'dm' for fruitfly (142,573,017), Default:hs. The effective genome size numbers for the above four species are collected from Deeptools https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html Please refer to deeptools to define the best genome size you plan to use.
#' @param format Format of tag file, \"AUTO\", \"BED\" or \"ELAND\" or \"ELANDMULTI\" or \"ELANDEXPORT\" or \"SAM\" or \"BAM\" or \"BOWTIE\" or \"BAMPE\" or \"BEDPE\". The default AUTO option will let MACS decide which format the file is. However, if you want to decide the average insertion size/fragment size from PE data such as BEDPE or BAMPE, please specify the format as BAMPE or BEDPE since MACS3 won't automatically recognize three two formats with -f AUTO. Please be aware that in PE mode, -g, -s, --bw, --d-min, -m, and --rfile have NO effect. DEFAULT: \"AUTO\"
#' @param plot PDF path of peak model and correlation plots.
#' @param tsize Tag size. This will override the auto detected tag
#'     size. DEFAULT: Not set
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
#' @return predicted fragment sizes.
#' @export
#' @examples
#' eh <- ExperimentHub::ExperimentHub()
#' CHIP <- eh[["EH4558"]]
#' predictd(CHIP, d_min = 10, gsize=5.2e+7, plot = NULL)
predictd <- function(ifile, gsize = "hs", format = "AUTO",
                     plot = normalizePath(tempdir(), "predictd_mode.pdf"),
                     tsize = NULL, bw = 300, d_min = 20, mfold = c(5, 50),
                     buffer_size = 100000, verbose = 2L, log = TRUE){
    if(is.character(ifile)){
        ifile <- as.list(normalizePath(ifile))
    }
    rfile <- tempfile()
    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.namespace, rfile){
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
            reticulate::py_capture_output(.predictd$run(opts))
        }else{
            .predictd$run(opts)
        }
    }, .namespace = .namespace, rfile = rfile)
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
