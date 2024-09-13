#' refinepeak
#'
#' Take raw reads alignment, refine peak summits. Inspired by SPP.
#'
#' @param bedfile Candidate peak file in BED format. REQUIRED.
#' @param ifile ChIP-seq alignment file. If multiple files are given
#'     as '-t A B C', then they will all be read and combined. Note
#'     that pair-end data is not supposed to work with this
#'     command. REQUIRED.
#' @param format Format of tag file, \"AUTO\", \"BED\" or \"ELAND\" or
#'     \"ELANDMULTI\" or \"ELANDEXPORT\" or \"SAM\" or \"BAM\" or
#'     \"BOWTIE\". The default AUTO option will let '%(prog)s' decide
#'     which format the file is. Please check the definition in README
#'     file if you choose
#'     ELAND/ELANDMULTI/ELANDEXPORT/SAM/BAM/BOWTIE. DEFAULT: \"AUTO\""
#' @param cutoff Cutoff. Regions with SPP wtd score lower than cutoff will not be considerred.  DEFAULT: 5
#' 
#' @param windowsize Scan window size on both side of the summit
#'     (default: 100bp)
#' @param buffer_size Buffer size for incrementally increasing
#'     internal array size to store reads alignment information. In
#'     most cases, you don't have to change this parameter. However,
#'     if there are large number of chromosomes/contigs/scaffolds in
#'     your alignment, it's recommended to specify a smaller buffer
#'     size in order to decrease memory usage (but it will take longer
#'     time to read alignment files). Minimum memory requested for
#'     reading an alignment file is about # of CHROMOSOME *
#'     BUFFER_SIZE * 8 Bytes. DEFAULT: 100000
#' @param verbose Set verbose level. 0: only show critical message, 1:
#'     show additional warning message, 2: show process information,
#'     3: show debug messages. If you want to know where are the
#'     duplicate reads, use 3. DEFAULT:2
#' @param outputfile Output bedGraph file name. If not specified, will
#'     write to standard output. REQUIRED.
#' @param oprefix Output file prefix. Mutually exclusive with -o/--ofile.
#' @param outdir Output file name. Mutually exclusive with --o-prefix.
#' @param log Whether to capture logs.
#' @return `macsList` object.
#' @export
#' @examples
#' eh <- ExperimentHub::ExperimentHub()
#' CHIP <- eh[["EH4558"]]
#' CTRL <- eh[["EH4563"]]
#' res <- callpeak(CHIP, CTRL, gsize = 5.2e7, cutoff_analysis = TRUE,
#'                 outdir = tempdir(), name = "callpeak_narrow0")
#' refinepeak(grep("narrowPeak", res$outputs, value = TRUE), CHIP,
#'            outdir = tempdir(), outputfile = "refine")
refinepeak <- function(bedfile, ifile,
                       format = c("AUTO","BAM","SAM","BED","ELAND",
                                  "ELANDMULTI","ELANDEXPORT","BOWTIE"),
                       cutoff = 5, windowsize = 200L,
                       buffer_size = 100000L, verbose = 2L,
                       outdir = "./", outputfile = character(), 
                       oprefix = character(), log =TRUE){
    format <- match.arg(format)
    if(is.character(ifile)){
        ifile <- as.list(normalizePath(ifile))
    }
    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.namespace, outdir){
        opts <- .namespace()$Namespace(bedfile = bedfile,
                                       ifile = ifile,
                                       format = format,
                                       cutoff = cutoff,
                                       windowsize = windowsize,
                                       buffer_size = buffer_size,
                                       verbose = verbose,
                                       ofile = outputfile,
                                       oprefix = oprefix,
                                       outdir = outdir)
        .refinepeak <- reticulate::import("MACS3.Commands.refinepeak_cmd")
        if(log){
            reticulate::py_capture_output(.refinepeak$run(opts))
        }else{
            .refinepeak$run(opts)
        }
    }, .namespace = .namespace, outdir = outdir)
    if(log){
        message(res)
    }

    ofile <- file.path(outdir, outputfile)
    args <- as.list(match.call())
    macsList(arguments = args, outputs = ofile, log = res)
}
