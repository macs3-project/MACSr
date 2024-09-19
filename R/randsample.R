#' randsample
#'
#' Randomly choose a number/percentage of total reads, then save in BED/BEDPE format file.
#'
#' @param ifile Alignment file. If multiple files are given as '-t A B C', then they will all be read and combined. REQUIRED.
#' @param percentage Percentage of tags you want to keep. Input 80.0 for 80%%. This option can't be used at the same time with -n/--num. If the setting is 100, it will keep all the reads and convert any format that MACS3 supports into BED or BEDPE (if input is BAMPE) format. REQUIRED
#' @param number Number of tags you want to keep. Input 8000000 or 8e+6 for 8 million. This option can't be used at the same time with -p/--percent. Note that the number of tags in output is approximate as the number specified here. REQUIRED
#' @param seed Set the random seed while down sampling data. Must be a non-negative integer in order to be effective. If you want more reproducible results, please specify a random seed and record it.DEFAULT: not set
#' @param tsize Tag size. This will override the auto detected tag
#'     size. DEFAULT: Not set
#' @param format Format of tag file, \"AUTO\", \"BED\" or \"ELAND\" or
#'     \"ELANDMULTI\" or \"ELANDEXPORT\" or \"SAM\" or \"BAM\" or
#'     \"BOWTIE\" or \"BAMPE\" or \"BEDPE\". The default AUTO option
#'     will %(prog)s decide which format the file is. Please check the
#'     definition in README file if you choose
#'     ELAND/ELANDMULTI/ELANDEXPORT/SAM/BAM/BOWTIE or
#'     BAMPE/BEDPE. DEFAULT: \"AUTO\""
#' 
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
#' @param outputfile Output BED file name. If not specified, will write to standard output. Note, if the input format is BAMPE or BEDPE, the output will be in BEDPE format. DEFAULT: stdout
#' @param outdir The output directory.
#' @param log Whether to capture logs.
#' @return `macsList` object.
#' @export
#' @examples
#' eh <- ExperimentHub::ExperimentHub()
#' CHIP <- eh[["EH4558"]]
#' randsample(CHIP, number = 1000, outdir = tempdir(), outputfile = "randsample.bed")
randsample <- function(ifile, outdir = ".", outputfile = character(),
                       percentage = numeric(), number = numeric(),
                       seed = -1L, tsize = NULL,
                       format = c("AUTO","BAM","SAM","BED","ELAND","ELANDMULTI",
                                  "ELANDEXPORT","BOWTIE","BAMPE","BEDPE"),
                       buffer_size = 100000L, verbose = 2L, log =TRUE){
    format <- match.arg(format)
    if(is.character(ifile)){
        ifile <- as.list(normalizePath(ifile))
    }

    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.namespace, outdir){
        opts <- .namespace()$Namespace(ifile = ifile,
                                       percentage = percentage,
                                       number = number,
                                       seed = seed,
                                       tsize = tsize,
                                       format = format,
                                       buffer_size = buffer_size,
                                       verbose = verbose,
                                       outputfile = outputfile,
                                       outdir = outdir)
        .randsample <- reticulate::import("MACS3.Commands.randsample_cmd")
        if(log){
            reticulate::py_capture_output(.randsample$run(opts))
        }else{
            .randsample$run(opts)
        }
    }, .namespace = .namespace, outdir = outdir)
    if(log){
        message(res)
    }

    ofile <- file.path(outdir, outputfile)
    args <- as.list(match.call())
    macsList(arguments = args, outputs = ofile, log = res)
}
