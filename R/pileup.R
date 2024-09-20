#' pileup
#'
#' Pileup aligned reads with a given extension size (fragment size or
#' d in MACS language). Note there will be no step for duplicate reads
#' filtering or sequencing depth scaling, so you may need to do
#' certain pre/post-processing.
#'
#' @param ifile Alignment file. If multiple files are given as '-t A B C', then they will all be read and combined. REQUIRED.
#' @param format Format of tag file, \"AUTO\", \"BED\", \"ELAND\",
#'     \"ELANDMULTI\", \"ELANDEXPORT\", \"SAM\", \"BAM\", \"BOWTIE\",
#'     \"BAMPE\", or \"BEDPE\". The default AUTO option will let
#'     '%(prog)s' decide which format the file is. DEFAULT: \"AUTO\",
#'     MACS3 will pick a format from \"AUTO\", \"BED\", \"ELAND\",
#'     \"ELANDMULTI\", \"ELANDEXPORT\", \"SAM\", \"BAM\" and
#'     \"BOWTIE\". If the format is BAMPE or BEDPE, please specify it
#'     explicitly. Please note that when the format is BAMPE or BEDPE,
#'     the -B and --extsize options would be ignored.
#' @param bothdirection By default, any read will be extended towards downstream direction by extension size. So it's (0,size-1) (1-based index system) for plus strand read and (-size+1,0) for minus strand read where position 0 is 5' end of the aligned read. Default behavior can simulate MACS3 way of piling up ChIP sample reads where extension size is set as fragment size/d. If this option is set as on, aligned reads will be extended in both upstream and downstream directions by extension size. It means (-size,size) where 0 is the 5' end of a aligned read. It can partially simulate MACS3 way of piling up control reads. However MACS3 local bias is calculated by maximizing the expected pileup over a ChIP fragment size/d estimated from 10kb, 1kb, d and whole genome background. This option will be ignored when the format is set as BAMPE or BEDPE. DEFAULT: False
#' @param extsize The extension size in bps. Each alignment read will
#'     become a EXTSIZE of fragment, then be piled up. Check
#'     description for -B for detail. It's twice the `shiftsize` in
#'     old MACSv1 language. This option will be ignored when the
#'     format is set as BAMPE or BEDPE. DEFAULT: 200
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
#' @param outputfile Output bedGraph file name. If not specified, will write to standard output. REQUIRED.
#' @param outdir The output directory.
#' @param log Whether to capture logs.
#' @return `macsList` object.
#' @export
#' @examples
#' eh <- ExperimentHub::ExperimentHub()
#' CHIP <- eh[["EH4558"]]
#' p <- pileup(CHIP, outdir = tempdir(), outputfile = "pileup_bed.bdg", format = "BED")
pileup <- function(ifile, outputfile = character(), outdir = ".",
                   format = c("AUTO","BAM","SAM","BED","ELAND","ELANDMULTI",
                       "ELANDEXPORT","BOWTIE","BAMPE","BEDPE"),
                   bothdirection = FALSE, extsize = 200L,
                   buffer_size = 100000L, verbose = 2L,
                   log = TRUE){
    format <- match.arg(format)
    if(is.character(ifile)){
        ifile <- as.list(normalizePath(ifile))
    }

    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.namespace, outdir){
        opts <- .namespace()$Namespace(ifile = ifile,
                                       format = format,
                                       bothdirection = bothdirection,
                                       extsize = extsize,
                                       buffer_size = buffer_size,
                                       verbose = verbose,
                                       outputfile = outputfile,
                                       outdir = outdir)
        .pileup <- reticulate::import("MACS3.Commands.pileup_cmd")
        if(log){
            reticulate::py_capture_output(.pileup$run(opts))
        }else{
            .pileup$run(opts)
        }
    }, .namespace = .namespace, outdir = outdir)
    if(log){
        message(res)
    }

    ofile <- file.path(outdir, outputfile)
    args <- as.list(match.call())
    macsList(arguments = args, outputs = ofile, log = res)
}
