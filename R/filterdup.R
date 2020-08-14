#' filterdup
#'
#' @param gsize Effective genome size. It can be 1.0e+9 or 1000000000,
#'     or shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9),
#'     'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8),
#'     Default:hs.
#' @param format Input file format.
#' @param keepduplicates It controls the behavior towards duplicate
#'     tags at the exact same location -- the same coordination and
#'     the same strand. The 'auto' option makes MACS calculate the
#'     maximum tags at the exact same location based on binomal
#'     distribution using 1e-5 as pvalue cutoff; and the 'all' option
#'     keeps every tags.  If an integer is given, at most this number
#'     of tags will be kept at the same location. Note, if you've used
#'     samtools or picard to flag reads as 'PCR/Optical duplicate' in
#'     bit 1024, MACS2 will still read them although the reads may be
#'     decided by MACS2 as duplicate later. If you plan to rely on
#'     samtools/picard/any other tool to filter duplicates, please
#'     remove those duplicate reads and save a new alignment file then
#'     ask MACS2 to keep all by '--keep-dup all'. The default is to
#'     keep one tag at the same location. Default: 1".
#' @param ifile Input file(s).
#' @param outputfile The output file.
#' @param outdir The output directory.
#' @param verbose Set verbose level of runtime message. 0: only show
#'     critical message, 1: show additional warning message, 2: show
#'     process information, 3: show debug messages.  DEFAULT: 2.
#' @param buffer_size Buffer size for incrementally increasing
#'     internal array size to store reads alignment information. In
#'     most cases, you don't have to change this parameter. However,
#'     if there are large number of chromosomes/contigs/scaffolds in
#'     your alignment, it's recommended to specify a smaller buffer
#'     size in order to decrease memory usage (but it will take longer
#'     time to read alignment files). Minimum memory requested for
#'     reading an alignment file is about # of CHROMOSOME *
#'     BUFFER_SIZE * 8 Bytes. DEFAULT: 100000.
#' @param dryrun When set, filterdup will only output numbers instead
#'     of writing output files, including maximum allowable
#'     duplicates, total number of reads before filtering, total
#'     number of reads after filtering, and redundant rate. Default:
#'     not set
#' @param intern Whether to load output file.
#' @export
#' @examples
#' \dontrun{
#' filterdup(ifile =
#'     list("tests/testthat/CTCF_PE_CTRL_chr22_50k.bedpe.gz"),
#'     outputfile = "test.bed", outdir = "/tmp"))
#' }
filterdup <- function(gsize = 2.7e+09, format = "BEDPE", keepduplicates = "auto",
                      ifile = list(), outputfile = character(),
                      outdir = character(), verbose = 1L,
                      buffer_size = 10000, dryrun = FALSE,
                      intern = FALSE){
    if(is.character(ifile)){
        ifile <- as.list(ifile)
    }
    opts <- .namespace$Namespace(gsize = gsize,
                                 format = format,
                                 keepduplicates = keepduplicates,
                                 verbose = verbose,
                                 outputfile = outputfile,
                                 outdir = outdir,
                                 ifile = ifile,
                                 buffer_size = buffer_size,
                                 dryrun = dryrun)
    .filterdup()$run(opts)
    ofile <- file.path(outdir, outputfile)
    if(intern == TRUE){
        
    }
    return(ofile)
}
