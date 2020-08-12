#' filterdup
#'
#' @param gsize
#' @param format
#' @param keepduplicates
#' @param outputfile
#' @param outdir
#' @param verbose
#' @param buffer_size
#' @param dryrun
#' @export
#' @examples
#' ## filterdup(ifile = list("tests/testthat/CTCF_PE_CTRL_chr22_50k.bedpe.gz", outputfile = "test.bed", outdir = "/tmp"))
filterdup <- function(gsize = 2.7e+09, format = "BEDPE", keepduplicates = "auto",
                      ifile = list(), outputfile = character(),
                      outdir = character(), verbose = 1L,
                      buffer_size = 10000, dryrun = FALSE,
                      intern = FALSE){
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
    if(intern = TRUE){
        ofile <- system.file(outdir, outputfile)
    }
    
}
