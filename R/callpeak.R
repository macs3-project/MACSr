#' callpeak
#'
#' Main MACS3 Function to call peaks from alignment results.
#' @param tfile ChIP-seq treatment files.
#' @param cfile Control files.
#' @param gsize Effective genome size. It can be 1.0e+9 or 1000000000,
#'     or shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9),
#'     'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8),
#'     Default:hs.
#' @param tsize Tag size/read length. This will override the auto
#'     detected tag size. DEFAULT: Not set
#' @param format Format of tag file, "AUTO", "BED" or "ELAND" or
#'     "ELANDMULTI" or "ELANDEXPORT" or "SAM" or "BAM" or "BOWTIE" or
#'     "BAMPE" or "BEDPE".
#' @param keepduplicates It controls the behavior towards duplicate
#'     tags at the exact same location -- the same coordination and
#'     the same strand.
#' @param outdir If specified all output files will be written to that
#'     directory.
#' @param name Experiment name, which will be used to generate output
#'     file names.
#' @param store_bdg Whether or not to save extended fragment pileup,
#'     and local lambda tracks (two files) at every bp into a bedGraph
#'     file.
#' @param do_SPMR If True, MACS will SAVE signal per million reads for
#'     fragment pileup profiles.
#' @param trackline Tells MACS to include trackline with bedGraph
#'     files.
#' @param nomodel Whether or not to build the shifting model.
#' @param shift The arbitrary shift in bp. Use discretion while
#'     setting it other than default value.
#' @param extsize The arbitrary extension size in bp.
#' @param bw Band width for picking regions to compute fragment size.
#' @param d_min Minimum fragment size in basepair. Any predicted
#'     fragment size less than this will be excluded.
#' @param mfold Select the regions within MFOLD range of
#'     high-confidence enrichment ratio against background to build
#'     model.
#' @param onauto Whether turn on the auto pair model process.
#' @param qvalue Minimum FDR (q-value) cutoff for peak detection.
#' @param pvalue Pvalue cutoff for peak detection. DEFAULT: not set.
#' @param tempdir Optional directory to store temp files.
#' @param nolambda If True, MACS will use fixed background lambda as
#'     local lambda for every peak region.
#' @param scaleto When set to 'small', scale the larger sample up to
#'     the smaller sample.
#' @param downsample When set, random sampling method will scale down
#'     the bigger sample. By default, MACS uses linear scaling.
#' @param slocal The small nearby region in basepairs to calculate
#'     dynamic lambda.
#' @param llocal The large nearby region in basepairs to calculate
#'     dynamic lambda.
#' @param broad If set, MACS will try to call broad peaks using the
#'     --broad-cutoff setting.
#' @param broadcutoff Cutoff for broad region. This option is not
#'     available unless --broad is set.
#' @param maxgap Maximum gap between significant sites to cluster them
#'     together. The DEFAULT value is the detected read length/tag
#'     size.
#' @param minlen Minimum length of a peak. The DEFAULT value is the predicted
#'     fragment size d.
#' @param cutoff_analysis While set, MACS2 will analyze number or
#'     total length of peaks that can be called by different p-value
#'     cutoff then output a summary table to help user decide a better
#'     cutoff.
#' @param fecutoff When set, the value will be used to filter out
#'     peaks with low fold-enrichment.
#' @param call_summits If set, MACS will use a more sophisticated
#'     signal processing approach to find subpeak summits in each
#'     enriched peak region.
#' @param buffer_size Buffer size for incrementally increasing
#'     internal array size to store reads alignment
#'     information. DEFAULT: 100000.
#' @param verbose Set verbose level of runtime message. 0: only show
#'     critical message, 1: show additional warning message, 2: show
#'     process information, 3: show debug messages. DEFAULT:2
#' @param log Whether to capture logs.
#' @param ... More options for macs2.
#' @importFrom reticulate py_capture_output
#' @export
#' @examples
#' \dontrun{
#' callpeak("inst/extdata/CTCF_SE_ChIP_chr22_50k.bed.gz",
#'         "inst/extdata/CTCF_SE_CTRL_chr22_50k.bed.gz",
#'         gsize = 5.2e7, cutoff_analysis = T)
#' }

callpeak <- function(tfile, cfile = NULL, gsize = "hs", tsize = NULL, format = "AUTO", keepduplicates = "1",
                     outdir = ".", name = "NA", store_bdg = FALSE, do_SPMR = FALSE, trackline = FALSE,
                     nomodel = FALSE, shift = 0, extsize = 200, bw = 300, d_min = 20,
                     mfold = c(5, 50), onauto = FALSE, qvalue = 0.05, pvalue = NULL,
                     tempdir = "/tmp", nolambda = FALSE, scaleto = "small", downsample = FALSE,
                     slocal = 1000, llocal = 10000, broad = FALSE, broadcutoff = 0.1,
                     maxgap = NULL, minlen = NULL,
                     cutoff_analysis = FALSE, fecutoff = 0.1, call_summits = FALSE,
                     buffer_size = 100000, verbose = 2L, log = TRUE, ...){
    if(is.character(tfile)){
        tfile <- as.list(tfile)
    }
    if(is.character(cfile)){
        cfile <- as.list(cfile)
    }
    opts <- .namespace()$Namespace(tfile = tfile,
                                   cfile = cfile,
                                   gsize = gsize,
                                   tsize = tsize,
                                   format = format,
                                   keepduplicates = keepduplicates,
                                   outdir = outdir,
                                   name = name,
                                   store_bdg = store_bdg,
                                   do_SPMR = do_SPMR,
                                   trackline = trackline,
                                   nomodel = nomodel,
                                   shift = shift,
                                   extsize = extsize,
                                   bw = bw,
                                   d_min = d_min,
                                   mfold = mfold,
                                   onauto = onauto,
                                   qvalue = qvalue,
                                   pvalue = pvalue,
                                   tempdir = tempdir,
                                   nolambda = nolambda,
                                   scaleto = scaleto,
                                   downsample = downsample,
                                   smalllocal = slocal,
                                   largelocal = llocal,
                                   broad = broad,
                                   broadcutoff = broadcutoff,
                                   maxgap = maxgap,
                                   minlen = minlen,
                                   cutoff_analysis = cutoff_analysis,
                                   fecutoff = fecutoff,
                                   call_summits = call_summits,
                                   buffer_size = buffer_size,
                                   verbose = verbose,
                                   ratio = NA, ...)
    if(log){
        .logging()$run()
        res <- py_capture_output(.callpeak()$run(opts))
        message(res)
    }else{
        res <- .callpeak()$run(opts)
    }
    outputs <- list.files(path = outdir, pattern = paste0(name, "_.*"), full.names = TRUE)
    args <- as.list(match.call())
    macsList(fun = args[[1]], arguments = args[-1], outputs = outputs, log = res)
}
