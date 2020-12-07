#' bdgdiff
#'
#' Differential peak detection based on paired four bedgraph
#' files. Note: All regions on the same chromosome in the bedGraph
#' file should be continuous so only bedGraph files from MACS3 are
#' accpetable.
#'
#' @param t1bdg MACS pileup bedGraph for condition 1. Incompatible
#'     with callpeak --SPMR output. REQUIRED
#' @param t2bdg MACS pileup bedGraph for condition 2. Incompatible
#'     with callpeak --SPMR output. REQUIRED
#' @param c1bdg MACS control lambda bedGraph for condition
#'     1. Incompatible with callpeak --SPMR output. REQUIRED
#' @param c2bdg MACS control lambda bedGraph for condition
#'     2. Incompatible with callpeak --SPMR output. REQUIRED
#' @param cutoff logLR cutoff. DEFAULT: 3 (likelihood
#'     ratio=1000)", default = 3
#' @param minlen Minimum length of differential region. Try bigger value to remove small regions. DEFAULT: 200",
#'     default = 200
#' @param maxgap Maximum gap to merge nearby differential
#'     regions. Consider a wider gap for broad marks. Maximum gap
#'     should be smaller than minimum length (-g). DEFAULT:
#'     100", default = 100
#' @param depth1 Sequencing depth (# of non-redundant reads in
#'     million) for condition 1. It will be used together with
#'     --d2. See description for --d2 below for how to assign
#'     them. Default: 1
#' @param depth2 Sequencing depth (# of non-redundant reads in
#'     million) for condition 2. It will be used together with
#'     --d1. DEPTH1 and DEPTH2 will be used to calculate scaling
#'     factor for each sample, to down-scale larger sample to the
#'     level of smaller one. For example, while comparing 10 million
#'     condition 1 and 20 million condition 2, use --d1 10 --d2 20,
#'     then pileup value in bedGraph for condition 2 will be divided
#'     by 2. Default: 1
#' 
#' @param oprefix Output file prefix. Actual files will be named as
#'     PREFIX_cond1.bed, PREFIX_cond2.bed and
#'     PREFIX_common.bed. Mutually exclusive with -o/--ofile.
#' @param outputfile Output filenames. Must give three arguments in
#'     order: 1. file for unique regions in condition 1; 2. file for
#'     unique regions in condition 2; 3. file for common regions in
#'     both conditions. Note: mutually exclusive with --o-prefix.
#' @param outdir The output directory.
#' @param log Whether to capture logs.
#' @examples
#' \dontrun{
#' }
bdgdiff <- function(t1bdg, t2bdg, c1bdg, c2bdg,
                    cutoff = 3, minlen = 200L, maxgap = 100L,
                    depth1 = 1, depth2 = 1,
                    outdir = ".",
                    oprefix = character(),
                    outputfile = list(),
                    log = TRUE){
    method <- lapply(method, function(x)match.arg(x, method))
    opts <- .namespace()$Namespace(t1bdg = t1bdg,
                                   t2bdg = t2bdg,
                                   c1bdg = c1bdg,
                                   c2bdg = c2bdg,
                                   cutoff = cutoff,
                                   minlen = minlen,
                                   maxgap = maxgap,
                                   depth1 = depth1,
                                   depth2 = depth2,
                                   oprefix = oprefix,
                                   ofile = outputfile,
                                   outdir = outdir)
    if(log){
        .logging()$run()
        res <- py_capture_output(.bdgdiff()$run(opts))
        message(res)
    }else{
        res <- .bdgdiff()$run(opts)
    }

    if(length(oprefix) > 0){
        ofile <- list.files(outdir, paste0("^", oprefix, "*"), full.names = TRUE)
    }else{
        ofile <- file.path(outdir, outputfile)
    }

    args <- as.list(match.call())
    macsList(fun = args[[1]], arguments = args[-1], outputs = ofile, log = res)
}
