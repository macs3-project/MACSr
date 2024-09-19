#' bdgcmp
#'
#' Deduct noise by comparing two signal tracks in bedGraph. Note: All
#' regions on the same chromosome in the bedGraph file should be
#' continuous so only bedGraph files from MACS3 are accpetable.
#'
#' @param tfile Treatment bedGraph file, e.g. *_treat_pileup.bdg from
#'     MACSv2. REQUIRED
#' @param cfile Control bedGraph file, e.g. *_control_lambda.bdg from
#'     MACSv2. REQUIRED
#' @param sfactor Scaling factor for treatment and control track. Keep
#'     it as 1.0 or default in most cases. Set it ONLY while you have
#'     SPMR output from MACS3 callpeak, and plan to calculate scores
#'     as MACS3 callpeak module. If you want to simulate 'callpeak'
#'     w/o '--to-large', calculate effective smaller sample size after
#'     filtering redudant reads in million (e.g., put 31.415926 if
#'     effective reads are 31,415,926) and input it for '-S'; for
#'     'callpeak --to-large', calculate effective reads in larger
#'     sample. DEFAULT: 1.0
#' @param pseudocount The pseudocount used for calculating logLR,
#'     logFE or FE. The count will be applied after normalization of
#'     sequencing depth. DEFAULT: 0.0, no pseudocount is applied.
#' @param method Method to use while calculating a score in any bin by comparing treatment value and control value. Available choices are: ppois, qpois, subtract, logFE, FE, logLR, slogLR, and max. They represent Poisson Pvalue (-log10(pvalue) form) using control as lambda and treatment as observation, q-value through a BH process for poisson pvalues, subtraction from treatment, linear scale fold enrichment, log10 fold enrichment(need to set pseudocount), log10 likelihood between ChIP-enriched model and open chromatin model(need to set pseudocount), symmetric log10 likelihood between two ChIP-enrichment models, or maximum value between the two tracks. Default option is ppois.
#' @param oprefix The PREFIX of output bedGraph file to write
#'     scores. If it is given as A, and method is 'ppois', output file
#'     will be A_ppois.bdg. Mutually exclusive with -o/--ofile.
#' @param outputfile Output filename. Mutually exclusive with
#'     --o-prefix. The number and the order of arguments for --ofile
#'     must be the same as for -m.
#' @param outdir The output directory.
#' @param log Whether to capture logs.
#' @param verbose Set verbose level. 0: only show critical message, 1:
#'     show additional warning message, 2: show process information,
#'     3: show debug messages. If you want to know where are the
#'     duplicate reads, use 3. DEFAULT:2
#' @return `macsList` object.
#' @export
#' @examples
#' eh <- ExperimentHub::ExperimentHub()
#' CHIP <- eh[["EH4558"]]
#' CTRL <- eh[["EH4563"]]
#' p1 <- pileup(CHIP, outdir = tempdir(),
#'              outputfile = "pileup_ChIP_bed.bdg", format = "BED")
#' p2 <- pileup(CTRL, outdir = tempdir(),
#'              outputfile = "pileup_CTRL_bed.bdg", format = "BED")
#' c1 <- bdgcmp(p1$outputs, p2$outputs, outdir = tempdir(),
#'              oprefix = "bdgcmp", pseudocount = 1, method = "FE")
bdgcmp <- function(tfile, cfile, sfactor = 1.0, pseudocount = 0.0,
                   method = c("ppois", "qpois", "subtract", "logFE", "FE", "logLR", "slogLR", "max"),
                   oprefix = character(), outputfile = list(),
                   outdir = ".", log = TRUE, verbose = 2L){
    method <- lapply(method, function(x)match.arg(x, method))
    tf <- normalizePath(tfile)
    cf <- normalizePath(cfile)
    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.namespace, outdir){
        opts <- .namespace()$Namespace(tfile = tf,
                                       cfile = cf,
                                       sfactor = sfactor,
                                       pseudocount = pseudocount,
                                       method = method,
                                       oprefix = oprefix,
                                       ofile = outputfile,
                                       outdir = outdir,
                                       verbose = verbose)
        .bdgcmp <- reticulate::import("MACS3.Commands.bdgcmp_cmd")
        if(log){
            reticulate::py_capture_output(.bdgcmp$run(opts))
        }else{
            .bdgcmp$run(opts)
        }
    }, .namespace = .namespace, outdir = outdir)
    if(log){
        message(res)
    }

    if(length(oprefix) > 0){
        ofile <- file.path(outdir, paste0(oprefix, "_", method, ".bdg"))
    }else{
        ofile <- file.path(outdir, outputfile)
    }

    args <- as.list(match.call())
    macsList(arguments = args, outputs = ofile, log = res)
}
