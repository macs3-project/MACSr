#' callvar
#'
#' Call variants in given peak regions from the alignment BAM files.
#'
#' @param peakbed Peak regions in BED format, sorted by
#'     coordinates. REQUIRED.
#' @param tfile ChIP-seq/ATAC-seq treatment file in BAM format,
#'     containing only records in peak regions, sorted by
#'     coordinates. Check instruction on how to make the file using
#'     samtools. REQUIRED.
#' @param cfile Control file in BAM format, containing only records in
#'     peak regions, sorted by coordinates. Check instruction on how
#'     to make the file using samtools.
#' @param GQCutoffHetero Genotype Quality score
#'     (-10log10((L00+L11)/(L01+L00+L11))) cutoff for Heterozygous
#'     allele type. Default:0, or there is no cutoff on GQ.
#' @param GQCutoffHomo Genotype Quality score
#'     (-10log10((L00+L01)/(L01+L00+L11))) cutoff for Homozygous
#'     allele (not the same as reference) type. Default:0, or ther is
#'     no cutoff on GQ.
#' @param Q Only consider bases with quality score greater than this
#'     value. Default: 20, which means Q20 or 0.01 error rate.
#' @param maxDuplicate Maximum duplicated reads allowed per mapping
#'     position, mapping strand and the same CIGAR code. Default:
#'     1. When sequencing depth is high, to set a higher value might
#'     help evaluate the correct allele ratio.
#' @param fermi Option to control when to apply local assembly through
#'     Fermi. By default (set as 'auto'), while SAPPER detects any
#'     INDEL variant in a peak region, it will utilize Fermi to
#'     recover the actual DNA sequences to refine the read
#'     alignments. If set as 'on', Fermi will be always invoked. It
#'     can increase specificity however sensivity and speed will be
#'     significantly lower. If set as 'off', Fermi won't be invoked at
#'     all. If so, speed and sensitivity can be higher but specificity
#'     will be significantly lower. Default: auto
#' @param fermiMinOverlap The minimal overlap for fermi to initially
#'     assemble two reads. Must be between 1 and read length. A longer
#'     fermiMinOverlap is needed while read length is small (e.g. 30
#'     for 36bp read, but 33 for 100bp read may work). Default:30
#' @param top2allelesMinRatio The reads for the top 2 most frequent
#'     alleles (e.g. a ref allele and an alternative allele) at a loci
#'     shouldn't be too few comparing to total reads mapped. The
#'     minimum ratio is set by this optoin. Must be a float between
#'     0.5 and 1.  Default:0.8 which means at least 80%% of reads
#'     contain the top 2 alleles.
#' @param altalleleMinCount The count of the alternative
#'     (non-reference) allele at a loci shouldn't be too few. By
#'     default, we require at least two reads support the alternative
#'     allele. Default:2
#' @param maxAR The maximum Allele-Ratio allowed while calculating
#'     likelihood for allele-specific binding. If we allow higher
#'     maxAR, we may mistakenly assign some homozygous loci as
#'     heterozygous. Default:0.95
#' @param np CPU used for mutliple processing. Please note that,
#'     assigning more CPUs does not guarantee the process being
#'     faster. Creating too many parrallel processes need memory
#'     operations and may negate benefit from multi
#'     processing. Default: 1
#' @param verbose Set verbose level of runtime message. 0: only show
#'     critical message, 1: show additional warning message, 2: show
#'     process information, 3: show debug messages. DEFAULT:2
#' @param outputfile Output VCF file name.
#' @param outdir The output directory.
#' @param log Whether to capture logs.
callvar <- function(peakbed, tfile, cfile,
                    outdir = ".", outputfile = character(),
                    GQCutoffHetero = 0,
                    GQCutoffHomo = 0,
                    Q = 20,
                    maxDuplicate = 1L,
                    fermi = "auto",
                    fermiMinOverlap = 30L,
                    top2allelesMinRatio = 0.8,
                    altalleleMinCount = 2L,
                    maxAR = 0.95,
                    np = 1L,
                    verbose = verbose,
                    log = log){
    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.logging, .namespace, outdir){
        opts <- .namespace()$Namespace(peakbed = peakbed,
                                       tfile = file.path(tfile),
                                       cfile = file.path(cfile),
                                       GQCutoffHetero = GQCutoffHetero,
                                       GQCutoffHomo = GQCutoffHomo,
                                       Q = Q,
                                       maxDuplicate = maxDuplicate,
                                       fermi = fermi,
                                       fermiMinOverlap = fermiMinOverlap,
                                       top2allelesMinRatio = top2allelesMinRatio,
                                       altalleleMinCount = altalleleMinCount,
                                       maxAR = maxAR,
                                       np = np,
                                       verbose = verbose,
                                       ofile = outputfile,
                                       outdir = outdir)
        .callvar <- reticulate::import("MACS3.Commands.callvar_cmd")
        if(log){
            .logging()$run()
            reticulate::py_capture_output(.callvar$run(opts))
        }else{
            .callvar$run(opts)
        }
    }, .logging = .logging, .namespace = .namespace, outdir = outdir)
    if(log){
        message(res)
    }

    ofile <- file.path(outdir, outputfile)
    args <- as.list(match.call())
    macsList(arguments = args, outputs = ofile, log = res)
}
