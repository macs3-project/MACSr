
datdir <- system.file("extdata", package = "MACSr")
eh <- ExperimentHub::ExperimentHub()
eh <- AnnotationHub::query(eh, "MACSdata")
CHIP <- eh[["EH4558"]]
CTRL <- eh[["EH4563"]]
CHIPPE <- eh[["EH4559"]]
CTRLPE <- eh[["EH4564"]]
CHIPBEDPE <- eh[["EH4560"]]
CTRLBEDPE <- eh[["EH4565"]]


# context("1. test callpeak")

# cp1 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
#                 name = "run_callpeak_narrow0", outdir = tempdir(),
#                 cutoff_analysis = TRUE, log = FALSE, format = "BED")
# cp2 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
#                 name = "run_callpeak_narrow1", outdir = tempdir(),
#                 dmin = 15, call_summits = TRUE)
# cp3 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
#                 name = "run_callpeak_narrow2", outdir = tempdir(),
#                 nomodel = TRUE, extsize = 100)
# cp4 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
#                 name = "run_callpeak_narrow3", outdir = tempdir(),
#                 nomodel = TRUE, extsize = 100, shift = -50)
# cp5 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
#                 name = "run_callpeak_narrow4", outdir = tempdir(),
#                 nomodel = TRUE, nolambda = TRUE, extsize = 100,
#                 shift = -50)
# cp6 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
#                 name = "run_callpeak_narrow5", outdir = tempdir(),
#                 scaleto = "large")
# test_that("test callpeak narrow", {
#     expect_identical(readLines(grep("narrowPeak", cp1$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_narrow0_peaks.narrowPeak")))
#     expect_identical(readLines(grep("narrowPeak", cp2$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_narrow1_peaks.narrowPeak")))
#     expect_identical(readLines(grep("narrowPeak", cp3$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_narrow2_peaks.narrowPeak")))
#     expect_identical(readLines(grep("narrowPeak", cp4$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_narrow3_peaks.narrowPeak")))
#     expect_identical(readLines(grep("narrowPeak", cp5$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_narrow4_peaks.narrowPeak")))
#     expect_identical(readLines(grep("narrowPeak", cp6$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_narrow5_peaks.narrowPeak")))
# })

# cp7 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
#                 name = "run_callpeak_broad", outdir = tempdir(),
#                 broad = TRUE)
# test_that("test callpeak broad", {
#     expect_identical(readLines(grep("broadPeak", cp7$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_broad_peaks.broadPeak")))
# })

# cp8 <- callpeak(CHIPPE, CTRLPE, gsize = 5.2e7, store_bdg = TRUE,
#                 format = "BAMPE", name = "run_callpeak_bampe_narrow", outdir = tempdir(),
#                 call_summits = TRUE)
# cp9 <- callpeak(CHIPBEDPE, CTRLBEDPE, gsize = 5.2e7, store_bdg = TRUE,
#                 format = "BEDPE", name = "run_callpeak_bedpe_narrow", outdir = tempdir(),
#                 call_summits = TRUE)
# cp10 <- callpeak(CHIPBEDPE, gsize = 5.2e7, store_bdg = TRUE,
#                 format = "BEDPE", name = "run_callpeak_pe_narrow_onlychip", outdir = tempdir())
# test_that("test callpeak on PE narrow", {
#     expect_identical(readLines(grep("narrowPeak", cp8$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_bampe_narrow_peaks.narrowPeak")))
#     expect_identical(readLines(grep("narrowPeak", cp9$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_bedpe_narrow_peaks.narrowPeak")))
#     expect_identical(readLines(grep("narrowPeak", cp10$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_pe_narrow_onlychip_peaks.narrowPeak")))
# })

# cp11 <- callpeak(CHIPPE, CTRLPE, gsize = 5.2e7, store_bdg = TRUE,
#                 format = "BAMPE", name = "run_callpeak_bampe_broad", outdir = tempdir(),
#                 broad = TRUE)
# cp12 <- callpeak(CHIPBEDPE, CTRLBEDPE, gsize = 5.2e7, store_bdg = TRUE,
#                 format = "BEDPE", name = "run_callpeak_bedpe_broad", outdir = tempdir(),
#                 broad = TRUE)
# test_that("test callpeak on PE broad", {
#     expect_identical(readLines(grep("broadPeak", cp11$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_bampe_broad_peaks.broadPeak")))
#     expect_identical(readLines(grep("broadPeak", cp12$outputs, value = TRUE)),
#                      readLines(file.path(datdir, "run_callpeak_bedpe_broad_peaks.broadPeak")))
# })

# context("2. test predictd")
# flen <- predictd(ifile = CHIP,
#                  d_min=10, gsize=5.2e+7, plot = NULL)
# test_that("test predictd", {
#     expect_equal(flen, 229)
# })


# context("3. test filterdup")
# fd <- filterdup(ifile = CHIP,
#                 gsize = 5.2e+7, outputfile = "test.bed", outdir = tempdir(), format = "AUTO")
# test_that("test filterdup", {
#     expect_identical(readLines(fd$outputs),
#                      readLines(file.path(datdir, "run_filterdup_result.bed.gz")))
# })

## ## test callvar
## if(file.exists(system.file(package="Rsamtools"))){
##     Rsamtools::indexBam(CHIPPE)
##     Rsamtools::indexBam(CTRLPE)
##     callvarpeak <- file.path(datdir, "callvar_testing.narrowPeak")
##     cv1 <- callvar(peakbed=callvarpeak, tfile=CHIPPE, cfile=CTRLPE, outputfile="PEsample.vcf")
##     test_that("callvar", {
##         expect_equal(nrow(read.table(cv1$outputs, sep="\t")), 16)
##     })
## }


## atac <- hmmratac(bam=ATACSEQBAM, outdir="/tmp/atac", name="hmmratac_yeast500k", save_train=TRUE)
# datdir <- system.file("extdata", package = "MACSr")

# CHIP <- file.path(datdir, "CTCF_SE_ChIP_chr22_50k.bed.gz")
# CTRL <- file.path(datdir, "CTCF_SE_CTRL_chr22_50k.bed.gz")

# CHIPPE <- file.path(datdir, "CTCF_PE_ChIP_chr22_50k.bam")
# CTRLPE <- file.path(datdir, "CTCF_PE_CTRL_chr22_50k.bam")

# CHIPBEDPE <- file.path(datdir, "CTCF_PE_ChIP_chr22_50k.bedpe.gz")
# CTRLBEDPE <- file.path(datdir, "CTCF_PE_CTRL_chr22_50k.bedpe.gz")


context("1. test callpeak")
    test_that("callpeak runs without error", {
        expect_no_error(cp1 <- callpeak(
            tfile = CHIP, 
            cfile = CTRL, 
            gsize = 5.2e7, 
            format = "BED", 
            name = "callpeak_output"))
    })

context("2. test pileup")
    test_that("pileup runs without error", {
        expect_no_error(pu1 <- pileup( ifile = CHIP,
            format = 'BED', 
            extsize = 200,
            outputfile = "pileup_output.bed.bdg"))
    })


context("3. test filterdup")
    test_that("filterdup runs without error", {
        expect_no_error(fd1 <- filterdup( ifile = CHIP,
            gsize = 52000000, 
            outputfile = "filterdup_output.bed"))
    })
    test_that("filterdup runs without error, format = bedpe", {
        expect_no_error(fd2 <- filterdup( ifile = CHIPPE,
            gsize = 52000000, 
            format = 'BAMPE',
            outputfile = "filterdup_output.bedpe"))
    })

context("4. test predictd") # this runs successfully, but has mismatch argument (plot vs rfile?)
    test_that("predictd runs without error", {
        expect_no_error(pd <- predictd(ifile = CHIP, 
        d_min = 10, 
        gsize=5.2e+7, 
        plot = NULL))
    })

context("5. test randsample")
    test_that("randsample runs without error", {
        expect_no_error(rs <- randsample(ifile = CHIP, 
        seed = 31415926, 
        number = 10000, 
        outputfile = "randsample_output.bed"))
    })

# context("6. test refinepeak")
#     test_that("refinepeak runs without error 1", {
#         expect_no_error(rp1 <- refinepeak(
#             bedfile = file.path(datdir, 'run_callpeak_narrow0_peaks.narrowPeak'), 
#             ifile = CHIP, 
#             oprefix="refinepeak_output_prefix"))
#     })
#     test_that("refinepeak runs without error 2", {
#         expect_no_error(rp2 <- refinepeak(
#             bedfile = file.path(datdir, 'run_callpeak_narrow0_peaks.narrowPeak'), 
#             ifile = CHIP, 
#             outputfile="refinepeak_output.bed"))
#     })

context("7. test bdgcmp")
    test_that("bdgcmp runs without error 1", {
    expect_no_error(bdgcmp1 <- bdgcmp(
        tfile = file.path(datdir, 'run_pileup_ChIP.bed.bdg'), 
        cfile =file.path(datdir, 'run_pileup_CTRL.bed.bdg'),
        method='ppois', 
        pseudocount=1, 
        oprefix="bdgcmp_output_prefix" ))
    })
    test_that("bdgcmp runs without error 2", {
    expect_no_error(bdgcmp2 <- bdgcmp(
        tfile = file.path(datdir, 'run_pileup_ChIP.bed.bdg'), 
        cfile =file.path(datdir, 'run_pileup_CTRL.bed.bdg'),
        method='ppois', 
        pseudocount=1, 
        outputfile= list("bdgcmp_output.bdg") )) # why is it defined as a list?
    })

# context("8. test bdgpeakcall")
#     test_that("bdgpeakcall runs without error 1", {
#         expect_no_error(bdgpc1 <- bdgpeakcall(
#             ifile = file.path(datdir, 'run_bdgcmp_FE.bdg'), 
#             cutoff = 2,
#             oprefix="bdgpeakcall_output_prefix"))
#     })
#     test_that("bdgpeakcall runs without error 2", {
#         expect_no_error(bdgpc2 <- bdgpeakcall(
#             ifile = file.path(datdir, 'run_bdgcmp_FE.bdg'), 
#             cutoff = 2,
#             outputfile="bdgpeakcall_output.txt"))
#     })

# context("9. test bdgbroadcall")
#     test_that("bdgbroadcall runs without error 1", {
#         expect_no_error(bdgbc1 <- bdgbroadcall(
#             ifile = file.path(datdir, 'run_bdgcmp_FE.bdg'), 
#             cutoffpeak=2, 
#             cutofflink=1.5,
#             oprefix="bdgbroadcall_output" ))
#     })
#     test_that("bdgbroadcall runs without error 2", {
#         expect_no_error(bdgbc2 <- bdgbroadcall(
#             ifile = file.path(datdir, 'run_bdgcmp_FE.bdg'), 
#             cutoffpeak=2, 
#             cutofflink=1.5,
#             outputfile="bdgbroadcall_output.bed" ))
#     })

# context("10. test bdgdiff")
#     test_that("bdgdiff runs without error 1", {
#         expect_no_error(bdgdiff1 <- bdgdiff(
#             t1bdg = file.path(datdir, 'run_callpeak_narrow0_treat_pileup.bdg'), 
#             c1bdg = file.path(datdir, 'run_callpeak_narrow0_control_lambda.bdg'), 
#             t2bdg = file.path(datdir, 'run_callpeak_narrow_revert_treat_pileup.bdg'), 
#             c2bdg = file.path(datdir, 'run_callpeak_narrow_revert_control_lambda.bdg'),
#             oprefix="bdgdiff_output_prefix"))
#         })
#     test_that("bdgdiff runs without error 2", {
#         expect_no_error(bdgdiff2 <- bdgdiff(
#             t1bdg = file.path(datdir, 'run_callpeak_narrow0_treat_pileup.bdg'), 
#             c1bdg = file.path(datdir, 'run_callpeak_narrow0_control_lambda.bdg'), 
#             t2bdg = file.path(datdir, 'run_callpeak_narrow_revert_treat_pileup.bdg'), 
#             c2bdg = file.path(datdir, 'run_callpeak_narrow_revert_control_lambda.bdg'),
#             outputfile= list("bdgdiff_output1.bed", "bdgdiff_output2.bed", "bdgdiff_output3.bed")))
#         })

# context("11. test cmbreps")
#      test_that("cmbreps runs without error", {
#         expect_no_error(cmbreps <- cmbreps(
#             ifile = list(file.path(datdir, 'run_callpeak_narrow0_treat_pileup.bdg'),
#                     file.path(datdir, 'run_callpeak_narrow0_control_lambda.bdg')),
#             method = "max",
#             outputfile = "cmbreps_output.bdg"
#             ))
#     })

context("12. test bdgopt")
    test_that("bdgopt runs without error", {
        expect_no_error(bdgopt <- bdgopt(
            ifile = file.path(datdir, 'run_callpeak_narrow0_treat_pileup.bdg'), 
            method = "min", 
            extraparam = 10,
            outputfile="bdgopt_output.bdg" ))
    })

# context("13. test callvar")
#     test_that("callvar runs without error", {
#         expect_no_error(cv1 <- callvar(
#             peakbed=file.path(datdir, "callvar_testing.narrowPeak"), 
#             tfile=CHIPPE, 
#             cfile=CTRLPE,
#             outputfile="callvar_output.vcf"
#             ))
#     })
    
# context("14. test 50k contigs with buffersize")

# context("15. test hmmratac")
#     test_that("hmmratac runs without error, defaults", {
#         expect_no_error(hm1 <- hmmratac(
#             input_file = file.path(datdir, 'yeast_500k_SRR1822137.bam'),
#             name = "hmmratac_output"))
#     })
#     test_that("hmmratac runs without error, format = BEDPE", {
#         expect_no_error(hm2 <- hmmratac(
#             input_file = file.path(datdir, 'yeast_500k_SRR1822137.bedpe.gz'), 
#             format = 'BEDPE',
#             name = "hmmratac_bedpe_output"))
#     })


