context("test filterdup")

datdir <- system.file("extdata", package = "MACSr")
eh <- ExperimentHub::ExperimentHub()
eh <- AnnotationHub::query(eh, "MACSdata")
CHIP <- eh[["EH4558"]]
CTRL <- eh[["EH4563"]]
CHIPPE <- eh[["EH4559"]]
CTRLPE <- eh[["EH4564"]]
CHIPBEDPE <- eh[["EH4560"]]
CTRLBEDPE <- eh[["EH4565"]]

fd <- filterdup(ifile = CHIP,
                gsize = 5.2e+7, outputfile = "test.bed", outdir = tempdir(), format = "AUTO")
test_that("test filterdup", {
    expect_identical(readLines(fd$outputs),
                     readLines(file.path(datdir, "run_filterdup_result.bed.gz")))
})

context("test predictd")

flen <- predictd(ifile = CHIP,
                 d_min=10, gsize=5.2e+7, plot = NULL)
test_that("test predictd", {
    expect_equal(flen, c(230, 239))
})

context("test callpeak")

cp1 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
                name = "run_callpeak_narrow0", outdir = tempdir(),
                cutoff_analysis = T, log = FALSE, format = "BED")
cp2 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
                name = "run_callpeak_narrow1", outdir = tempdir(),
                dmin = 15, call_summits = TRUE)
cp3 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
                name = "run_callpeak_narrow2", outdir = tempdir(),
                nomodel = TRUE, extsize = 100)
cp4 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
                name = "run_callpeak_narrow3", outdir = tempdir(),
                nomodel = TRUE, extsize = 100, shift = -50)
cp5 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
                name = "run_callpeak_narrow4", outdir = tempdir(),
                nomodel = TRUE, nolambda = TRUE, extsize = 100,
                shift = -50)
cp6 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
                name = "run_callpeak_narrow5", outdir = tempdir(),
                scaleto = "large")
test_that("test callpeak narrow", {
    expect_identical(readLines(grep("narrowPeak", cp1$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_narrow0_peaks.narrowPeak")))
    expect_identical(readLines(grep("narrowPeak", cp2$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_narrow1_peaks.narrowPeak")))
    expect_identical(readLines(grep("narrowPeak", cp3$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_narrow2_peaks.narrowPeak")))
    expect_identical(readLines(grep("narrowPeak", cp4$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_narrow3_peaks.narrowPeak")))
    expect_identical(readLines(grep("narrowPeak", cp5$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_narrow4_peaks.narrowPeak")))
    expect_identical(readLines(grep("narrowPeak", cp6$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_narrow5_peaks.narrowPeak")))
})

cp7 <- callpeak(CHIP, CTRL, gsize = 5.2e7, store_bdg = TRUE,
                name = "run_callpeak_broad", outdir = tempdir(),
                broad = TRUE)
test_that("test callpeak broad", {
    expect_identical(readLines(grep("broadPeak", cp7$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_broad_peaks.broadPeak")))
})

cp8 <- callpeak(CHIPPE, CTRLPE, gsize = 5.2e7, store_bdg = TRUE,
                format = "BAMPE", name = "run_callpeak_bampe_narrow", outdir = tempdir(),
                call_summits = TRUE)
cp9 <- callpeak(CHIPBEDPE, CTRLBEDPE, gsize = 5.2e7, store_bdg = TRUE,
                format = "BEDPE", name = "run_callpeak_bedpe_narrow", outdir = tempdir(),
                call_summits = TRUE)
cp10 <- callpeak(CHIPBEDPE, gsize = 5.2e7, store_bdg = TRUE,
                format = "BEDPE", name = "run_callpeak_pe_narrow_onlychip", outdir = tempdir())
test_that("test callpeak on PE narrow", {
    expect_identical(readLines(grep("narrowPeak", cp8$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_bampe_narrow_peaks.narrowPeak")))
    expect_identical(readLines(grep("narrowPeak", cp9$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_bedpe_narrow_peaks.narrowPeak")))
    expect_identical(readLines(grep("narrowPeak", cp10$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_pe_narrow_onlychip_peaks.narrowPeak")))
})

cp11 <- callpeak(CHIPPE, CTRLPE, gsize = 5.2e7, store_bdg = TRUE,
                format = "BAMPE", name = "run_callpeak_bampe_broad", outdir = tempdir(),
                broad = TRUE)
cp12 <- callpeak(CHIPBEDPE, CTRLBEDPE, gsize = 5.2e7, store_bdg = TRUE,
                format = "BEDPE", name = "run_callpeak_bedpe_broad", outdir = tempdir(),
                broad = TRUE)
test_that("test callpeak on PE broad", {
    expect_identical(readLines(grep("broadPeak", cp11$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_bampe_broad_peaks.broadPeak")))
    expect_identical(readLines(grep("broadPeak", cp12$outputs, value = TRUE)),
                     readLines(file.path(datdir, "run_callpeak_bedpe_broad_peaks.broadPeak")))
})

