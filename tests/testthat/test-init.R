context("test filterdup")

datdir <- system.file("extdata", package = "MACSr")
print(datdir)
ofile <- filterdup(ifile = file.path(datdir, "CTCF_SE_ChIP_chr22_50k.bed.gz"),
                   gsize = 5.2e+7, outputfile = "test.bed", outdir = tempdir())

test_that("test filterdup", {
    expect_identical(readLines(ofile),
                     readLines(file.path(datdir, "run_filterdup_result.bed.gz")))
})

context("test predictd")

flen <- predictd(ifile = file.path(datdir, "CTCF_SE_ChIP_chr22_50k.bed.gz"),
                 d_min=10, gsize=5.2e+7, plot = NULL)
test_that("test filterdup", {
    expect_equal(flen, 244)
})
