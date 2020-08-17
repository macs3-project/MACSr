context("test filterdup")

datdir <- system.file("extdata", package = "MACSr")
print(datdir)
ofile <- filterdup(ifile = file.path(datdir, "CTCF_PE_CTRL_chr22_50k.bedpe.gz"),
                   outputfile = "test.bed", outdir = tempdir(), keepduplicates = "1")

test_that("test filterdup", {
    expect_identical(readLines(ofile),
                     readLines(file.path(datdir, "CTCF_PE_CTRL_chr22_50k.filterdup.bed")))
})
