test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

## library(reticulate)
## use_condaenv("MACS")
## datdir <- system.file("extdata", package = "MACSr")
##filterdup(ifile = file.path(datdir, "CTCF_PE_CTRL_chr22_50k.bedpe.gz"),
##          outputfile = "test.bed", outdir = tempdir(), keepduplicates = "1")
