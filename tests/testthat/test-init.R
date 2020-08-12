test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

library(reticulate)
use_condaenv("MACS")
filterdup(ifile = list("tests/testthat/CTCF_PE_CTRL_chr22_50k.bedpe.gz"), outputfile = "test1.bed", outdir = "/tmp", keepduplicates = "1")
