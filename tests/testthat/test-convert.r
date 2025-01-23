
test_that("convert_to_SE", {
  
  skip_on_cran()
  skip_if_not_installed('SummarizedExperiment')

  x <- suppressMessages(convert_to_SE(hmp5))

  expect_s4_class(x, 'SummarizedExperiment')
  expect_identical(colnames(x), hmp5$samples)
  expect_identical(rownames(x), hmp5$otus)
  expect_equal(x$Age, hmp5$metadata$Age)
})


test_that("convert_to_TSE", {
  
  skip_on_cran()
  skip_if_not_installed('TreeSummarizedExperiment')
  
  rowTree      <- getFromNamespace('rowTree',      'TreeSummarizedExperiment')
  referenceSeq <- getFromNamespace('referenceSeq', 'TreeSummarizedExperiment')
  
  x <- suppressMessages(convert_to_TSE(hmp5))
  
  expect_s4_class(x, 'TreeSummarizedExperiment')
  expect_identical(colnames(x), hmp5$samples)
  expect_identical(rownames(x), hmp5$otus)
  expect_equal(x$Age, hmp5$metadata$Age)
  expect_s3_class(rowTree(x), 'phylo')
  expect_s4_class(referenceSeq(x), 'DNAStringSet')
})
