

test_that("convert_to_phyloseq", {
  
  skip_on_cran()
  skip_if_not_installed('phyloseq')
  
  x <- suppressMessages(convert_to_phyloseq(hmp5))
  
  expect_s4_class(x, 'phyloseq')
  expect_identical(colnames(x@otu_table), hmp5$samples)
  expect_setequal(rownames(x@tax_table), hmp5$otus)
  expect_equal(x@sam_data$Age, as.character(hmp5$metadata$Age))
  expect_s3_class(x@phy_tree, 'phylo')
  expect_s4_class(x@refseq, 'DNAStringSet')
  
  expect_s3_class(as_rbiom(x), 'rbiom')
})



test_that("convert_to_SE", {
  
  skip_on_cran()
  skip_if_not_installed('SummarizedExperiment')

  x <- suppressMessages(convert_to_SE(hmp5))

  expect_s4_class(x, 'SummarizedExperiment')
  expect_identical(colnames(x), hmp5$samples)
  expect_identical(rownames(x), hmp5$otus)
  expect_equal(x$Age, hmp5$metadata$Age)
  
  expect_s3_class(as_rbiom(x), 'rbiom')
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
  
  expect_s3_class(as_rbiom(x), 'rbiom')
})
