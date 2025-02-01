test_that("as_rbiom", {
  old_hmp5 <- as_rbiom(as.list(hmp5), pkg_version = "1")
  expect_equal_rbiom(as_rbiom(hmp5), hmp5)
  expect_equal_rbiom(as_rbiom(old_hmp5), hmp5)
  expect_equal_rbiom(as_rbiom(hmp5$counts), min5)
  expect_equal_rbiom(as_rbiom(as.matrix(hmp5$counts)), min5)
  expect_error(as_rbiom(NULL))
})

test_that("as_rbiom.default", {
  json <- write_biom(hmp5, file = NULL, format = "json")
  expect_equal_rbiom(as_rbiom(json), hmp5)
})

test_that("as_rbiom.phyloseq", {
  skip_on_cran()
  
  setClass('phyloseq', slots = c(
    otu_table = 'matrix', sam_data = 'NULL', tax_table = 'NULL', refseq = 'NULL', phy_tree = 'NULL' ))
  expect_s3_class(as_rbiom(new('phyloseq', otu_table = as.matrix(hmp5$counts))), 'rbiom')
})

test_that("as_rbiom.SummarizedExperiment", {
  skip_on_cran()
  
  setClass('SimpleAssaysZZ', slots = c(data = 'list'))
  setClass('SummarizedExperiment', slots = c(assays = 'SimpleAssaysZZ', colData = 'NULL', elementMetadata = 'NULL', NAMES = 'NULL' ))
  expect_s3_class(as_rbiom(new('SummarizedExperiment', assays = new('SimpleAssaysZZ', data = list(as.matrix(hmp5$counts))))), 'rbiom')
})

test_that("as_rbiom.TreeSummarizedExperiment", {
  skip_on_cran()
  
  setClass('SimpleAssaysZZ', slots = c(data = 'list'))
  setClass('PartitioningByEndZZ', slots = c(NAMES = 'NULL'))
  setClass('CompressedGRangesListZZ', slots = c(elementMetadata = 'NULL', partitioning = 'PartitioningByEndZZ'))
  setClass('TreeSummarizedExperiment', slots = c(
    assays = 'SimpleAssaysZZ', colData = 'NULL', rowRanges = 'CompressedGRangesListZZ', rowTree = 'NULL', referenceSeq = 'NULL' ))
  expect_s3_class(as_rbiom(new('TreeSummarizedExperiment', assays = new('SimpleAssaysZZ', data = list(as.matrix(hmp5$counts))))), 'rbiom')
})
