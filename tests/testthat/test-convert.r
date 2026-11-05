

test_that("convert_to_animalcules", {

  skip_on_cran()
  skip_if_not_installed('S4Vectors')
  skip_if_not_installed('SummarizedExperiment')
  skip_if_not_installed('MultiAssayExperiment')

  x <- suppressMessages(convert_to_animalcules(hmp5))

  expect_s4_class(x, 'MultiAssayExperiment')
  expect_s3_class(as_rbiom(x), 'rbiom')

  x <- x@ExperimentList[[1]]
  expect_s4_class(x, 'SummarizedExperiment')
  expect_identical(colnames(x), hmp5$samples)
  expect_identical(rownames(x), hmp5$otus)
  expect_equal(x$Age, hmp5$metadata$Age)

  expect_s3_class(as_rbiom(x), 'rbiom')
})


test_that("convert_to_biomformat", {

  skip_on_cran()
  skip_if_not_installed('biomformat')

  x <- suppressMessages(convert_to_biomformat(hmp5))

  expect_s4_class(x, 'biom')
  b <- as_rbiom(x)
  expect_s3_class(b, 'rbiom')
  expect_identical(b$samples, hmp5$samples)
  expect_identical(b$fields, hmp5$fields)
  expect_identical(b$ranks,  hmp5$ranks)

  # Without metadata or taxonomy.
  x <- suppressMessages(convert_to_biomformat(min5))
  expect_s4_class(x, 'biom')
  expect_equal(as.matrix(as_rbiom(x)$counts), as.matrix(min5$counts))
})


test_that("convert_to_phyloseq", {

  skip_on_cran()
  skip_if_not_installed('phyloseq')

  x <- suppressMessages(convert_to_phyloseq(hmp5))

  expect_s4_class(x, 'phyloseq')
  expect_identical(colnames(x@otu_table), hmp5$samples)
  expect_setequal(rownames(x@tax_table), hmp5$otus)
  expect_identical(colnames(x@tax_table), hmp5$ranks[-1])
  expect_equal(x@sam_data$Age, hmp5$metadata$Age)
  expect_true(is.numeric(x@sam_data$Age))
  expect_true("Body Site" %in% colnames(x@sam_data))
  expect_s3_class(x@phy_tree, 'phylo')
  if (requireNamespace('Biostrings', quietly = TRUE))
    expect_s4_class(x@refseq, 'DNAStringSet')

  b <- as_rbiom(x)
  expect_s3_class(b, 'rbiom')
  expect_identical(b$fields, hmp5$fields)
  expect_identical(b$ranks,  hmp5$ranks)
  expect_equal(b$metadata$Age, hmp5$metadata$Age)

  # Without taxonomy, metadata, tree, or sequences.
  x <- suppressMessages(convert_to_phyloseq(min5))
  expect_s4_class(x, 'phyloseq')
  expect_equal(as.matrix(as_rbiom(x)$counts), as.matrix(min5$counts))
})


test_that("convert_to_SE", {

  skip_on_cran()
  skip_if_not_installed('SummarizedExperiment')

  x <- suppressMessages(convert_to_SE(hmp5))

  expect_s4_class(x, 'SummarizedExperiment')
  expect_identical(colnames(x), hmp5$samples)
  expect_identical(rownames(x), hmp5$otus)
  expect_equal(x$Age, hmp5$metadata$Age)
  expect_true("Body Site" %in% names(SummarizedExperiment::colData(x)))

  b <- as_rbiom(x)
  expect_s3_class(b, 'rbiom')
  expect_identical(b$fields, hmp5$fields)
  expect_identical(b$ranks,  hmp5$ranks)

  expect_s3_class(as_rbiom(convert_to_SE(min5)), 'rbiom')
})


test_that("convert_to_TSE", {

  skip_on_cran()
  skip_if_not_installed('TreeSummarizedExperiment')
  skip_if_not_installed('Biostrings')

  x <- suppressMessages(convert_to_TSE(hmp5))

  expect_s4_class(x, 'TreeSummarizedExperiment')
  expect_identical(colnames(x), hmp5$samples)
  expect_identical(rownames(x), hmp5$otus)
  expect_equal(x$Age, hmp5$metadata$Age)
  expect_s3_class(TreeSummarizedExperiment::rowTree(x), 'phylo')
  expect_s4_class(TreeSummarizedExperiment::referenceSeq(x), 'DNAStringSet')

  b <- as_rbiom(x)
  expect_s3_class(b, 'rbiom')
  expect_identical(b$fields, hmp5$fields)
  expect_identical(b$sequences, hmp5$sequences)

  # Without sequences or a tree.
  x <- suppressMessages(convert_to_TSE(min5))
  expect_s4_class(x, 'TreeSummarizedExperiment')
  expect_null(TreeSummarizedExperiment::referenceSeq(x))
  expect_s3_class(as_rbiom(x), 'rbiom')
})
