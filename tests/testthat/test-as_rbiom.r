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
  phy_file <- system.file("extdata", "rich_sparse_otu_table.biom", package="phyloseq")
  phy      <- phyloseq::import_biom(phy_file)
  phy_mtx  <- structure(attr(phy, 'otu_table'), class=NULL, taxa_are_rows = NULL)
  biom_mtx <- as.matrix(as_rbiom(phy)$counts)
  expect_equal(biom_mtx, phy_mtx)
  expect_error(as_rbiom(phy, no_phyloseq = TRUE))
})
