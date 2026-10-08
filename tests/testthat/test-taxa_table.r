test_that("taxa_table", {
  
  skip_on_cran()
  
  expect_silent(taxa_table(hmp5, transform = "percent"))
  expect_silent(taxa_table(hmp5, transform = "rank"))
  expect_silent(taxa_table(hmp5, transform = "log1p"))
  expect_silent(taxa_table(hmp5, taxa = 0.01))
  expect_error(taxa_table(hmp5, taxa = FALSE))
  expect_error(taxa_table(hmp5, taxa = 'doesnotexist'))
  
  expect_silent(subset_taxa(hmp5, Phylum == 'Bacteroidetes'))

})


test_that("taxa_matrix ranks taxa when a sample has none at the rank", {

  # S3's only OTU has no Genus, so unc = "drop" leaves it an empty column.
  biom <- as_rbiom(list(
    counts = matrix(
      data     = c(5, 3, 0,  2, 4, 0,  1, 1, 9),
      nrow     = 3,
      byrow    = TRUE,
      dimnames = list(c("OTU1", "OTU2", "OTU3"), c("S1", "S2", "S3")) ),
    taxonomy = data.frame(
      .otu   = c("OTU1", "OTU2", "OTU3"),
      Phylum = c("Firmicutes", "Bacteroidetes", "Firmicutes"),
      Genus  = c("Lactobacillus", "Bacteroides", NA) )))

  # Mean shares over S1 and S2: Lactobacillus 8/14, Bacteroides 6/14.
  expect_identical(
    rownames(taxa_matrix(biom, 'Genus', taxa = 1, unc = 'drop')),
    'Lactobacillus' )

  expect_setequal(
    rownames(taxa_matrix(biom, 'Genus', taxa = 0.1, unc = 'drop')),
    c('Lactobacillus', 'Bacteroides') )

  mtx <- taxa_matrix(biom, 'Genus', taxa = 1, other = TRUE, unc = 'drop')
  expect_identical(nrow(mtx), 2L)
  expect_equal(unname(colSums(mtx)), c(7, 7, 0))

})
