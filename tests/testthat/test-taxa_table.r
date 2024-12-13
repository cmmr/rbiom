test_that("taxa_table", {
  
  expect_silent(taxa_table(hmp5, transform = "percent"))
  expect_silent(taxa_table(hmp5, transform = "rank"))
  expect_silent(taxa_table(hmp5, transform = "log1p"))
  expect_silent(taxa_table(hmp5, taxa = 0.01))
  expect_error(taxa_table(hmp5, taxa = FALSE))
  expect_error(taxa_table(hmp5, taxa = 'doesnotexist'))
  
  expect_silent(subset_taxa(hmp5, Phylum == 'Bacteroidetes'))
  
})
