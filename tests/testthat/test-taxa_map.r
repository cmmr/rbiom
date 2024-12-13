test_that("taxa_map", {
  
  skip_on_cran()
  
  expect_silent(taxa_map(rare5, unc = 'grouped'))
  expect_silent(taxa_map(rare5, unc = 'drop'))
  expect_silent(taxa_map(rare5, rank = 'Genus', taxa = 4, other = TRUE))
  expect_silent(taxa_map(rare5, rank = 'Genus', taxa = 4, other = FALSE))
  
})
