test_that("taxa_stacked", {
  
  expect_silent(taxa_stacked(rare5))
  
  skip_on_cran()
  
  expect_silent(taxa_stacked(rare5, rank = c("Genus", "Phylum")))
  expect_silent(taxa_stacked(rare5, label.by = 'Body Site'))
  expect_silent(taxa_stacked(rare5, label.by = 'Body Site', order.by = 'Sex'))
  expect_silent(taxa_stacked(hmp5,  patterns = TRUE, order.by = 'Sex'))
  
})
