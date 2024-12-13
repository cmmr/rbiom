test_that("taxa_heatmap", {
  
  expect_silent(taxa_heatmap(hmp5, rank = 'Phylum'))
  
  skip_on_cran()
  
  expect_silent(taxa_heatmap(hmp5, rank = 'Phylum', legend = 'bottom', label = 'rows', trees = 'both'))
  expect_silent(taxa_heatmap(hmp5, rank = 'Phylum', label = 'cols', trees = 'none'))
  expect_silent(taxa_heatmap(hmp5, rank = c('Phylum', 'Class'), tracks = c('Sex', 'Age')))
  
})
