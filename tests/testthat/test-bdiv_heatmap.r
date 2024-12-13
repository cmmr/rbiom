test_that("bdiv_heatmap", {
  
  p <- expect_silent(bdiv_heatmap(rare5, weighted=c(TRUE,FALSE), tracks=c("Body Site", "Age")))
  
  skip_on_cran()
  
  p <- expect_silent(bdiv_heatmap(rare5, weighted=c(TRUE,FALSE), tracks=c("Body Site", "Age")))
  
  expect_error(bdiv_heatmap(rare5[1]))
  
})
