test_that("bdiv_heatmap", {
  
  expect_silent(bdiv_heatmap(rare5, weighted=c(TRUE,FALSE), tracks=c("Body Site", "Age")))
  
  skip_on_cran()
  
  expect_silent(bdiv_heatmap(rare5, weighted=c(TRUE,FALSE), tracks=c("Body Site", "Age")))
  expect_silent(bdiv_heatmap(rare5[1:2]))
  expect_error(bdiv_heatmap(rare5[1]))
  expect_error(bdiv_heatmap(rare5, title = NA))
})
