test_that("write_xlsx", {
  
  skip_on_cran()
  
  tfile <- tempfile(fileext = '.xlsx')
  on.exit(unlink(tfile), add = TRUE)
  
  biom <- hmp5$clone()
  attr(biom, 'mtcars') <- head(mtcars)
  
  expect_silent(write_xlsx(biom, file = tfile))
  unlink(tfile)
  
  expect_silent(write_xlsx(rare5, file = tfile))
  
})
