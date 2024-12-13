test_that("rarefy", {
  
  skip_on_cran()
  
  mtx <- min5$counts
  
  expect_error(rarefy_cols(mtx, depth = NULL))
  expect_error(rarefy_cols(mtx, depth = 1:2))
  expect_error(rarefy_cols(mtx, depth = 1.5))
  expect_error(rarefy_cols(mtx, n = 1:2))
  expect_error(rarefy_cols(mtx, n = 1.5))
  expect_error(rarefy_cols(mtx, seed = 1.5))
  
  expect_silent(rarefy_cols(mtx, n=0))
  expect_silent(rarefy_cols(mtx, n=0.5))
  expect_silent(rarefy_cols(mtx, n=-2))
  
  expect_silent(rescale_rows(mtx))
  expect_silent(rescale_rows(min5))
  
  expect_silent(rare_suggest(mtx))
  
  mtx <- expect_silent(rescale_cols(min5))
  expect_silent(rarefy_cols(mtx))
})
