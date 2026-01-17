test_that("rarefy", {
  
  skip_on_cran()
  
  mtx <- min5$counts
  
  expect_error(mtx_rarefy(mtx, margin = 2L, depth = NULL))
  expect_error(mtx_rarefy(mtx, margin = 2L, depth = 1:2))
  expect_error(mtx_rarefy(mtx, margin = 2L, depth = 1.5))
  expect_error(mtx_rarefy(mtx, margin = 2L, n = 1:2))
  expect_error(mtx_rarefy(mtx, margin = 2L, n = 1.5))
  expect_error(mtx_rarefy(mtx, margin = 2L, seed = 1.5))
  
  expect_silent(mtx_rarefy(mtx, margin = 2L, n = 0))
  expect_silent(mtx_rarefy(mtx, margin = 2L, n = 0.5))
  expect_silent(mtx_rarefy(mtx, margin = 2L, n = -2))
  
  expect_silent(mtx_rarefy(mtx, margin = 1L))
  expect_silent(mtx_rarefy(mtx, margin = 2L))
  
  expect_silent(rare_suggest(mtx))
})
