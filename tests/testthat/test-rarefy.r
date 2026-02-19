test_that("rarefy", {
  
  skip_on_cran()
  
  mtx <- min5$counts
  
  expect_error(suppressWarnings(mtx_rarefy(mtx, margin = 2L, depth = 1:2)))
  expect_error(suppressWarnings(mtx_rarefy(mtx, margin = 2L, depth = 1.5)))
  expect_error(suppressWarnings(mtx_rarefy(mtx, margin = 2L, seed = 1.5)))
  
  expect_warning(mtx_rarefy(mtx, margin = 1L))
  expect_warning(mtx_rarefy(mtx, margin = 2L))
  
  expect_silent(suggest_rarefy_depth(mtx))
})
