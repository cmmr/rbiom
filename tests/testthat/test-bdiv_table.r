test_that("bdiv_table", {
  
  expect_equal(   # UniFrac
    object    = round(bdiv_table(hmp5, 'u_unifrac')$.distance, 2), 
    expected  = c(0.43, 0.46, 0.47, 0.37, 0.33, 0.35, 0.35, 0.27, 0.3, 0.31) )
  expect_equal(
    object    = round(bdiv_table(hmp5, 'w_unifrac')$.distance, 2), 
    expected  = c(0.2, 0.37, 0.42, 0.14, 0.21, 0.32, 0.25, 0.3, 0.41, 0.4) )
  expect_equal(
    object    = round(bdiv_table(hmp5, 'n_unifrac')$.distance, 2), 
    expected  = c(0.23, 0.43, 0.52, 0.17, 0.24, 0.39, 0.3, 0.36, 0.48, 0.5) )
  expect_error(bdiv_table(min5, 'UniFrac'))
  
  expect_equal(   # Within
    tolerance = 0.02,
    object    = round(bdiv_table(hmp5, within = 'Sex')$.distance, 2), 
    expected  = c(0.25, 0.40, 0.59, 0.52) )
  expect_equal(
    object   = bdiv_table(hmp5, within = 'Sex')$.distance, 
    expected = bdiv_table(hmp5, md = '==Sex')$.distance )
  expect_error(bdiv_table(hmp5, md = '==Sex', between = 'Sex'))
  expect_error(bdiv_table(hmp5, within = 'Sex', between = 'Sex'))
  
  expect_equal(   # Between
    tolerance = 0.02,
    object    = round(bdiv_table(hmp5, 'jac', between = 'Body Site')$.distance, 2), 
    expected  = c(0.69, 0.68, 0.51, 0.55, 0.51, 0.53) )
  expect_equal(
    object   = bdiv_table(hmp5, 'jac', between = 'Body Site')$.distance, 
    expected = bdiv_table(hmp5, 'jac', md = '!=Body Site')$.distance )
  expect_error(bdiv_table(hmp5, md = '==Body Site', between = 'Body Site'))
  expect_error(bdiv_table(hmp5, within = 'Body Site', between = 'Body Site'))
  expect_error(bdiv_table(hmp5, md = c('==Body Site', '!=Body Site')))
  
  expect_equal(   # Transform
    object   = bdiv_table(hmp5, trans = 'rank', tree = hmp5$tree)$.distance, 
    expected = c(4, 9, 10, 1, 2, 6, 3, 5, 7, 8) )
  expect_equal(
    object   = bdiv_table(hmp5, c('man', 'bray'), trans = 'rank')$.distance, 
    expected = c(4, 9, 10, 1, 2, 6, 3, 5, 7, 8, 4, 9, 10, 1, 2, 6, 3, 5, 7, 8) )
  expect_equal(
    object   = round(bdiv_table(hmp5, 'euc', trans = 'sqrt')$.distance), 
    expected = c(1, 1, 1, 0, 0, 1, 1, 0, 1, 1) )
})

test_that("bdiv_distmat", {
  dm <- bdiv_distmat(min5)
  expect_equal(attr(dm, 'Labels'), min5$samples)
  expect_equal(
    tolerance = 0.02,
    object    = round(as.vector(dm), 2), 
    expected  = c(0.48, 0.69, 0.76, 0.25, 0.4, 0.59, 0.44, 0.52, 0.62, 0.64) )
})
