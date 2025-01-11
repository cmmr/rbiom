test_that("bdiv_table", {
  
  expect_equal(   # UniFrac
    tolerance = 0.02,
    object    = round(bdiv_table(hmp5, 'unifrac', weighted = FALSE)$.distance, 2), 
    expected  = c(0.43, 0.46, 0.47, 0.37, 0.33, 0.35, 0.35, 0.27, 0.3, 0.31) )
  expect_equal(
    tolerance = 0.02,
    object    = round(bdiv_table(hmp5, 'unifrac', weighted = TRUE)$.distance, 2), 
    expected  = c(0.2, 0.37, 0.42, 0.14, 0.21, 0.32, 0.25, 0.3, 0.41, 0.4) )
  expect_error(bdiv_table(min5, 'unifrac'))
  
  expect_equal(   # Within
    tolerance = 0.02,
    object    = round(bdiv_table(hmp5, within = 'Sex', weighted = c(T, F))$.distance, 2), 
    expected  = c(0.54, 0.4, 0.57, 0.51, 0.4, 0.35, 0.38, 0.3) )
  expect_equal(
    object   = bdiv_table(hmp5, within = 'Sex', weighted = c(T, F))$.distance, 
    expected = bdiv_table(hmp5, md = '==Sex', weighted = c(T, F))$.distance )
  expect_error(bdiv_table(hmp5, md = '==Sex', between = 'Sex'))
  expect_error(bdiv_table(hmp5, within = 'Sex', between = 'Sex'))
  
  expect_equal(   # Between
    tolerance = 0.02,
    object    = round(bdiv_table(hmp5, 'jac', between = 'Body Site')$.distance, 2), 
    expected  = c(0.83, 0.86, 0.57, 0.73, 0.85, 0.8) )
  expect_equal(
    object   = bdiv_table(hmp5, 'jac', between = 'Body Site')$.distance, 
    expected = bdiv_table(hmp5, 'jac', md = '!=Body Site')$.distance )
  expect_error(bdiv_table(hmp5, md = '==Body Site', between = 'Body Site'))
  expect_error(bdiv_table(hmp5, within = 'Body Site', between = 'Body Site'))
  expect_error(bdiv_table(hmp5, md = c('==Body Site', '!=Body Site')))
  
  expect_equal(   # Transform
    object   = bdiv_table(hmp5, trans = 'rank', tree = hmp5$tree)$.distance, 
    expected = c(2, 8, 10, 4, 1, 5, 6, 3, 9, 7) )
  expect_equal(
    object   = bdiv_table(hmp5, c('man', 'bray'), trans = 'rank')$.distance, 
    expected = c(2, 5, 6, 7, 1, 4, 8, 3, 10, 9, 2, 8, 10, 4, 1, 5, 6, 3, 9, 7) )
  expect_equal(
    object   = round(bdiv_table(hmp5, 'euc', trans = 'sqrt')$.distance), 
    expected = c(24, 29, 31, 41, 17, 22, 47, 20, 49, 50) )
})

test_that("bdiv_distmat", {
  dm <- bdiv_distmat(min5)
  expect_equal(attr(dm, 'Labels'), min5$samples)
  expect_equal(
    tolerance = 0.02,
    object    = round(as.vector(dm), 2), 
    expected  = c(0.51, 0.72, 0.75, 0.54, 0.4, 0.57, 0.64, 0.51, 0.74, 0.67) )
})
