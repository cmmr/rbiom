test_that("bdiv_table", {
  
  expect_equal(   # UniFrac
    object    = round(bdiv_table(hmp5, 'UniFrac', weighted = FALSE)$.distance, 2), 
    expected  = c(0.43, 0.46, 0.47, 0.37, 0.33, 0.35, 0.35, 0.27, 0.3, 0.31) )
  expect_equal(
    object    = round(bdiv_table(hmp5, 'UniFrac', weighted = TRUE, normalized = TRUE)$.distance, 2), 
    expected  = c(0.23, 0.43, 0.52, 0.17, 0.24, 0.39, 0.3, 0.36, 0.48, 0.5) )
  expect_equal(
    object    = round(bdiv_table(hmp5, 'UniFrac', weighted = TRUE, normalized = FALSE)$.distance, 2), 
    expected  = c(0.2, 0.37, 0.42, 0.14, 0.21, 0.32, 0.25, 0.3, 0.41, 0.4) )
  expect_error(bdiv_table(min5, 'UniFrac'))
  
  expect_equal(   # Bray-Curtis
    object    = round(bdiv_table(hmp5, 'Bray-Curtis', weighted = FALSE)$.distance, 2), 
    expected  = c(0.44, 0.53, 0.52, 0.4, 0.35, 0.38, 0.34, 0.3, 0.34, 0.36) )
  expect_equal(
    object    = round(bdiv_table(hmp5, 'Bray-Curtis', weighted = TRUE)$.distance, 2), 
    expected  = c(0.51, 0.72, 0.75, 0.54, 0.4, 0.57, 0.64, 0.51, 0.74, 0.67) )
  
  expect_equal(   # Manhattan
    object    = round(bdiv_table(hmp5, 'Manhattan', weighted = FALSE)$.distance, 2), 
    expected  = c(54, 66, 68, 46, 52, 60, 48, 48, 48, 54) )
  expect_equal(
    object    = round(bdiv_table(hmp5, 'Manhattan', weighted = TRUE)$.distance), 
    expected  = c(1545, 2155, 2679, 3015, 1082, 1858, 3404, 1662, 3892, 3890) )
  
  expect_equal(   # Euclidean
    object    = round(bdiv_table(hmp5, 'Euclidean', weighted = FALSE)$.distance, 2), 
    expected  = c(7.35, 8.12, 8.25, 6.78, 7.21, 7.75, 6.93, 6.93, 6.93, 7.35) )
  expect_equal(
    object    = round(bdiv_table(hmp5, 'Euclidean', weighted = TRUE)$.distance), 
    expected  = c(594, 852, 962, 1647, 292, 506, 2169, 383, 2415, 2478) )
  
  expect_equal(   # Jaccard
    object    = round(bdiv_table(hmp5, 'Jaccard', weighted = FALSE)$.distance, 2), 
    expected  = c(0.61, 0.69, 0.68, 0.57, 0.51, 0.55, 0.51, 0.47, 0.51, 0.53) )
  expect_equal(
    object    = round(bdiv_table(hmp5, 'Jaccard', weighted = TRUE)$.distance, 2), 
    expected  = c(0.68, 0.83, 0.86, 0.7, 0.57, 0.73, 0.78, 0.68, 0.85, 0.8) )
  
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
  
  expect_no_error(pthreads())
})

test_that("bdiv_distmat", {
  dm <- bdiv_distmat(min5)
  expect_equal(attr(dm, 'Labels'), min5$samples)
  expect_equal(
    tolerance = 0.02,
    object    = round(as.vector(dm), 2), 
    expected  = c(0.51, 0.72, 0.75, 0.54, 0.4, 0.57, 0.64, 0.51, 0.74, 0.67) )
})
