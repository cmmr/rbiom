test_that("bdply", {
  
  expect_error(bdply(min5, NULL, NULL))
  expect_error(bdply(NULL, 'Sex', `$`, 'n_samples'))
  
  expect_equal(bdply(min5,  NULL,  `$`, 'n_samples')$value, 5)
  expect_equal(bdply(hmp50, "Sex", `$`, 'n_samples')$V1, c(30, 20))
  
  skip_on_cran()
  
  
  iters <- list(w = c(TRUE, FALSE), d = c("bray", "euclid"))
  fun   <- function (b, w, d) {
    r <- range(bdiv_distmat(biom = b, bdiv = d, weighted = w))
    round(data.frame(min = r[[1]], max = r[[2]]))
  }
  
  x <- bdply(hmp50, NULL, iters = iters, FUN = fun)
  expect_equal(x, as_rbiom_tbl(tibble::tibble(
    'w'   = c(T, F, T, F),
    'd'   = c('bray','bray','euclid','euclid'),
    'min' = c(0, 0, 123, 3),
    'max' = c(1, 1, 19611, 14) )))
  
  x <- bdply(hmp50, "Sex", iters = iters, prefix = TRUE, FUN = fun)
  expect_equal(x, as_rbiom_tbl(tibble::tibble(
    'Sex' = factor(rep(c('Female', 'Male'), each = 4)),
    '.w'  = c(T, F, T, F, T, F, T, F),
    '.d'  = c('bray','bray','euclid','euclid','bray','bray','euclid','euclid'),
    'min' = c(0, 0, 123, 3, 0, 0, 292, 5),
    'max' = c(1, 1, 17185, 12, 1, 1, 11855, 14) )))
  
})


test_that("blply", {
  
  expect_error(blply(hmp50, NULL, NULL))
  expect_error(blply(NULL, 'Sex', `$`, 'n_samples'))
  
  
  x <- blply(hmp50, NULL, `$`, 'n_samples')
  expect_identical(x, list(50L))
  
  x <- blply(hmp50, "Sex", `$`, 'n_samples') %>% unlist()
  expect_identical(x, c(Female = 30L, Male = 20L))
  
  skip_on_cran()
  
  
  iters <- list(w = c(TRUE, FALSE), d = c("bray", "euclid"))
  fun   <- function (b, w, d) {
    r <- range(bdiv_distmat(biom = b, bdiv = d, weighted = w))
    round(data.frame(min = r[[1]], max = r[[2]]))
  }
  
  
  x <- blply(hmp50, NULL, iters = iters, FUN = fun)
  expect_equal(unlist(x), c(
    `1.min`=0,   `1.max`=1,     `2.min`=0, `2.max`=1, 
    `3.min`=123, `3.max`=19611, `4.min`=3, `4.max`=14 ))
  expect_equal(attr(x, 'split_labels')$w, c(T, F, T, F))
  expect_equal(attr(x, 'split_labels')$d, c('bray','bray','euclid','euclid'))
  
  x <- blply(hmp50, "Sex", iters = iters, prefix = TRUE, FUN = fun)
  expect_equal(unlist(x), c(
    Female.1.min=0,   Female.1.max=1,     Female.2.min=0, Female.2.max=1, 
    Female.3.min=123, Female.3.max=17185, Female.4.min=3, Female.4.max=12, 
    Male.1.min=0,     Male.1.max=1,       Male.2.min=0,   Male.2.max=1, 
    Male.3.min=292,   Male.3.max=11855,   Male.4.min=5,   Male.4.max=14 ))
  expect_equal(attr(x, 'split_labels')$Sex, factor(rep(c('Female', 'Male'), each = 4)))
  expect_equal(attr(x, 'split_labels')$`.w`, rep(c(T, F), 4))
  expect_equal(attr(x, 'split_labels')$`.d`, rep(c('bray','bray','euclid','euclid'), 2))
  
})
