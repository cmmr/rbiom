test_that("bdply", {
  
  expect_error(bdply(min5, NULL, NULL))
  expect_error(bdply(NULL, 'Sex', `$`, 'n_samples'))
  
  expect_equal(bdply(min5,  NULL,  `$`, 'n_samples')$value, 5)
  expect_equal(bdply(hmp50, "Sex", `$`, 'n_samples')$V1, c(30, 20))
  
  skip_on_cran()
  
  
  iters <- list(tr = c('none', 'log1p'), d = c('bray', 'ham'))
  fun   <- function (b, tr, d) {
    r <- range(bdiv_distmat(biom = b, bdiv = d, transform = tr))
    round(data.frame(min = r[[1]], max = r[[2]]))
  }
  
  x <- bdply(hmp50, NULL, iters = iters, FUN = fun)
  expect_equal(x, as_rbiom_tbl(tibble::tibble(
    'tr'  = c('none','log1p','none','log1p'),
    'd'   = c('bray','bray','ham','ham'),
    'min' = c(0, 0, 9, 2),
    'max' = c(1, 1, 206, 5) )))
  
  x <- bdply(hmp50, 'Sex', iters = iters, prefix = TRUE, FUN = fun)
  expect_equal(x, as_rbiom_tbl(tibble::tibble(
    'Sex' = factor(rep(c('Female', 'Male'), each = 4)),
    '.tr' = c('none', 'log1p', 'none', 'log1p', 'none', 'log1p', 'none', 'log1p'),
    '.d'  = c('bray','bray','ham','ham','bray','bray','ham','ham'),
    'min' = c(0, 0, 9, 2, 0, 0, 26, 3),
    'max' = c(1, 1, 154, 5, 1, 1, 206, 5) )))
  
})


test_that("blply", {
  
  expect_error(blply(hmp50, NULL, NULL))
  expect_error(blply(NULL, 'Sex', `$`, 'n_samples'))
  
  
  x <- blply(hmp50, NULL, `$`, 'n_samples')
  expect_identical(x, list(50L))
  
  x <- blply(hmp50, "Sex", `$`, 'n_samples') %>% unlist()
  expect_identical(x, c(Female = 30L, Male = 20L))
  
  skip_on_cran()
  
  
  iters <- list(tr = c('none', 'log1p'), d = c('bray', 'ham'))
  fun   <- function (b, tr, d) {
    r <- range(bdiv_distmat(biom = b, bdiv = d, transform = tr))
    round(data.frame(min = r[[1]], max = r[[2]]))
  }
  
  
  x <- blply(hmp50, NULL, iters = iters, FUN = fun)
  expect_equal(unlist(x), c(
    `1.min` = 0, `1.max` = 1,   `2.min` = 0, `2.max` = 1, 
    `3.min` = 9, `3.max` = 206, `4.min` = 2, `4.max` = 5 ))
  expect_equal(attr(x, 'split_labels')$tr, c('none','log1p','none','log1p'))
  expect_equal(attr(x, 'split_labels')$d, c('bray','bray','ham','ham'))
  
  x <- blply(hmp50, "Sex", iters = iters, prefix = TRUE, FUN = fun)
  expect_equal(unlist(x), c(
    Female.1.min = 0,  Female.1.max = 1,   Female.2.min = 0, Female.2.max = 1, 
    Female.3.min = 9,  Female.3.max = 154, Female.4.min = 2, Female.4.max = 5, 
    Male.1.min   = 0,  Male.1.max   = 1,   Male.2.min   = 0, Male.2.max   = 1, 
    Male.3.min   = 26, Male.3.max   = 206, Male.4.min   = 3, Male.4.max   = 5 ))
  expect_equal(attr(x, 'split_labels')$Sex, factor(rep(c('Female', 'Male'), each = 4)))
  expect_equal(attr(x, 'split_labels')$`.tr`, rep(c('none','log1p'), 4))
  expect_equal(attr(x, 'split_labels')$`.d`, rep(c('bray','bray','ham','ham'), 2))
  
})
