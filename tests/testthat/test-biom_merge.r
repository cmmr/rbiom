test_that("biom_merge", {
  
  actual        <- biom_merge(hmp5[1:2], hmp5[3:5], tree = hmp5$tree)
  actual$counts <- actual$counts[sort(actual$otus), sort(actual$samples)]
  
  expected        <- hmp5$clone()
  expected$counts <- expected$counts[sort(expected$otus), sort(expected$samples)]
  
  expect_equal_rbiom(actual, expected)
  
  skip_on_cran()
  
  expect_error(biom_merge())
  expect_error(biom_merge(NA))
  expect_error(biom_merge(list()))
  expect_error(biom_merge(list(NA, NA)))
  expect_error(biom_merge(min5, min5))
  
  expect_warning(biom_merge(min5[1:2], {x <- min5[3:5]; x$otus %<>% paste0('_'); x}))
  
  expect_equal_rbiom(biom_merge(hmp5), hmp5)
  expect_s3_class(biom_merge(list(min5)), 'rbiom')
  expect_s3_class(biom_merge(system.file("extdata", "hmp50.bz2", package = "rbiom")), 'rbiom')
})
