test_that("taxa_clusters", {
  x <- taxa_clusters(min5, k = 2)
  x <- setNames(x == x[[1]], names(x))
  expect_equal(x, setNames(c(T,T,T,T,F), min5$samples))
})

test_that("bdiv_clusters", {
  x <- bdiv_clusters(min5, k = 2)
  x <- setNames(x == x[[1]], names(x))
  expect_equal(x, setNames(c(T,F,F,F,T), min5$samples))
})
