test_that("taxa_clusters", {
  x <- taxa_clusters(min5, k = 2)
  expect_s3_class(x, 'factor')
  expect_identical(names(x), min5$samples)
  expect_identical(nlevels(x), 2L)
  expect_identical(x, taxa_clusters(min5, k = 2))

  # Seeded, with several random starts: the same seed gives the same partition.
  expect_identical(taxa_clusters(rare50, k = 3, seed = 1), taxa_clusters(rare50, k = 3, seed = 1))
  set.seed(11); before <- .Random.seed
  taxa_clusters(rare50, k = 3, seed = 4)
  expect_identical(.Random.seed, before)
})

test_that("bdiv_clusters", {
  x <- bdiv_clusters(min5, k = 2)
  expect_s3_class(x, 'factor')
  expect_identical(names(x), min5$samples)
  expect_identical(nlevels(x), 2L)
  expect_identical(x, bdiv_clusters(min5, k = 2))

  # k-medoids on the distance matrix itself.
  dm <- bdiv_distmat(min5)
  expect_identical(unname(as.integer(x)), unname(cluster::pam(dm, k = 2)$clustering))
})
