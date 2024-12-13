
test_that("adiv_matrix", {
  df <- as.data.frame(adiv_matrix(min5))
  expect_equal(df$Depth, c(1660, 1371, 1353, 1895, 3939))
  expect_equal(df$OTUs,  c(49, 75, 75, 83, 67))
  expect_equal(round(df$Shannon, 3), c(1.741, 2.587, 2.951, 3.256, 1.463))
  expect_equal(round(df$Chao1, 3), c(59.083, 104.389,  90.042, 93.321, 127.750))
  expect_equal(round(df$Simpson, 3), c(0.567, 0.813, 0.894, 0.932, 0.525))
  expect_equal(round(df$InvSimpson, 3), c(2.309, 5.359, 9.404, 14.785, 2.106))
})

test_that("adiv_table", {
  df <- adiv_table(min5, trans = "rank")
  df <- adiv_table(min5, trans = "rank") # tests caching
  expect_equal(as.character(df$.sample), min5$samples)
  expect_equal(df$.depth, c(1660, 1371, 1353, 1895, 3939))
  expect_equal(df$.adiv, factor(rep('Shannon', 5)))
  expect_equal(df$.diversity, c(2, 3, 4, 5, 1))
  
  df <- adiv_table(min5, adiv = c("OTUs", "InvSimpson"))
  expect_equal(nrow(df), 10)
  expect_equal(round(df$.diversity), c(49, 75, 75, 83, 67, 2, 5, 9, 15, 2))
})

test_that("sample_sums", {
  expect_equal(sample_sums(min5, 0), setNames(c(1660, 1371, 1353, 1895, 3939), min5$samples))
})

test_that("sample_apply", {
  nnz   <- sample_apply(min5, FUN = function (x) sum(x > 0))
  means <- sample_apply(min5, base::mean, sort = 'asc')
  expect_equal(nnz, setNames(c(49, 75, 75, 83, 67), min5$samples))
  expect_equal(unname(round(means)), c(10, 10, 13, 14, 30))
})
