
# Evaluate `x$code` in a clean environment and return `data`.
replay_data <- function (x) {
  env <- new.env(parent = globalenv())
  eval(parse(text = x$code), envir = env)
  env$data
}


test_that("adiv_matrix", {
  df <- as.data.frame(adiv_matrix(min5))
  expect_equal(df$depth, c(1660, 1371, 1353, 1895, 3939))
  expect_equal(df$observed,  c(49, 75, 75, 83, 67))
  expect_equal(round(df$shannon, 3), c(1.741, 2.587, 2.951, 3.256, 1.463))
  expect_equal(round(df$simpson, 3), c(0.567, 0.813, 0.894, 0.932, 0.525))
})

test_that("adiv_table", {
  df <- adiv_table(min5, adiv = 'shannon', trans = "rank")
  df <- adiv_table(min5, adiv = 'shannon', trans = "rank") # tests caching
  expect_equal(as.character(df$.sample), min5$samples)
  expect_equal(df$.depth, c(1660, 1371, 1353, 1895, 3939))
  expect_equal(df$.adiv, factor(rep('shannon', 5)))
  expect_equal(df$.diversity, c(2, 3, 4, 5, 1))

  df <- adiv_table(min5, adiv = c("observed", "inv_simpson"))
  expect_equal(nrow(df), 10)
  expect_equal(round(df$.diversity), c(49, 75, 75, 83, 67, 2, 5, 9, 15, 2))

  # "percent" is not a diversity transform.
  expect_error(adiv_vector(min5, transform = "percent"))
})

test_that("adiv provenance replays", {

  skip_on_cran()

  b <- rarefy(hmp5, depth = 1000, seed = 3)

  tbl <- adiv_table(b, adiv = c("shannon", "simpson"), md = "Sex")
  expect_s3_class(tbl$code, 'rbiom_code')
  expect_true(grepl("rarefy(biom, depth = 1000L, seed = 3)", tbl$code, fixed = TRUE))
  expect_equal(as.data.frame(replay_data(tbl)), as.data.frame(tbl), ignore_attr = TRUE)

  mtx <- adiv_matrix(b)
  expect_s3_class(mtx, 'rbiom_obj')
  expect_true(is.matrix(mtx))
  expect_equal(unclass(replay_data(mtx)), unclass(mtx), ignore_attr = TRUE)

  vec <- adiv_vector(b, "simpson", transform = "rank")
  expect_equal(unname(replay_data(vec)), unname(vec), ignore_attr = TRUE)

  # Printing hides the provenance attributes.
  expect_output(print(adiv_vector(rare50)), "HMP01")
  expect_no_match(paste(capture.output(print(adiv_vector(rare50))), collapse = "\n"), "cmd")
})

test_that("adiv works with the cache disabled", {
  skip_on_cran()
  rlang::with_options(rbiom.cache_dir = "FALSE", {
    expect_s3_class(adiv_table(min5), 'rbiom_tbl')
    expect_equal(length(adiv_vector(min5)), 5L)
  })
})

test_that("rank transforms do not disturb the RNG", {
  skip_on_cran()
  set.seed(42)
  before <- .Random.seed
  adiv_vector(min5, transform = "rank", seed = 1)
  expect_identical(.Random.seed, before)
  rm(list = '.Random.seed', envir = globalenv())
  adiv_vector(min5, transform = "rank", seed = 2)
  expect_false(exists('.Random.seed', envir = globalenv()))
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
