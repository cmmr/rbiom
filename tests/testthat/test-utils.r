test_that("utils", {

  skip_on_cran()

  expect_error(eval_envir(new.env(), x = 1, y = 1))
  expect_true(eq(structure(NA, hash = 1), structure(NA, hash = 2)))
  expect_silent(drop_cols(data.frame()))
  expect_null(keep_cols(data.frame()))
  expect_identical(coan(NULL), '')
  expect_identical(vw(NULL), '<none>')
  expect_identical(vw('a'), 'a')
  expect_identical(vw(c('a', 'b')), 'a and b')
  expect_error(require_package('doesnotexist', 'test'))
})


test_that("with_seed restores the caller's RNG state", {

  set.seed(42); a <- runif(1)
  set.seed(42)
  x <- with_seed(1, runif(3))
  y <- with_seed(1, runif(3))
  b <- runif(1)

  expect_identical(a, b)
  expect_identical(x, y)

  # Also when no seed existed beforehand.
  rm(list = '.Random.seed', envir = globalenv())
  with_seed(1, runif(1))
  expect_false(exists('.Random.seed', envir = globalenv()))
})


test_that("current_cmd builds commands from evaluated arguments", {

  f <- function (biom, adiv = "shannon", md = ".all", cpus = 3L, ...) current_cmd('adiv_table')
  environment(f) <- asNamespace('rbiom')

  expect_identical(as.character(f(hmp5)), 'adiv_table(biom)')
  expect_identical(as.character(f(hmp5, "simpson", md = "Sex")), 'adiv_table(biom, "simpson", "Sex")')

  # Arguments are never re-evaluated.
  n <- 0
  g <- function () { n <<- n + 1; hmp5 }
  f(g())
  expect_equal(n, 1)

  # `cpus` is omitted; `.dots` are flattened; unknown args are dropped
  # unless the function has `...`.
  h <- function (biom, cpus = 1, ...) NULL
  expect_false(grepl('cpus', as.args(list(biom = hmp5, cpus = 6, .dots = list(z = 1)), fun = h)))
  expect_true(grepl('z = 1',  as.args(list(biom = hmp5, cpus = 6, .dots = list(z = 1)), fun = h)))
  expect_false(grepl('z = 1', as.args(list(biom = hmp5, .dots = list(z = 1)), fun = adiv_table)))
})


test_that("default_value evaluates defaults in context", {

  f <- function (x, y = x, z = n_cpus(), w) NULL

  expect_identical(default_value(f, 'y', list(x = 5)), 5)
  expect_true(is.integer(default_value(f, 'z')))
  expect_s3_class(default_value(f, 'w'), 'rbiom_no_default')
  expect_s3_class(default_value(f, 'x'), 'rbiom_no_default')
  expect_s3_class(default_value(f, 'q'), 'rbiom_no_default')
})


test_that("split_apply and rows_apply", {

  df <- data.frame(g = c('b', 'a', 'b', NA), h = c(1, 1, 2, 2), v = 1:4)

  r <- split_apply(df, 'g', function (d) data.frame(n = nrow(d), s = sum(d$v)))
  expect_equal(r$g, c('a', 'b'))
  expect_equal(r$n, c(1L, 2L))
  expect_equal(r$s, c(2L, 4L))

  r <- split_apply(df, c('g', 'h'), function (d) c(s = sum(d$v)))
  expect_equal(nrow(r), 3L)
  expect_equal(r$s, c(2, 1, 3))

  r <- split_apply(df, 'g', function (d) nrow(d), as = 'list')
  expect_named(r, c('a', 'b'))
  expect_equal(attr(r, 'split_labels')$g, c('a', 'b'))

  r <- split_apply(df, NULL, function (d) data.frame(n = nrow(d)))
  expect_equal(r$n, 4L)

  r <- split_apply(df, 'g', function (d) NULL)
  expect_equal(nrow(r), 0L)

  r <- rows_apply(data.frame(a = 1:2), function (row) data.frame(b = row$a * 2))
  expect_equal(r$b, c(2, 4))

  expect_true(is_empty_df(NULL))
  expect_true(is_empty_df(data.frame()))
  expect_false(is_empty_df(df))
})
