
# A stand-in for a multi-threaded BLAS that rbiom can control, so that
# blas_limit() can be tested anywhere: R's BLAS is single-threaded on
# Windows and CRAN. `fake$sets` records every count that was set.
local_fake_blas <- function (n = 6L, env = parent.frame()) {

  fake <- new.env()
  fake$n    <- n
  fake$sets <- integer(0)

  local_mocked_bindings(
    blas_ctl     = function () TRUE,
    blas_threads = function (n = NULL) {
      if (is.null(n)) return (fake$n)
      fake$n    <- as.integer(n)
      fake$sets <- c(fake$sets, fake$n)
      invisible() },
    .env = env )

  fake
}

# The thread count while a function that called blas_limit(cpus) runs.
count_in <- function (fake, cpus) {
  fun <- function (cpus) { blas_limit(cpus); fake$n }
  fun(cpus)
}


test_that("blas_limit uses one thread while the function runs", {

  fake <- local_fake_blas(6L)

  expect_identical(count_in(fake, 1L), 1L)
  expect_identical(fake$n, 6L)

  expect_identical(count_in(fake, 4), 1L)
  expect_identical(fake$n, 6L)
  expect_identical(fake$sets, c(1L, 6L, 1L, 6L))
})


test_that("blas_limit leaves a count that is within cpus", {

  fake <- local_fake_blas(4L)

  expect_identical(count_in(fake, 4L), 4L)
  expect_identical(count_in(fake, 6L), 4L)

  fake$n <- 1L
  expect_identical(count_in(fake, 1L), 1L)
  expect_identical(fake$sets, integer(0))
})


test_that("blas_limit restores the count after an error", {

  fake <- local_fake_blas(6L)
  fun  <- function () { blas_limit(2L); stop("boom") }

  expect_error(fun(), "boom")
  expect_identical(fake$n, 6L)
  expect_identical(fake$sets, c(1L, 6L))
})


test_that("blas_limit restores the count after an interrupt", {

  fake <- local_fake_blas(6L)
  fun  <- function () { blas_limit(2L); rlang::interrupt() }

  tryCatch(fun(), interrupt = function (e) NULL)
  expect_identical(fake$n, 6L)
  expect_identical(fake$sets, c(1L, 6L))
})


test_that("nested limits restore in order", {

  fake  <- local_fake_blas(6L)
  inner <- function (cpus) { blas_limit(cpus); fake$n }
  outer <- function (cpus) { blas_limit(cpus); c(fake$n, inner(2L), fake$n) }

  expect_identical(outer(6L), c(6L, 1L, 6L))
  expect_identical(outer(3L), c(1L, 1L, 1L))
  expect_identical(fake$n, 6L)
  expect_identical(fake$sets, c(1L, 6L, 1L, 6L))
})


test_that("blas_limit does nothing it can't undo, or can't do", {

  fake <- local_fake_blas(6L)

  # No function frame to restore the count in.
  blas_limit(1L, env = new.env())
  expect_identical(fake$sets, integer(0))

  # Not a usable count; validate_cpus() reports these.
  for (cpus in list(NULL, NA, 0, -1, 1.5, "1", 1:2))
    expect_identical(count_in(fake, cpus), 6L)
  expect_identical(fake$sets, integer(0))

  # The library's count can't be read.
  fake$n <- NA_integer_
  expect_identical(count_in(fake, 1L), NA_integer_)
  expect_identical(fake$sets, integer(0))
})


test_that("blas_limit is silent and does nothing without control", {

  fake <- local_fake_blas(6L)
  local_mocked_bindings(blas_ctl = function () FALSE)

  fun <- function () { blas_limit(1L); fake$n }
  expect_silent(expect_identical(fun(), 6L))
  expect_identical(fake$sets, integer(0))
})


test_that("the real control is silent, whatever ecodive and the BLAS provide", {

  before <- blas_threads()
  if (!blas_ctl()) expect_identical(before, NA_integer_)

  fun <- function () { blas_limit(1L); blas_threads() }
  expect_silent(n <- fun())
  expect_true(is.na(n) || n == 1L)
  expect_identical(blas_threads(), before)

  # Setting the count it already has, or a count that isn't one, changes
  # nothing.
  expect_silent(blas_threads(before))
  for (bad in list(NA, NA_integer_, 0L, -1, 1.5, "2", 1:2))
    expect_silent(blas_threads(bad))
  expect_identical(blas_threads(), before)
})


test_that("every exported function with `cpus` applies it before computing", {

  skip_on_cran()
  skip_if_not_installed('openxlsx')
  rlang::local_options(rbiom.cache_dir = FALSE, lifecycle_verbosity = "quiet")

  counts <- as.matrix(hmp50$counts)
  dm     <- bdiv_distmat(hmp50)
  df     <- adiv_table(hmp50)

  # Stop each call where it limits the threads, and say which function's
  # frame the limit would last for.
  local_mocked_bindings(blas_limit = function (cpus, env = parent.frame()) {
    frame <- which(vapply(sys.frames(), identical, logical(1), env))
    rlang::abort("limited", class = "rbiom_test_blas_limit", cpus = cpus, fun = sys.function(frame[[1]]))
  })

  calls  <- list(
    adiv_boxplot      = function () adiv_boxplot(hmp50, cpus = 1),
    adiv_corrplot     = function () adiv_corrplot(hmp50, x = "Age", cpus = 1),
    adiv_matrix       = function () adiv_matrix(hmp50, cpus = 1),
    adiv_stats        = function () adiv_stats(hmp50, stat.by = "Sex", cpus = 1),
    adiv_table        = function () adiv_table(hmp50, cpus = 1),
    adiv_vector       = function () adiv_vector(hmp50, cpus = 1),
    bdiv_boxplot      = function () bdiv_boxplot(hmp50, cpus = 1),
    bdiv_clusters     = function () bdiv_clusters(hmp50, cpus = 1),
    bdiv_corrplot     = function () bdiv_corrplot(hmp50, x = "Age", cpus = 1),
    bdiv_distmat      = function () bdiv_distmat(hmp50, cpus = 1),
    bdiv_heatmap      = function () bdiv_heatmap(hmp50, cpus = 1),
    bdiv_matrix       = function () bdiv_matrix(hmp50, cpus = 1),
    bdiv_ord_plot     = function () bdiv_ord_plot(hmp50, cpus = 1),
    bdiv_ord_table    = function () bdiv_ord_table(hmp50, cpus = 1),
    bdiv_stats        = function () bdiv_stats(hmp50, stat.by = "Sex", cpus = 1),
    bdiv_table        = function () bdiv_table(hmp50, cpus = 1),
    distmat_ord_table = function () distmat_ord_table(dm, cpus = 1),
    distmat_stats     = function () distmat_stats(dm, groups = pull(hmp50, "Sex"), cpus = 1),
    mtx_rarefy        = function () mtx_rarefy(counts, cpus = 1),
    rare_corrplot     = function () rare_corrplot(hmp50, cpus = 1),
    rare_multiplot    = function () rare_multiplot(hmp50, cpus = 1),
    rarefy            = function () rarefy(hmp50, cpus = 1),
    rarefy_cols       = function () rarefy_cols(counts, cpus = 1),
    stats_boxplot     = function () stats_boxplot(df, x = "Sex", cpus = 1),
    stats_corrplot    = function () stats_corrplot(df, x = "Age", cpus = 1),
    stats_table       = function () stats_table(df, stat.by = "Sex", cpus = 1),
    taxa_boxplot      = function () taxa_boxplot(hmp50, cpus = 1),
    taxa_corrplot     = function () taxa_corrplot(hmp50, x = "Age", cpus = 1),
    taxa_da           = function () taxa_da(hmp50, stat.by = "Sex", cpus = 1),
    taxa_forest       = function () taxa_forest(hmp50, stat.by = "Sex", cpus = 1),
    taxa_stats        = function () taxa_stats(hmp50, stat.by = "Sex", cpus = 1),
    write_xlsx        = function () write_xlsx(hmp50, tempfile(fileext = ".xlsx"), cpus = 1) )

  # A new function with `cpus` needs a line here.
  fns <- getNamespaceExports('rbiom')
  fns <- fns[vapply(fns, function (f) {
    f <- get(f, envir = asNamespace('rbiom'))
    is.function(f) && 'cpus' %in% names(formals(f)) }, logical(1))]
  expect_setequal(names(calls), fns)

  for (fn in names(calls)) {
    e <- tryCatch(calls[[fn]](), rbiom_test_blas_limit = identity)
    expect_s3_class(e, "rbiom_test_blas_limit")
    expect_equal(e$cpus, 1, label = fn)
    expect_identical(e$fun, get(fn, envir = asNamespace('rbiom')), label = fn)
  }
})


test_that("cpus limits a multi-threaded BLAS during an ordination", {

  skip_if_not(blas_ctl(), "ecodive can't control the BLAS")
  skip_if(is.na(blas_threads()) || blas_threads() < 2, "BLAS is single-threaded")
  rlang::local_options(rbiom.cache_dir = FALSE)

  before <- blas_threads()
  seen   <- NULL
  orig   <- distmat_ord_table

  local_mocked_bindings(distmat_ord_table = function (...) {
    seen <<- c(seen, blas_threads())
    orig(...) })

  bdiv_ord_table(hmp50, cpus = 1)
  expect_identical(seen, 1L)
  expect_identical(blas_threads(), before)

  # The default, all available cores, leaves the library's count as it
  # is unless fewer cores are available.
  seen <- NULL
  bdiv_ord_table(hmp50)
  expect_identical(seen, if (n_cpus() < before) 1L else before)


  # Restored after an error, too.
  local_mocked_bindings(distmat_ord_table = function (...) {
    seen <<- blas_threads()
    stop("boom") })

  seen <- NULL
  expect_error(bdiv_ord_table(hmp50, cpus = 1), "boom")
  expect_identical(seen, 1L)
  expect_identical(blas_threads(), before)
})
