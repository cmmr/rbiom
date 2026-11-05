
# Replays `p$code` in a clean environment and compares what ggplot2 draws.
expect_replays <- function (p) {
  code <- p$code
  expect_s3_class(code, 'rbiom_code')
  env   <- new.env(parent = globalenv())
  r     <- eval(parse(text = code), envir = env)
  built <- function (x) lapply(ggplot2::ggplot_build(x)$data, function (d) { d <- as.data.frame(d); rownames(d) <- NULL; d })
  expect_equal(built(r), built(p), tolerance = 1e-6)
  invisible(r)
}

# The ggplot2 function names used by a plot's layers, in order.
layer_fns <- function (p) {
  regmatches(p$code, gregexpr("[A-Za-z_:.]+\\(", p$code))[[1]] %>% sub("\\($", "", .)
}


test_that("boxplot_build layers", {

  skip_on_cran()

  params <- list2env(list(df = structure(list(metadata = NA), class = 'rbiom')))
  expect_error(boxplot_build(params))

  df <- adiv_table(rare50)
  expect_warning(stats_boxplot(df, x = 'Sex', facet.by = 'Sex'))
  expect_warning(stats_boxplot(df, stat.by = 'Sex', facet.by = 'Sex'))

  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'vbx', outliers = FALSE, test = 'none'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'crossbar',   ci = 'range', test = 'none'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'errorbar',   ci = 'mad',   test = 'none'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'linerange',  ci = 'sd',    test = 'none'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'pointrange', ci = 'se',    test = 'none'))
  expect_silent(stats_boxplot(df, x = 'Sex', facet.by = 'Body Site', layers = 'pointrange', scales = 'free_y', test = 'none'))
  expect_silent(stats_boxplot(df, x = 'Sex', facet.by = 'Body Site', layers = 'x', scales = 'fixed', test = 'none'))
  expect_error(stats_boxplot(df, x = 'Sex', layers = 'x', scales = 'free', test = 'none'), "scales")
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'cv', test = 'none'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', flip = TRUE, test = 'none'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'b', yaxis.transform = 'sqrt', test = 'none'))
  expect_error(stats_boxplot(df, x = 'Sex', layers = 'b', yaxis.transform = 'log10', test = 'none'), "sqrt")

  # Every CI type replays (intervals are recomputed from `data` by the code).
  for (ci in c('ci', 'range', 'sd', 'se', 'mad'))
    expect_replays(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'xe', ci = ci, test = 'none'))
})


test_that("boxplot_build intervals are not clamped", {

  df <- data.frame(g = rep(c("a", "b"), each = 4), v = c(-3, -2, -1, -4, 1, 2, 3, 0.5))
  p  <- stats_boxplot(df, x = "g", y = "v", layers = "e", test = "none")
  iv <- attr(p$data, 'vline')
  expect_true(all(iv$.ymin[iv$g == "a"] < 0))
  expect_true(all(iv$v[iv$g == "a"] < 0))
  
  # A plain data.frame has no provenance; the code says so.
  expect_true(grepl("# data <- <data.frame; 8 x 2>   # not reproducible", p$code, fixed = TRUE))
  env <- new.env(parent = globalenv())
  env$data <- df
  r <- eval(parse(text = p$code), envir = env)
  expect_equal(ggplot2::layer_data(r, 1), ggplot2::layer_data(p, 1))
})


test_that("boxplot_build x = NULL", {

  skip_on_cran()

  df <- adiv_table(rare50)

  p <- stats_boxplot(df, flip = TRUE, test = 'none')
  expect_true(grepl('aes(x = "", y = .diversity)', p$code, fixed = TRUE))
  expect_true(grepl("coord_flip", p$code, fixed = TRUE))
  expect_replays(p)

  # A single group of samples compared across stat.by, labeled at x = 1.
  p <- stats_boxplot(df, stat.by = 'Sex', p.label = 1)
  expect_equal(attr(p$data, 'stat_labels')$.x, 1)
})


test_that("boxplot_build xlab.angle", {

  df <- adiv_table(rare50)

  p <- stats_boxplot(df, x = 'Body Site', test = 'none', xlab.angle = 90)
  expect_true(grepl("angle = -90", p$code, fixed = TRUE))
  p <- stats_boxplot(df, x = 'Body Site', test = 'none', xlab.angle = 0)
  expect_false(grepl("angle", p$code, fixed = TRUE))
  p <- stats_boxplot(df, x = 'Body Site', test = 'none', xlab.angle = 'auto')
  expect_true(grepl("angle = -30", p$code, fixed = TRUE))
  p <- stats_boxplot(df, x = 'Sex', test = 'none', xlab.angle = 'auto')
  expect_false(grepl("angle", p$code, fixed = TRUE))
})


test_that("boxplot_build dot and strip layers", {

  skip_on_cran()

  df <- adiv_table(rare50)

  # Without ggbeeswarm, jittered points are used instead.
  old <- options(rbiom.disable_pkgs = "ggbeeswarm")
  on.exit(options(old), add = TRUE)

  p <- stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'ds', test = 'none')
  expect_false(grepl("ggbeeswarm", p$code, fixed = TRUE))
  expect_true(grepl("position_jitterdodge", p$code, fixed = TRUE))
  expect_replays(p)

  p <- stats_boxplot(df, x = 'Sex', layers = 'd', test = 'none')
  expect_true(grepl("position_jitter(", p$code, fixed = TRUE))
  expect_replays(p)

  options(rbiom.disable_pkgs = NULL)
  skip_if_not_installed("ggbeeswarm")

  p <- expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'xds', test = 'none'))
  expect_true(grepl("ggbeeswarm::geom_beeswarm", p$code, fixed = TRUE))
  expect_true(grepl("ggbeeswarm::geom_quasirandom", p$code, fixed = TRUE))
  expect_true(grepl("round(data[[\".diversity\"]], 12)", p$code, fixed = TRUE))
  expect_replays(p)
})


test_that("boxplot_build patterns", {

  skip_on_cran()

  df <- adiv_table(rare50)

  old <- options(rbiom.disable_pkgs = "fillpattern")
  on.exit(options(old), add = TRUE)
  expect_error(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', patterns = TRUE, test = 'none'), "fillpattern")

  options(rbiom.disable_pkgs = NULL)
  skip_if_not_installed("fillpattern")

  p <- expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', patterns = TRUE, layers = 'dbsl', test = 'none'))
  expect_true(grepl("fillpattern::scale_fill_pattern", p$code, fixed = TRUE))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', patterns = TRUE, layers = 'cv', test = 'none'))
})


test_that("boxplot_build intervals fit the y-axis transform", {

  skip_on_cran()

  # CI lower bounds below zero can't be drawn on a sqrt axis: they're drawn
  # at zero, and the replayed code does the same.
  p <- taxa_boxplot(rare50, layers = "p", x = "Body Site", p.label = 0, taxa = 3)
  iv <- attr(p$data, 'vline')
  expect_true(all(iv$.ymin >= 0, na.rm = TRUE))
  expect_true(grepl("pmax(v, 0)", p$code, fixed = TRUE))
  expect_no_warning(ggplot2::ggplot_build(p))
  expect_replays(p)

  # On a log axis, bounds at or below zero go to the smallest positive y.
  df <- data.frame(g = rep(c("a", "b"), each = 4), v = c(0.01, 0.02, 5, 0.03, 1, 2, 3, 0.5))
  p  <- stats_boxplot(df, x = "g", y = "v", layers = "e", test = "none", yaxis.transform = "log10")
  expect_equal(min(attr(p$data, 'vline')$.ymin), 0.01)
  expect_no_warning(ggplot2::ggplot_build(p))

  # An untransformed axis shows negative bounds as computed.
  p <- stats_boxplot(df, x = "g", y = "v", layers = "e", test = "none")
  expect_lt(min(attr(p$data, 'vline')$.ymin), 0)
})
