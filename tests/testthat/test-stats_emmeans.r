
replay_stats <- function (x) {
  expect_false(any(grepl("not reproducible", x$code, fixed = TRUE)))
  env <- new.env(parent = globalenv())
  eval(parse(text = x$code), envir = env)
  env$stats
}

strip <- function (df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (a in c('cmd', 'code', 'history', 'data_cmd', 'tbl_sum')) attr(df, a) <- NULL
  rownames(df) <- NULL
  df
}

expect_replays <- function (x) {
  expect_equal(strip(replay_stats(x)), strip(x), tolerance = 1e-6)
}


test_that("stats_emmeans", {

  df <- expect_silent(adiv_table(rare50))

  expect_silent(stats_emmeans(df))
  expect_silent(stats_emtrends(df, regr = "Age"))
  expect_error(stats_emtrends(df), "regr")

  skip_on_cran()

  expect_silent(stats_emmeans(df, regr = "Age"))
  expect_silent(stats_emmeans(df, stat.by = "Sex"))
  expect_silent(stats_emmeans(df, regr = "Age", stat.by = "Sex"))
  expect_silent(stats_emmeans(df, regr = "Age", split.by = "Sex"))
  expect_silent(stats_emmeans(df, regr = "Age", stat.by = "Sex", fit = "gam", at = c(30, 90, 150)))

  expect_silent(stats_emtrends(df, regr = "Age", stat.by = "Sex"))
  expect_silent(stats_emtrends(df, regr = "Age", split.by = "Sex"))
  expect_silent(stats_emtrends(df, regr = "Age", fit = "gam", at = c(30, 90, 150)))
})


test_that("every regression case replays from $code", {

  skip_on_cran()

  df <- adiv_table(rare50)

  for (x in list(
    stats_emmeans(df),
    stats_emmeans(df, stat.by = "Sex"),
    stats_emmeans(df, stat.by = "Body Site", split.by = "Sex"),
    stats_emmeans(df, regr = "Age"),
    stats_emmeans(df, regr = "Age", fit = "lm"),
    stats_emmeans(df, regr = "Age", fit = "log"),
    stats_emmeans(df, regr = "Age", stat.by = "Sex"),
    stats_emmeans(df, regr = "Age", stat.by = "Sex", fit = "lm", at = c(25, 30)),
    stats_emtrends(df, regr = "Age"),
    stats_emtrends(df, regr = "Age", fit = "lm"),
    stats_emtrends(df, regr = "Age", stat.by = "Sex"),
    stats_emtrends(df, regr = "Age", stat.by = "Sex", fit = "log"),
    suppressWarnings(stats_emtrends(df, regr = "Age", stat.by = "Sex", split.by = "Body Site", fit = "lm")) )) {

    expect_s3_class(x, 'rbiom_tbl')
    expect_true(all(c('.p.val', '.adj.p') %in% names(x)))
    suppressWarnings(expect_replays(x))  # the split.by case warns about a single-group stratum
  }
})


test_that("the model header names the model that was fitted", {
  df <- adiv_table(rare50)
  expect_match(attr(stats_emmeans(df, stat.by = "Sex"), 'tbl_sum')[['Model']], "^lm\\(")
  expect_match(attr(stats_emmeans(df, regr = "Age"), 'tbl_sum')[['Model']], "^gam\\(")
  expect_match(attr(stats_emmeans(df, regr = "Age", fit = "lm"), 'tbl_sum')[['Model']], "^lm\\(.resp ~ .regr\\)")
  expect_match(attr(stats_emmeans(df, regr = "Age", fit = "log"), 'tbl_sum')[['Model']], "log\\(.regr\\)")
})


test_that("generated code uses no retired packages", {
  df <- adiv_table(rare50)
  for (x in list(
    stats_emmeans(df, regr = "Age", stat.by = "Sex"),
    suppressWarnings(stats_emtrends(df, regr = "Age", stat.by = "Sex", split.by = "Body Site")) )) {
    expect_false(grepl("(^|[^d])plyr::|broom::|scales::|labeling::|\\.\\(", x$code))
  }
})
