
# Replays `p$code` in a clean environment and compares what ggplot2 draws.
expect_replays <- function (p) {
  code <- p$code
  expect_s3_class(code, 'rbiom_code')
  expect_false(any(grepl("not reproducible", code, fixed = TRUE)))
  env   <- new.env(parent = globalenv())
  r     <- eval(parse(text = code), envir = env)
  built <- function (x) lapply(ggplot2::ggplot_build(x)$data, function (d) { d <- as.data.frame(d); rownames(d) <- NULL; d })
  if (inherits(p, 'patchwork')) {
    expect_true(inherits(r, 'patchwork'))
    expect_equal(length(r$patches$plots), length(p$patches$plots))
    for (i in seq_along(p$patches$plots))
      expect_equal(built(r$patches$plots[[i]]), built(p$patches$plots[[i]]), tolerance = 1e-6)
  }
  expect_equal(built(r), built(p), tolerance = 1e-6)
  invisible(r)
}


test_that("corrplot smoke", {

  df <- expect_silent(adiv_table(rare50))

  expect_silent(stats_corrplot(df,    'age', fit = 'lm'))
  expect_silent(adiv_corrplot(rare50, 'age', fit = 'lm'))
  expect_silent(bdiv_corrplot(rare50[1:20], 'age', fit = 'lm', permutations = 99))
  expect_silent(taxa_corrplot(rare50, 'age', fit = 'lm', taxa = 2))
  expect_silent(rare_corrplot(rare50[1:10]))

  skip_on_cran()

  expect_silent(adiv_corrplot(rare50, x = 'Age', stat.by = 'Sex', layers = 'ptr', colors = 'okabe', transform = 'rank', fit = 'lm'))
  expect_silent(adiv_corrplot(rare50, x = 'Age', stat.by = 'Sex', layers = 'tcn', fit = 'log'))
  expect_silent(adiv_corrplot(rare50, x = 'Age', facet.by = 'Sex', layers = 't', test = 'emtrends', fit = 'lm'))
  expect_silent(adiv_corrplot(rare50, x = 'Age', adiv = c('shannon', 'simpson'), layers = 'p', test = 'none'))
  expect_silent(taxa_corrplot(rare50, x = 'Age', taxa = 2, layers = 'ptc', fit = 'lm', p.top = 1))
  expect_silent(taxa_corrplot(rare50, x = 'Age', taxa = 1, fit = 'lm', test = 'none'))
  expect_silent(rare_corrplot(rare50[1:10], stat.by = 'Sex', adiv = c('sh', 'o'), rline = 100))
  expect_silent(rare_corrplot(rare50[1:10], rline = FALSE, layers = 'p'))
})


test_that("corrplot code replays", {

  skip_on_cran()

  expect_replays(adiv_corrplot(rare50, x = "Age", stat.by = "Sex", layers = "tcp"))
  expect_replays(adiv_corrplot(rare50, x = "Age", fit = "lm", layers = "tc"))
  expect_replays(adiv_corrplot(rare50, x = "Age", fit = "log", layers = "tr", test = "none"))
  expect_replays(adiv_corrplot(rare50, x = "Age", stat.by = "Sex", test = "emtrends", layers = "ptr"))
  expect_replays(bdiv_corrplot(rare50[1:20], x = "Age", stat.by = "Sex", fit = "lm", permutations = 99))
  expect_replays(taxa_corrplot(rare50, x = "Age", taxa = 2, layers = "ptc", fit = "lm"))
  expect_replays(taxa_corrplot(rare50, x = "Age", taxa = 6, p.top = 2, layers = "tc", fit = "lm"))
  expect_replays(rare_corrplot(rare50[1:10], stat.by = "Sex"))
  expect_replays(stats_corrplot(adiv_table(rare50), "Age", stat.by = "Sex", facet.by = "Body Site", layers = "t", fit = "lm", test = "none"))
})


test_that("corrplot model diagnostics", {

  skip_on_cran()

  p <- adiv_corrplot(rare50, x = 'Age', stat.by = 'Sex', check = TRUE, fit = 'lm')
  expect_s3_class(p, 'patchwork')
  expect_s3_class(p, 'rbiom_plot')
  expect_true(grepl("Normal Q-Q", p$code, fixed = TRUE))
  expect_true(grepl("data <- adiv_table(biom", p$code, fixed = TRUE))
  expect_true(grepl("residuals <- data", p$code, fixed = TRUE))
  expect_s3_class(p$stats, 'data.frame')
  expect_true(is.data.frame(p$data))
  expect_equal(nrow(p$data), nrow(adiv_table(rare50)))
  expect_replays(p)
})


test_that("corrplot statistics", {

  skip_on_cran()

  p <- adiv_corrplot(rare50, x = "Age", stat.by = "Sex")
  expect_s3_class(p$stats, 'data.frame')
  expect_true(hasName(p$stats, "Age"))
  expect_true(grepl("stat_vline <-", p$code, fixed = TRUE))
  expect_true(grepl("labels <-", p$code, fixed = TRUE))

  # Global tests have a label but no vertical line.
  p <- adiv_corrplot(rare50, x = "Age", fit = "lm")
  expect_false(grepl("stat_vline", p$code, fixed = TRUE))
  expect_true(grepl("labels <-", p$code, fixed = TRUE))

  # Beta diversity: Mantel test.
  p <- bdiv_corrplot(rare50[1:20], x = "Age", fit = "lm", permutations = 99)
  expect_true(grepl("Mantel", p$code, fixed = TRUE))
  expect_true(grepl("mantel", p$stats$code, fixed = TRUE))
  expect_error(bdiv_corrplot(rare50[1:20], x = "Age", test = "emmeans"), "test")

  # Strata.
  p <- adiv_corrplot(babies, x = "Age", stat.by = "Delivery", strata = "Subject", fit = "lm")
  expect_true(grepl("random intercept per Subject", p$code, fixed = TRUE))
  expect_true(grepl("Subject", p$stats$code, fixed = TRUE))

  # Empty statistics: label layer dropped, no error.
  df <- adiv_table(rare50)
  df$.diversity[df$Sex == "Male"] <- df$.diversity[df$Sex == "Male"][[1]]
  expect_null(stats_corrplot(df, "Age", test = "none")$stats)

  # A table with only NA p-values: kept as $stats, but nothing is drawn.
  df <- adiv_table(rare50)
  df <- do.call(rbind, lapply(split(df, df$`Body Site`), function (d) d[1:2, ]))
  p <- suppressWarnings(stats_corrplot(df, "Age", facet.by = "Body Site", fit = "lm", layers = "p"))
  expect_true(all(is.na(p$stats$.p.val)))
  expect_false(grepl("geom_label", p$code, fixed = TRUE))
  expect_silent(ggplot2::ggplot_build(p))

  # Residual segments are lines: no deprecated `size` aesthetic.
  expect_silent(ggplot2::ggplot_build(taxa_corrplot(rare50, "Age", taxa = 2, layers = "ptcr", fit = "lm")))

  # No statistics means no label layers at all.
  p <- adiv_corrplot(rare50, x = "Age", test = "none")
  expect_false(grepl("geom_label", p$code, fixed = TRUE))
})


test_that("corrplot provenance", {

  skip_on_cran()

  p <- adiv_corrplot(rare50, x = "Age", stat.by = "Sex", layers = "tcp", fit = "lm")
  expect_identical(as.character(p$cmd), 'adiv_corrplot(biom, "Age", layers = "tcp", stat.by = "Sex", fit = "lm")')
  expect_true(grepl("fit <- data", p$code, fixed = TRUE))
  expect_true(grepl("stats::lm(.resp ~ .regr", p$code, fixed = TRUE))
  expect_identical(as.character(p$stats$history), as.character(p$history))

  p <- rare_corrplot(rare50[1:10], stat.by = "Sex")
  expect_identical(as.character(p$cmd), 'rare_corrplot(biom, stat.by = "Sex")')
  expect_false(grepl("data <- adiv_table", p$code, fixed = TRUE))
  expect_true(grepl("depths <- ", p$code, fixed = TRUE))
  expect_true(grepl("Shannon Diversity Index", p$code, fixed = TRUE))
})


test_that("corrplot errors", {

  df <- adiv_table(rare50)

  expect_error(stats_corrplot(df, "Sex"), "not numeric")
  expect_error(stats_corrplot(df, "Age", stat.by = "Age"), "must differ")
  expect_error(stats_corrplot(df, "Age", size = 2), "not an argument")
  expect_error(adiv_corrplot(rare50, "Age", shape.by = "Sex"), "shape.by")
  expect_error(stats_corrplot(df[, setdiff(names(df), '.sample')], "Age", layers = "n", test = "none"), ".sample")
  expect_error(stats_corrplot(df, "Age", fit = "nope"))
  expect_error(stats_corrplot(1:10, "Age"), "data.frame")

  # A model that can't be fitted drops the trendline with a warning.
  tiny <- df[1:2, ]
  expect_warning(stats_corrplot(tiny, "Age", layers = "tp", test = "none"), "could not be fitted")
})


test_that("line layers use linewidth, not size", {

  skip_on_cran()

  # Soft deprecations only surface when called from the user's environment
  # (e.g. a vignette); force them here.
  old <- options(lifecycle_verbosity = "warning")
  on.exit(options(old), add = TRUE)

  p <- suppressWarnings(taxa_corrplot(rare50, layers = "tcrp", x = "Age", taxa = 1, rank = "Genus", stat.by = "Body Site", t.size = 2))
  expect_true(grepl("linewidth = 2", p$code, fixed = TRUE))
  expect_false(grepl("size = 2", p$code, fixed = TRUE))
  expect_no_warning(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p)), message = "size")
  suppressWarnings(expect_replays(p))   # mgcv notes on the tiny groups

  p <- adiv_boxplot(rare50, x = "Sex", layers = "xe", e.size = 1, box.size = 0.3, stats.size = 4)
  expect_true(grepl("linewidth += 1,", p$code))
  expect_true(grepl("linewidth += 0.3", p$code))
  expect_no_warning(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p)), message = "size")
})
