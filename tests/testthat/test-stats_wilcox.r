test_that("stats_wilcox", {

  skip_on_cran()

  df <- expect_silent(adiv_table(rare50))
  expect_silent(stats_wilcox(df))
  expect_error(stats_kruskal(df, stat.by = NULL))

  # exact = NULL: R picks the exact test for small samples without ties.
  d <- df[1:12, ]
  x <- stats_wilcox(d, stat.by = "Sex")
  w <- stats::wilcox.test(d$.diversity[d$Sex == levels(droplevels(d$Sex))[[1]]], d$.diversity[d$Sex == levels(droplevels(d$Sex))[[2]]], conf.int = TRUE)
  expect_equal(x$.p.val, w$p.value)
  expect_equal(x$.loc.shift, unname(w$estimate))

  # One-sample test against mu, honoring level and alt.
  x <- stats_wilcox(df, mu = 1, alt = ">", level = 0.9)
  w <- stats::wilcox.test(df$.diversity, mu = 1, alternative = "greater", conf.int = TRUE, conf.level = 0.9)
  expect_equal(x$.p.val, w$p.value)
  expect_equal(x$.pseudomedian, unname(w$estimate))
  expect_equal(x$.lower, w$conf.int[[1]])
  expect_identical(as.character(x$.h1), "> 1")

  # Kruskal-Wallis within strata.
  x <- stats_kruskal(df, stat.by = "Body Site", split.by = "Sex")
  expect_equal(nrow(x), 2L)
  k <- stats::kruskal.test(.diversity ~ `Body Site`, data = df[df$Sex == "Female", ])
  expect_equal(x$.p.val[x$Sex == "Female"], k$p.value)

  # Failures become NA rows with one warning.
  d <- df[df$Sex == "Male" | df$`Body Site` == "Saliva", ]
  expect_warning(stats_kruskal(d, stat.by = "Sex", split.by = "Body Site"), "could not be computed")
})
