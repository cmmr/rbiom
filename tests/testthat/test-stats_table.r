
# Evaluate `x$code` in a clean environment and return `stats`.
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

b  <- rarefy(hmp50, depth = 1000, seed = 7)
df <- adiv_table(b)


test_that("stats_table", {

  expect_silent(adiv_stats(min5, fit = 'lm'))
  expect_silent(taxa_stats(hmp5, fit = 'lm', rank = 'p', taxa = 2))
  expect_silent(adiv_stats(rare50, stat.by='Sex', regr='Age'))

  expect_error(bdiv_stats(min5), "stat.by")
  expect_error(stats_table(df, stat.by = "Sex", test = "nope"))
  expect_error(stats_table(df, stat.by = "Age"), "categorical")
  expect_error(stats_table(df, regr = "Sex"), "numeric")
  expect_error(stats_table(df, stat.by = "Sex", split.by = "Sex"), "unique")
})


test_that("kruskal / wilcox match base R", {

  x <- adiv_stats(b, stat.by = "Body Site", test = "kruskal")
  k <- stats::kruskal.test(df$.diversity ~ df$`Body Site`)
  expect_equal(x$.stat,  unname(k$statistic))
  expect_equal(x$.p.val, k$p.value)
  expect_equal(x$.df,    unname(k$parameter))
  expect_equal(x$.n,     nrow(df))

  x <- adiv_stats(b, stat.by = "Sex", test = "wilcox")
  w <- stats::wilcox.test(df$.diversity[df$Sex == "Female"], df$.diversity[df$Sex == "Male"], conf.int = TRUE)
  expect_equal(x$.p.val,     w$p.value)
  expect_equal(x$.stat,      unname(w$statistic))
  expect_equal(x$.loc.shift, unname(w$estimate))
  expect_equal(x$.lower,     w$conf.int[[1]])
  expect_equal(x$.upper,     w$conf.int[[2]])
  expect_identical(x$Sex, "Female - Male")

  # `level` is honored.
  x90 <- adiv_stats(b, stat.by = "Sex", test = "wilcox", level = 0.9)
  w90 <- stats::wilcox.test(df$.diversity[df$Sex == "Female"], df$.diversity[df$Sex == "Male"], conf.int = TRUE, conf.level = 0.9)
  expect_equal(x90$.lower, w90$conf.int[[1]])

  # Pairwise rows and the omnibus test.
  x <- adiv_stats(b, stat.by = "Body Site", test = "wilcox")
  expect_equal(nrow(x), choose(5, 2))
  expect_equal(x$.adj.p, p.adjust(x$.p.val, 'fdr'))

  expect_error(adiv_stats(b, stat.by = "Sex", test = "wilcox", strata = "Body Site"), "strata")
  expect_error(adiv_stats(b, stat.by = "Sex", test = "kruskal", strata = "Body Site"), "strata")
})


test_that("regression slopes match lm()", {

  x <- adiv_stats(b, regr = "Age", fit = "lm")
  m <- stats::lm(.diversity ~ Age, data = df)
  s <- summary(m)$coefficients
  expect_equal(x$.slope,   s['Age', 'Estimate'])
  expect_equal(x$.se,      s['Age', 'Std. Error'])
  expect_equal(x$.t.ratio, s['Age', 't value'])
  expect_equal(x$.p.val,   s['Age', 'Pr(>|t|)'])
  expect_equal(x$.r.sqr,   summary(m)$r.squared)
  expect_equal(x$.aic,     AIC(m))
  expect_equal(x$.n,       nrow(df))
  expect_false(hasName(x, 'Age'))

  # Same slope via emtrends; log fit gives the coefficient of log(Age).
  expect_equal(adiv_stats(b, regr = "Age", fit = "lm", test = "emtrends")$.slope, x$.slope)
  ml <- stats::lm(.diversity ~ log(Age), data = df)
  expect_equal(adiv_stats(b, regr = "Age", fit = "log")$.slope, unname(coef(ml)[['log(Age)']]))

  # Pairwise slope differences match the interaction coefficient.
  x <- adiv_stats(b, regr = "Age", stat.by = "Sex", fit = "lm", test = "emtrends")
  mi <- stats::lm(.diversity ~ Age * Sex, data = df)
  expect_equal(x$.slope.diff, -unname(coef(mi)[['Age:SexMale']]))
  expect_equal(x$.p.val, summary(mi)$coefficients['Age:SexMale', 'Pr(>|t|)'])

  # gam: the smooth term's significance from summary.gam.
  x <- adiv_stats(b, regr = "Age")
  g <- mgcv::gam(.diversity ~ s(Age, bs = 'cs', k = 9), data = df, method = 'REML')
  expect_equal(x$.p.val, unname(summary(g)$s.table['s(Age)', 'p-value']))
  expect_equal(x$.edf,   unname(summary(g)$s.table['s(Age)', 'edf']))
})


test_that("at defaults to the median and accepts a vector", {

  x <- adiv_stats(b, regr = "Age", stat.by = "Sex")
  expect_equal(x$Age, median(df$Age))
  expect_equal(nrow(x), 1L)

  x <- adiv_stats(b, regr = "Age", stat.by = "Sex", at = c(25, 30))
  expect_setequal(x$Age, c(25, 30))
  expect_equal(x$.adj.p, p.adjust(x$.p.val, 'fdr'))

  # gam local slopes differ between positions (the old code recycled one value).
  x <- adiv_stats(b, regr = "Age", stat.by = "Sex", test = "emtrends", at = c(25, 30))
  expect_equal(nrow(x), 2L)
  expect_false(isTRUE(all.equal(x$.slope.diff[[1]], x$.slope.diff[[2]])))

  expect_message(adiv_stats(b, regr = "Age", fit = "lm", at = 30), "ignored")
})


test_that("strata fits a random intercept", {

  skip_on_cran()

  bb <- rarefy(babies)
  x  <- adiv_stats(bb, regr = "Age", stat.by = "Delivery", strata = "Subject", fit = "lm")
  expect_true(grepl("s(.strata, bs = \"re\")", attr(x, 'tbl_sum')[['Model']], fixed = TRUE))
  expect_true(grepl("mgcv::gam", x$code, fixed = TRUE))
  expect_equal(strip(replay_stats(x)), strip(x), tolerance = 1e-6)

  x <- adiv_stats(bb, stat.by = "Delivery", strata = "Subject")
  expect_equal(strip(replay_stats(x)), strip(x), tolerance = 1e-6)

  # Ignoring repeated measures gives a different answer.
  y <- adiv_stats(bb, stat.by = "Delivery")
  expect_false(isTRUE(all.equal(x$.p.val, y$.p.val)))
})


test_that("NA rows and small groups are reported", {

  d <- df
  d$Age[1:3] <- NA
  expect_message(stats_table(d, regr = "Age", fit = "lm"), "Dropped 3")

  d <- df[df$Sex == "Male", ]
  expect_error(stats_table(d, stat.by = "Sex"), "only one group")

  d <- df[c(which(df$Sex == "Male"), which(df$Sex == "Female")[1]), ]
  expect_warning(stats_table(d, stat.by = "Sex", test = "wilcox"), "fewer than two")
})


test_that("default regression test has nominal type-I error", {

  skip_on_cran()

  set.seed(1)
  reps <- 100
  p <- vapply(seq_len(reps), FUN.VALUE = 0, function (i) {
    d <- data.frame(x = runif(40), y = rnorm(40))
    attr(d, 'response') <- 'y'
    suppressMessages(stats_table(d, regr = "x"))$.p.val
  })
  rate <- mean(p < 0.05)
  expect_gte(rate, 0.01)
  expect_lte(rate, 0.12)
})


test_that("bdiv_stats runs permutation tests on samples", {

  skip_on_cran()

  dm <- bdiv_distmat(b, "bray")

  # PERMANOVA: omnibus plus pairwise rows, matching vegan::adonis2.
  x <- bdiv_stats(b, stat.by = "Body Site", permutations = 199)
  expect_equal(nrow(x), 1 + choose(5, 2))
  expect_identical(sort(unique(x$.groups)), c(2L, 5L))
  omni <- x[x$.groups == 5, ]
  grp  <- pull(b, "Body Site")[attr(dm, 'Labels')]
  set.seed(0); a <- vegan::adonis2(dm ~ grp, permutations = 199)
  expect_equal(omni$.stat,  a$F[[1]])
  expect_equal(omni$.r.sqr, a$R2[[1]])
  expect_equal(omni$.p.val, a$`Pr(>F)`[[1]])
  expect_equal(x$.adj.p, p.adjust(x$.p.val, 'fdr'))
  expect_equal(strip(replay_stats(x)), strip(x), tolerance = 1e-6)

  # Two groups: a single row.
  x <- bdiv_stats(b, stat.by = "Sex", permutations = 99)
  expect_equal(nrow(x), 1L)
  expect_identical(x$Sex, "Female vs Male")

  # Mantel test for a numeric covariate.
  x   <- bdiv_stats(b, regr = "Age", permutations = 199)
  age <- pull(b, "Age")[attr(dm, 'Labels')]
  set.seed(0); m <- vegan::mantel(dm, dist(age), method = 'spearman', permutations = 199)
  expect_equal(x$.stat,  m$statistic)
  expect_equal(x$.p.val, m$signif)
  expect_equal(strip(replay_stats(x)), strip(x), tolerance = 1e-6)

  # Distance-based regressions are no longer offered.
  expect_error(bdiv_stats(b, stat.by = "Sex", test = "wilcox"))
  expect_error(bdiv_stats(b, stat.by = "Sex", test = "emmeans"))
  expect_error(bdiv_stats(b, regr = "Age", test = "mrpp"))

  # split.by and several metrics. One stratum (Mid vagina) has a single Sex
  # level, which is reported in a summary warning.
  expect_warning(
    x <- bdiv_stats(b, stat.by = "Sex", split.by = "Body Site", bdiv = c("bray", "jaccard"), permutations = 49),
    "could not be computed" )
  expect_true(all(c('Body Site', '.bdiv', 'Sex') %in% names(x)))
  expect_equal(strip(suppressWarnings(replay_stats(x))), strip(x), tolerance = 1e-6)

  # Repeated measures: delivery mode is constant within each infant, so
  # whole infants are permuted (not within-infant blocks, which would
  # always give p = 1).
  bb <- rarefy(subset(babies, `Age (days)` <= 60))
  x  <- bdiv_stats(bb, stat.by = "Delivery mode", strata = "Subject ID")
  expect_match(attr(x, 'tbl_sum')[['Test']], "whole `Subject ID` strata permuted .*all 220 arrangements, exact")
  expect_lt(x$.p.val, 1)
  expect_equal(x$.stat, bdiv_stats(bb, stat.by = "Delivery mode", permutations = 9)$.stat)
  expect_equal(strip(replay_stats(x)), strip(x), tolerance = 1e-6)
  expect_error(
    bdiv_stats(bb, stat.by = "Delivery mode", strata = "Subject ID", test = "mrpp"),
    "constant within each" )

  # A within-infant variable is permuted within infants.
  bb <- mutate(bb, Early = `Age (days)` <= 30)
  x  <- bdiv_stats(bb, stat.by = "Early", strata = "Subject ID", permutations = 99)
  expect_match(attr(x, 'tbl_sum')[['Test']], "99 permutations within `Subject ID`")
  expect_true(grepl("strata", x$code, fixed = TRUE))
  expect_equal(strip(replay_stats(x)), strip(x), tolerance = 1e-6)
})


test_that("stats tables work with the cache disabled", {
  skip_on_cran()
  x <- adiv_stats(b, stat.by = "Sex")
  rlang::with_options(rbiom.cache_dir = "FALSE", {
    y <- adiv_stats(b, stat.by = "Sex")
  })
  expect_equal(strip(x), strip(y))
  expect_identical(x$code, y$code)
})
