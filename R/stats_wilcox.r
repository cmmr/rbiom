#________________________________________________________
# Rank-based tests: Wilcoxon / Mann-Whitney and Kruskal-Wallis.
#
# As in stats_emmeans.r, each test is written twice, side by side: the
# closure that computes it and the code text recorded in `$code`. Keep
# them in step. See stats_table() for the user-facing documentation.
#________________________________________________________


stats_wilcox <- function (
    df, resp = attr(df, 'response'), stat.by = NULL, split.by = NULL, strata = NULL,
    level = 0.95, alt = "!=", mu = 0, p.adj = 'fdr' ) {

  test <- "wilcox"
  validate_stats_test()

  stats_validate(df, regr = NULL, resp, stat.by, split.by, strata = NULL)

  alternative <- switch(alt, '!=' = "two.sided", '>' = "greater", '<' = "less")
  pairwise    <- !is.null(stat.by)
  h1          <- paste(alt, mu)


  if (!pairwise) {

    #________________________________________________________
    # One-sample test: is the (pseudo)median different from `mu`?
    #________________________________________________________
    test <- list(

      fun = function (data, pair, at) {
        test <- wilcox_quiet(stats::wilcox.test(data$.resp, mu = mu, alternative = alternative, conf.int = TRUE, conf.level = level))
        ci   <- if (isTRUE(attr(test$conf.int, "conf.level") < level)) c(-Inf, Inf) else test$conf.int
        tibble::tibble(
          .pseudomedian = unname(test$estimate),
          .h1           = factor(h1),
          .p.val        = test$p.value,
          .lower        = ci[[1]],
          .upper        = ci[[2]],
          .n            = nrow(data),
          .stat         = unname(test$statistic) )
      },

      code = paste(sep = "\n",
        sprintf("test <- stats::wilcox.test(data$.resp, mu = %s, alternative = %s, conf.int = TRUE, conf.level = %s)", deparse_vec(mu), double_quote(alternative), deparse_vec(level)),
        wilcox_ci_code(level),
        tibble_code(c(
          .pseudomedian = "unname(test$estimate)",
          .h1           = sprintf("factor(%s)", double_quote(h1)),
          .p.val        = "test$p.value",
          .lower        = "ci[[1]]",
          .upper        = "ci[[2]]",
          .n            = "nrow(data)",
          .stat         = "unname(test$statistic)" )) ))

    model <- "wilcox.test(.resp)"

  } else {

    #________________________________________________________
    # Two-sample test for each pair of groups. The estimate is the
    # Hodges-Lehmann location shift (first group minus second).
    #________________________________________________________
    test <- list(

      fun = function (pdata, pair, at) {
        pdata$.stat.by <- factor(pdata$.stat.by, levels = pair)
        test <- wilcox_quiet(stats::wilcox.test(.resp ~ .stat.by, data = pdata, conf.int = TRUE, conf.level = level))
        ci   <- if (isTRUE(attr(test$conf.int, "conf.level") < level)) c(-Inf, Inf) else test$conf.int
        tibble::tibble(
          !!stat.by := paste(pair, collapse = ' - '),
          .loc.shift = unname(test$estimate),
          .h1        = factor('!= 0'),
          .p.val     = test$p.value,
          .lower     = ci[[1]],
          .upper     = ci[[2]],
          .n         = nrow(pdata),
          .stat      = unname(test$statistic) )
      },

      code = paste(sep = "\n",
        "pdata$.stat.by <- factor(pdata$.stat.by, levels = pair)",
        sprintf("test <- stats::wilcox.test(.resp ~ .stat.by, data = pdata, conf.int = TRUE, conf.level = %s)", deparse_vec(level)),
        wilcox_ci_code(level),
        tibble_code(c(
          stats::setNames(list("paste(pair, collapse = ' - ')"), stat.by),
          .loc.shift = "unname(test$estimate)",
          .h1        = "factor('!= 0')",
          .p.val     = "test$p.value",
          .lower     = "ci[[1]]",
          .upper     = "ci[[2]]",
          .n         = "nrow(pdata)",
          .stat      = "unname(test$statistic)" )) ))

    model <- "wilcox.test(.resp ~ .stat.by)"
  }


  stats_finalize(
    df       = df,
    pre      = NULL,
    test     = test,
    regr     = NULL,
    resp     = resp,
    stat.by  = stat.by,
    split.by = split.by,
    strata   = NULL,
    p.adj    = p.adj,
    pairwise = pairwise,
    model    = model )
}



#________________________________________________________
# stats::wilcox.test() changed in R 4.6. Earlier versions fall back to the
# normal approximation when there are ties or zeros, and when the requested
# confidence level is out of reach (very small groups) they return a
# narrower interval at whatever level is achievable - [-1, -1] at level 0
# for one observation per group. R >= 4.6 computes exact results with ties
# and returns an unbounded interval instead.
#
# rbiom reports the unbounded interval on every R version (the code below,
# which is also written into `$code`), so the warnings about it, and about
# the fallback to the normal approximation, carry no information for the
# user and are muffled. Any other warning is passed on.
#________________________________________________________
wilcox_quiet <- function (expr) {

  msgs <- c(
    "Requested conf.level not achievable",
    "requested conf.level not achievable",
    "cannot compute exact p-value with ties",
    "cannot compute exact p-value with zeroes",
    "cannot compute exact confidence interval with ties",
    "cannot compute exact confidence intervals with ties",
    "cannot compute exact confidence interval with zeroes",
    "cannot compute exact confidence intervals with zeroes" )
  msgs <- unique(c(msgs, gettext(msgs, domain = "R-stats")))

  withCallingHandlers(
    expr    = expr,
    warning = function (w) {
      if (conditionMessage(w) %in% msgs)
        invokeRestart('muffleWarning') })
}

wilcox_ci_code <- function (level) {
  sprintf(
    "ci   <- if (isTRUE(attr(test$conf.int, \"conf.level\") < %s)) c(-Inf, Inf) else test$conf.int # unbounded when %s isn't achievable",
    deparse_vec(level), deparse_vec(level) )
}



stats_kruskal <- function (
    df, resp = attr(df, 'response'), stat.by, split.by = NULL, strata = NULL, p.adj = 'fdr') {

  test <- "kruskal"
  validate_stats_test()

  stats_validate(df, regr = NULL, resp, stat.by, split.by, strata = NULL)


  test <- list(

    fun = function (data, pair, at) {
      test <- stats::kruskal.test(.resp ~ .stat.by, data = data)
      tibble::tibble(
        .stat  = unname(test$statistic),
        .df    = unname(test$parameter),
        .p.val = test$p.value,
        .n     = nrow(data) )
    },

    code = paste(sep = "\n",
      "test <- stats::kruskal.test(.resp ~ .stat.by, data = data)",
      tibble_code(c(
        .stat  = "unname(test$statistic)",
        .df    = "unname(test$parameter)",
        .p.val = "test$p.value",
        .n     = "nrow(data)" )) ))


  stats_finalize(
    df       = df,
    pre      = NULL,
    test     = test,
    regr     = NULL,
    resp     = resp,
    stat.by  = stat.by,
    split.by = split.by,
    strata   = NULL,
    p.adj    = p.adj,
    pairwise = FALSE,
    model    = "kruskal.test(.resp ~ .stat.by)" )
}
