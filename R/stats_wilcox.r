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

  if (!is.null(strata))
    cli_abort(c(
      'x' = "`strata` is not supported by rank tests.",
      'i' = "Use `test = \"emmeans\"` for a mixed model with a random intercept per stratum." ))

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
        test <- stats::wilcox.test(data$.resp, mu = mu, alternative = alternative, conf.int = TRUE, conf.level = level)
        tibble::tibble(
          .pseudomedian = unname(test$estimate),
          .h1           = factor(h1),
          .p.val        = test$p.value,
          .lower        = test$conf.int[[1]],
          .upper        = test$conf.int[[2]],
          .n            = nrow(data),
          .stat         = unname(test$statistic) )
      },

      code = paste(sep = "\n",
        sprintf("test <- stats::wilcox.test(data$.resp, mu = %s, alternative = %s, conf.int = TRUE, conf.level = %s)", deparse_vec(mu), double_quote(alternative), deparse_vec(level)),
        tibble_code(c(
          .pseudomedian = "unname(test$estimate)",
          .h1           = sprintf("factor(%s)", double_quote(h1)),
          .p.val        = "test$p.value",
          .lower        = "test$conf.int[[1]]",
          .upper        = "test$conf.int[[2]]",
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
        test <- stats::wilcox.test(.resp ~ .stat.by, data = pdata, conf.int = TRUE, conf.level = level)
        tibble::tibble(
          !!stat.by := paste(pair, collapse = ' - '),
          .loc.shift = unname(test$estimate),
          .h1        = factor('!= 0'),
          .p.val     = test$p.value,
          .lower     = test$conf.int[[1]],
          .upper     = test$conf.int[[2]],
          .n         = nrow(pdata),
          .stat      = unname(test$statistic) )
      },

      code = paste(sep = "\n",
        "pdata$.stat.by <- factor(pdata$.stat.by, levels = pair)",
        sprintf("test <- stats::wilcox.test(.resp ~ .stat.by, data = pdata, conf.int = TRUE, conf.level = %s)", deparse_vec(level)),
        tibble_code(c(
          stats::setNames(list("paste(pair, collapse = ' - ')"), stat.by),
          .loc.shift = "unname(test$estimate)",
          .h1        = "factor('!= 0')",
          .p.val     = "test$p.value",
          .lower     = "test$conf.int[[1]]",
          .upper     = "test$conf.int[[2]]",
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



stats_kruskal <- function (
    df, resp = attr(df, 'response'), stat.by, split.by = NULL, strata = NULL, p.adj = 'fdr') {

  if (!is.null(strata))
    cli_abort(c(
      'x' = "`strata` is not supported by rank tests.",
      'i' = "Use `test = \"emmeans\"` for a mixed model with a random intercept per stratum." ))

  if (is.null(stat.by))
    cli_abort("`stat.by` is required when `test` = 'kruskal'.")

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
