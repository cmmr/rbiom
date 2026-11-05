
#' Regression statistics: estimated marginal means and trends.
#'
#' `stats_emmeans()` and `stats_emtrends()` fit a regression model to `df`
#' and summarize it with the [`emmeans`][emmeans::emmeans] package. They are
#' the workhorses behind `stats_table(test = "emmeans")` and
#' `stats_table(test = "emtrends")`; see [stats_table()] for the user-facing
#' documentation of the arguments and the returned columns.
#'
#' Each case below is written twice, side by side: once as the closure that
#' computes the result, and once as the code text recorded in `$code`. Keep
#' them in step. The text is a record only and is never evaluated.
#'
#' @inherit documentation_default
#'
#' @noRd
#' @keywords internal
#'

stats_emmeans <- function (
    df, regr = NULL, resp = attr(df, 'response'), stat.by = NULL, split.by = NULL,
    strata = NULL, fit = "gam", at = NULL, level = 0.95, alt = "!=", mu = 0, p.adj = 'fdr' ) {

  stats_regression(
    df = df, regr = regr, resp = resp, stat.by = stat.by, split.by = split.by,
    strata = strata, test = "emmeans", fit = fit, at = at, level = level,
    alt = alt, mu = mu, p.adj = p.adj )
}



stats_emtrends <- function (
    df, regr = NULL, resp = attr(df, 'response'), stat.by = NULL, split.by = NULL,
    strata = NULL, fit = "gam", at = NULL, level = 0.95, alt = "!=", mu = 0, p.adj = 'fdr' ) {

  if (is.null(regr))
    cli_abort('`regr` is required when `test` = "emtrends".')

  stats_regression(
    df = df, regr = regr, resp = resp, stat.by = stat.by, split.by = split.by,
    strata = strata, test = "emtrends", fit = fit, at = at, level = level,
    alt = alt, mu = mu, p.adj = p.adj )
}




#________________________________________________________
# Build and run the regression for one (test, fit, regr, stat.by)
# combination. See stats_table() for the semantics of each case.
#________________________________________________________
stats_regression <- function (
    df, regr, resp, stat.by, split.by, strata, test, fit, at, level, alt, mu, p.adj ) {

  stats_validate(df, regr, resp, stat.by, split.by, strata)

  has_regr    <- !is.null(regr)
  has_stat.by <- !is.null(stat.by)
  has_strata  <- !is.null(strata)

  # Without `regr` there is no trend to fit: the model is `.resp ~ .stat.by`
  # (or `.resp ~ 1`), fitted with lm(), or with gam() only to add `strata`.
  if (!has_regr) fit <- "lm"

  type      <- stats_model_type(fit, has_strata)
  pairwise  <- has_stat.by
  var       <- if (pairwise) "pdata" else "data"
  trend_var <- if (fit == "log") "log(.regr)" else ".regr"
  trend_col <- paste0(trend_var, ".trend")
  h1        <- paste(alt, mu)


  #________________________________________________________
  # Which cases evaluate the model at `at`?
  #________________________________________________________
  uses_at <- has_regr && (
    (test == "emmeans"  && has_stat.by) ||
    (test == "emtrends" && fit == "gam") )

  if (!is.null(at) && has_regr && !uses_at)
    cli::cli_inform(c(
      'i' = "`at` is ignored: with `fit = \"{fit}\"` the slope{?s} {?is/are} the same at every `regr` value." ),
      .frequency = "once", .frequency_id = "rbiom_at_ignored_slope" )

  if (has_regr && !has_stat.by && test == "emmeans" && fit == "gam" && (!eq(alt, "!=") || !eq(mu, 0)))
    cli::cli_inform(c(
      'i' = "`alt` and `mu` are ignored for the smooth-term test (`fit = \"gam\"`, `test = \"emmeans\"`)." ),
      .frequency = "once", .frequency_id = "rbiom_alt_mu_ignored_gam" )

  pre <- if (uses_at) {
    if (is.null(at)) {
      list(fun = function (data) median(data$.regr), code = "at <- median(data$.regr)")
    } else {
      list(fun = function (data) at, code = sprintf("at <- %s", deparse_vec(at)))
    }
  }


  #________________________________________________________
  # Model fitting, shared by the cases below.
  #________________________________________________________
  fit_model <- function (data) {
    stats_fit_model(
      df      = data,
      fit     = fit,
      stat.by = has_stat.by,
      regr    = if (has_regr)   ".regr",
      strata  = if (has_strata) ".strata" )
  }

  fit_code <- paste(sep = "\n",
    stats_model_code(var, fit, stat.by = has_stat.by, regr = has_regr, strata = has_strata),
    "gof   <- summary(model)" )

  gof_vals <- function (model, gof) stats_gof_vals(type, model, gof)
  gof_code <- stats_gof_code(type)

  level_code <- deparse_vec(level)
  side_code  <- sprintf("side = %s, null = %s", double_quote(alt), deparse_vec(mu))
  h1_code    <- sprintf("factor(%s)", double_quote(h1))



  #________________________________________________________
  # Case A: no regr, no stat.by. Is the mean different from `mu`?
  #________________________________________________________
  if (!has_regr && !has_stat.by) {

    test <- list(

      fun = function (data, pair, at) {
        model <- fit_model(data)
        emm   <- emmeans::emmeans(model, specs = ~ 1, infer = TRUE, level = level)
        summ  <- summary(emm, side = alt, null = mu)
        tibble::tibble(
          .mean    = summ$emmean,
          .h1      = factor(h1),
          .p.val   = summ$p.value,
          .lower   = summ$lower.CL,
          .upper   = summ$upper.CL,
          .se      = summ$SE,
          .n       = nrow(data),
          .df      = summ$df,
          .t.ratio = summ$t.ratio )
      },

      code = paste(sep = "\n",
        stats_model_code("data", fit, stat.by = FALSE, regr = FALSE, strata = has_strata),
        sprintf("emm   <- emmeans::emmeans(model, specs = ~ 1, infer = TRUE, level = %s)", level_code),
        sprintf("summ  <- summary(emm, %s)", side_code),
        tibble_code(c(
          .mean    = "summ$emmean",
          .h1      = h1_code,
          .p.val   = "summ$p.value",
          .lower   = "summ$lower.CL",
          .upper   = "summ$upper.CL",
          .se      = "summ$SE",
          .n       = "nrow(data)",
          .df      = "summ$df",
          .t.ratio = "summ$t.ratio" )) ))
  }


  #________________________________________________________
  # Case B: regr only.
  #________________________________________________________
  else if (has_regr && !has_stat.by) {

    if (test == "emmeans" && fit == "gam") {

      # The smooth term's approximate significance (mgcv::summary.gam).
      test <- list(

        fun = function (data, pair, at) {
          model <- fit_model(data)
          gof   <- summary(model)
          sm    <- gof$s.table['s(.regr)', ]
          tibble::tibble(
            .edf      = sm[['edf']],
            .ref.df   = sm[['Ref.df']],
            .f        = sm[['F']],
            .p.val    = sm[['p-value']],
            .r.sqr    = gof$r.sq,
            .dev.expl = gof$dev.expl,
            .aic      = stats::AIC(model),
            .bic      = stats::BIC(model),
            .loglik   = as.numeric(stats::logLik(model)),
            .n        = nrow(data) )
        },

        code = paste(sep = "\n",
          fit_code,
          "sm    <- gof$s.table['s(.regr)', ]",
          tibble_code(c(
            .edf      = "sm[['edf']]",
            .ref.df   = "sm[['Ref.df']]",
            .f        = "sm[['F']]",
            .p.val    = "sm[['p-value']]",
            .r.sqr    = "gof$r.sq",
            .dev.expl = "gof$dev.expl",
            .aic      = "stats::AIC(model)",
            .bic      = "stats::BIC(model)",
            .loglik   = "as.numeric(stats::logLik(model))",
            .n        = "nrow(data)" )) ))

    } else if (fit == "gam") {

      # Local slope of the smooth at `at`.
      test <- list(

        fun = function (data, pair, at) {
          model <- fit_model(data)
          gof   <- summary(model)
          emt   <- emmeans::emtrends(model, specs = ~ .regr, var = '.regr', at = list(.regr = at), infer = TRUE, level = level)
          summ  <- summary(emt, side = alt, null = mu)
          tibble::tibble(
            !!regr  := summ$.regr,
            .slope   = summ$.regr.trend,
            .h1      = factor(h1),
            .p.val   = summ$p.value,
            .lower   = summ$lower.CL,
            .upper   = summ$upper.CL,
            .se      = summ$SE,
            .n       = nrow(data),
            .df      = summ$df,
            .t.ratio = summ$t.ratio,
            !!!gof_vals(model, gof) )
        },

        code = paste(sep = "\n",
          fit_code,
          sprintf("emt   <- emmeans::emtrends(model, specs = ~ .regr, var = '.regr', at = list(.regr = at), infer = TRUE, level = %s)", level_code),
          sprintf("summ  <- summary(emt, %s)", side_code),
          tibble_code(c(
            stats::setNames(list("summ$.regr"), regr),
            .slope   = "summ$.regr.trend",
            .h1      = h1_code,
            .p.val   = "summ$p.value",
            .lower   = "summ$lower.CL",
            .upper   = "summ$upper.CL",
            .se      = "summ$SE",
            .n       = "nrow(data)",
            .df      = "summ$df",
            .t.ratio = "summ$t.ratio",
            gof_code )) ))

    } else {

      # The (global) slope's t-test: d .resp / d .regr, or d .resp / d log(.regr).
      test <- list(

        fun = function (data, pair, at) {
          model <- fit_model(data)
          gof   <- summary(model)
          emt   <- emmeans::emtrends(model, specs = ~ 1, var = trend_var, infer = TRUE, level = level)
          summ  <- summary(emt, side = alt, null = mu)
          tibble::tibble(
            .slope   = summ[[trend_col]],
            .h1      = factor(h1),
            .p.val   = summ$p.value,
            .lower   = summ$lower.CL,
            .upper   = summ$upper.CL,
            .se      = summ$SE,
            .n       = nrow(data),
            .df      = summ$df,
            .t.ratio = summ$t.ratio,
            !!!gof_vals(model, gof) )
        },

        code = paste(sep = "\n",
          fit_code,
          sprintf("emt   <- emmeans::emtrends(model, specs = ~ 1, var = %s, infer = TRUE, level = %s)", double_quote(trend_var), level_code),
          sprintf("summ  <- summary(emt, %s)", side_code),
          tibble_code(c(
            .slope   = sprintf("summ[[%s]]", double_quote(trend_col)),
            .h1      = h1_code,
            .p.val   = "summ$p.value",
            .lower   = "summ$lower.CL",
            .upper   = "summ$upper.CL",
            .se      = "summ$SE",
            .n       = "nrow(data)",
            .df      = "summ$df",
            .t.ratio = "summ$t.ratio",
            gof_code )) ))
    }
  }


  #________________________________________________________
  # Case C: stat.by only. Pairwise differences in group means.
  #________________________________________________________
  else if (!has_regr && has_stat.by) {

    test <- list(

      fun = function (pdata, pair, at) {
        model <- fit_model(pdata)
        emm   <- emmeans::emmeans(model, specs = ~ .stat.by)
        summ  <- summary(emmeans::contrast(emm, method = 'pairwise'), infer = TRUE)
        eff   <- as.data.frame(emmeans::eff_size(emm, sigma = stats::sigma(model), edf = stats::df.residual(model)))
        tibble::tibble(
          !!stat.by   := as.character(summ$contrast),
          .mean.diff   = summ$estimate,
          .h1          = factor('!= 0'),
          .p.val       = summ$p.value,
          .lower       = summ$lower.CL,
          .upper       = summ$upper.CL,
          .effect.size = eff$effect.size,
          .se          = summ$SE,
          .n           = nrow(pdata),
          .df          = summ$df,
          .t.ratio     = summ$t.ratio )
      },

      code = paste(sep = "\n",
        stats_model_code("pdata", fit, stat.by = TRUE, regr = FALSE, strata = has_strata),
        "emm   <- emmeans::emmeans(model, specs = ~ .stat.by)",
        "summ  <- summary(emmeans::contrast(emm, method = 'pairwise'), infer = TRUE)",
        "eff   <- as.data.frame(emmeans::eff_size(emm, sigma = stats::sigma(model), edf = stats::df.residual(model)))",
        tibble_code(c(
          stats::setNames(list("as.character(summ$contrast)"), stat.by),
          .mean.diff   = "summ$estimate",
          .h1          = "factor('!= 0')",
          .p.val       = "summ$p.value",
          .lower       = "summ$lower.CL",
          .upper       = "summ$upper.CL",
          .effect.size = "eff$effect.size",
          .se          = "summ$SE",
          .n           = "nrow(pdata)",
          .df          = "summ$df",
          .t.ratio     = "summ$t.ratio" )) ))
  }


  #________________________________________________________
  # Case D: regr and stat.by.
  #________________________________________________________
  else {

    if (test == "emmeans") {

      # Pairwise differences in group means at `at`.
      test <- list(

        fun = function (pdata, pair, at) {
          model <- fit_model(pdata)
          gof   <- summary(model)
          emm   <- emmeans::emmeans(model, specs = c('.regr', '.stat.by'), at = list(.regr = at))
          summ  <- summary(emmeans::contrast(emm, method = 'pairwise', simple = '.stat.by'), infer = TRUE)
          eff   <- as.data.frame(emmeans::eff_size(emm, sigma = stats::sigma(model), edf = stats::df.residual(model), simple = '.stat.by'))
          tibble::tibble(
            !!regr      := summ$.regr,
            !!stat.by   := as.character(summ$contrast),
            .mean.diff   = summ$estimate,
            .h1          = factor('!= 0'),
            .p.val       = summ$p.value,
            .lower       = summ$lower.CL,
            .upper       = summ$upper.CL,
            .effect.size = eff$effect.size,
            .se          = summ$SE,
            .n           = nrow(pdata),
            .df          = summ$df,
            .t.ratio     = summ$t.ratio,
            !!!gof_vals(model, gof) )
        },

        code = paste(sep = "\n",
          fit_code,
          "emm   <- emmeans::emmeans(model, specs = c('.regr', '.stat.by'), at = list(.regr = at))",
          "summ  <- summary(emmeans::contrast(emm, method = 'pairwise', simple = '.stat.by'), infer = TRUE)",
          "eff   <- as.data.frame(emmeans::eff_size(emm, sigma = stats::sigma(model), edf = stats::df.residual(model), simple = '.stat.by'))",
          tibble_code(c(
            stats::setNames(list("summ$.regr", "as.character(summ$contrast)"), c(regr, stat.by)),
            .mean.diff   = "summ$estimate",
            .h1          = "factor('!= 0')",
            .p.val       = "summ$p.value",
            .lower       = "summ$lower.CL",
            .upper       = "summ$upper.CL",
            .effect.size = "eff$effect.size",
            .se          = "summ$SE",
            .n           = "nrow(pdata)",
            .df          = "summ$df",
            .t.ratio     = "summ$t.ratio",
            gof_code )) ))

    } else if (fit == "gam") {

      # Pairwise differences in local slopes at `at`.
      test <- list(

        fun = function (pdata, pair, at) {
          model <- fit_model(pdata)
          gof   <- summary(model)
          emt   <- emmeans::emtrends(model, specs = c('.regr', '.stat.by'), var = '.regr', at = list(.regr = at))
          summ  <- summary(emmeans::contrast(emt, method = 'pairwise', simple = '.stat.by'), infer = TRUE)
          tibble::tibble(
            !!regr      := summ$.regr,
            !!stat.by   := as.character(summ$contrast),
            .slope.diff  = summ$estimate,
            .h1          = factor('!= 0'),
            .p.val       = summ$p.value,
            .lower       = summ$lower.CL,
            .upper       = summ$upper.CL,
            .se          = summ$SE,
            .n           = nrow(pdata),
            .df          = summ$df,
            .t.ratio     = summ$t.ratio,
            !!!gof_vals(model, gof) )
        },

        code = paste(sep = "\n",
          fit_code,
          "emt   <- emmeans::emtrends(model, specs = c('.regr', '.stat.by'), var = '.regr', at = list(.regr = at))",
          "summ  <- summary(emmeans::contrast(emt, method = 'pairwise', simple = '.stat.by'), infer = TRUE)",
          tibble_code(c(
            stats::setNames(list("summ$.regr", "as.character(summ$contrast)"), c(regr, stat.by)),
            .slope.diff  = "summ$estimate",
            .h1          = "factor('!= 0')",
            .p.val       = "summ$p.value",
            .lower       = "summ$lower.CL",
            .upper       = "summ$upper.CL",
            .se          = "summ$SE",
            .n           = "nrow(pdata)",
            .df          = "summ$df",
            .t.ratio     = "summ$t.ratio",
            gof_code )) ))

    } else {

      # Pairwise differences in (global) slopes.
      test <- list(

        fun = function (pdata, pair, at) {
          model <- fit_model(pdata)
          gof   <- summary(model)
          emt   <- emmeans::emtrends(model, specs = ~ .stat.by, var = trend_var)
          summ  <- summary(emmeans::contrast(emt, method = 'pairwise'), infer = TRUE)
          tibble::tibble(
            !!stat.by   := as.character(summ$contrast),
            .slope.diff  = summ$estimate,
            .h1          = factor('!= 0'),
            .p.val       = summ$p.value,
            .lower       = summ$lower.CL,
            .upper       = summ$upper.CL,
            .se          = summ$SE,
            .n           = nrow(pdata),
            .df          = summ$df,
            .t.ratio     = summ$t.ratio,
            !!!gof_vals(model, gof) )
        },

        code = paste(sep = "\n",
          fit_code,
          sprintf("emt   <- emmeans::emtrends(model, specs = ~ .stat.by, var = %s)", double_quote(trend_var)),
          "summ  <- summary(emmeans::contrast(emt, method = 'pairwise'), infer = TRUE)",
          tibble_code(c(
            stats::setNames(list("as.character(summ$contrast)"), stat.by),
            .slope.diff  = "summ$estimate",
            .h1          = "factor('!= 0')",
            .p.val       = "summ$p.value",
            .lower       = "summ$lower.CL",
            .upper       = "summ$upper.CL",
            .se          = "summ$SE",
            .n           = "nrow(pdata)",
            .df          = "summ$df",
            .t.ratio     = "summ$t.ratio",
            gof_code )) ))
    }
  }


  stats_finalize(
    df       = df,
    pre      = pre,
    test     = test,
    regr     = regr,
    resp     = resp,
    stat.by  = stat.by,
    split.by = split.by,
    strata   = strata,
    p.adj    = p.adj,
    pairwise = pairwise,
    model    = stats_model_str(fit, has_stat.by, has_regr, has_strata) )
}
