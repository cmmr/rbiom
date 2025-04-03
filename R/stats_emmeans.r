
#' Test for x/y numeric relationships using regression models.
#' 
#' Calculates estimated marginal means (aka least-squares means) using the 
#' [`emmeans`][emmeans::emmeans] package. The `stats_emmeans()` function will 
#' estimate means at any location along the x-axis, and can test for 
#' differences between categories. The `stats_emtrends()` function examines 
#' the slope of the trendlines, also with support for testing differences 
#' between categories.
#' 
#' For further reading on this topic, see 
#' <https://stats.oarc.ucla.edu/r/seminars/interactions-r/>.
#' 
#' 
#' @inherit documentation_default
#' 
#' @noRd
#' 
#' @param df   A data.frame with any columns named by `regr`, `resp`, 
#'        `stat.by`, and `split.by.` Required.
#' 
#' @param regr   The predictive (independent; x-axis) variable in `df`.
#'        Must be a numeric column. Required.
#' 
#' @param resp   The response (dependent; y-axis) variable in `df`, such as 
#'        taxa abundance or alpha diversity. Must be numeric. 
#'        Default: `attr(df, 'response')`
#' 
#' @param fit   How to fit the trendline. Options are `'lm'`, `'log'`, and 
#'        `'gam'`. Default: `'lm'`
#' 
#' @param stat.by   The variable in `df` defining the statistical groups. 
#'        Must be categorical. Default: `NULL`
#'        
#' @param split.by   The variable(s) in `df` for splitting the data by prior to 
#'        any calculations. Must all be categorical. Default: `NULL`
#'        
#' @param alt   Alternative hypothesis direction. Options are `'!='` 
#'        (two-sided; not equal to `mu`), `'<'` (less than `mu`), or `'>'` 
#'        (greater than `mu`). Default: `'!='`
#'        
#' @param mu   Reference value to test against. Default: `0`
#'        
#' @param at   Calculate means or slopes at this `regr` value. Default: `NULL`, 
#'        which sets `at` to the median of `regr`'s values.
#' 
#' 
#' @return A tibble data.frame with fields from the following table:
#' 
#' | **Field**    | **Description**                                      |
#' | ------------ | ---------------------------------------------------- |
#' | .mean        | Estimated marginal mean. See [emmeans::emmeans()].   |
#' | .mean.diff   | Difference in means.                                 |
#' | .slope       | Trendline slope. See [emmeans::emtrends()].          |
#' | .slope.diff  | Difference in slopes.                                |
#' | .h1          | Alternate hypothesis.                                |
#' | .p.val       | Probability that null hypothesis is correct.         |
#' | .adj.p       | `.p.val` after adjusting for multiple comparisons.   |
#' | .effect.size | Effect size. See [emmeans::eff_size()].              |
#' | .lower       | Confidence interval lower bound.                     |
#' | .upper       | Confidence interval upper bound.                     |
#' | .se          | Standard error.                                      |
#' | .n           | Number of samples.                                   |
#' | .df          | Degrees of freedom.                                  |
#' | .stat        | Wilcoxon or Kruskal-Wallis rank sum statistic.       |
#' | .t.ratio     | `.mean` / `.se`                                      |
#' | .r.sqr       | Percent of variation explained by the model.         |
#' | .adj.r       | `.r.sqr`, taking degrees of freedom into account.    |
#' | .aic         | Akaike Information Criterion (predictive models).    |
#' | .bic         | Bayesian Information Criterion (descriptive models). |
#' | .loglik      | Log-likelihood goodness-of-fit score.                |
#' | .fit.p       | P-value for observing this fit by chance.            |
#' 
#' 
#' 
#' 
#' @examples
#'     library(rbiom)
#'     library(ggplot2)
#'     
#'     df <- adiv_table(babies)
#'     
#'     # Does Shannon diversity increase with Age? ------------------------
#'     stats_emtrends(df, regr = "Age") %>% select(!.n:.df)
#'     
#'     
#'     # Does Sex make a difference? ------------------------
#'     stats_emtrends(df, "Age", split.by = "Sex") %>% select(!.n:.upper)
#'     
#'     stats_emtrends(df, "Age", stat.by = "Sex") %>% select(!.n:.t.ratio)
#'     
#'     
#'     # Run a generalized additive model at three time points. ----------------------
#'     stats_emmeans(df, "Age", fit = "gam", at = c(30, 90, 150))
#'     
#'     stats_emtrends(df, "Age", fit = "gam", at = c(30, 90, 150))
#'     
#'     
#'     # Pipe results into a figure. ------------------------
#'     stats_emmeans(df, "Age", fit = "gam", at = c(1:40) * 5, split.by = "Sex") %>% 
#'       ggplot(aes(x = `Age (days)`)) + 
#'       geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = Sex), alpha = 0.4) + 
#'       geom_line(aes(y = .mean, color = Sex))
#'     
#'     
#'     # See the underlying stats calculations: ------------------------
#'     stats_emtrends(df, regr = "Age")$code
#'     

stats_emmeans <- function (
    df, regr = NULL, resp = attr(df, 'response'), stat.by = NULL, split.by = NULL, 
    fit = "gam", at = NULL, level = 0.95, alt = "!=", mu = 0, p.adj = 'fdr' ) {
  
  stats_validate()
  
  # for CRAN check only
  emmean <- contrast <- estimate <- NULL
  p.value <- lower.CL <- upper.CL <- SE <- t.ratio <- NULL
  .regr <- .p.val <- .t.ratio <- .se <- .stat.by <- NULL
  
  
  func     <- NULL
  pairwise <- NULL
  .pairs   <- P('graphics::pairs')
  
  if (is.null(regr) && is.null(stat.by)) {
    
    pairwise <- FALSE # CRAN check requires consistent func signature
    func <- function (data, pair) {
      
      data %>% 
        stats::lm(.resp ~ 1, .) %>% 
        emmeans::emmeans(
          specs = '1', 
          infer = TRUE, 
          level = level ) %>% 
        summary(side = alt, null = mu) %>% 
        with(tibble(
          .mean    = emmean,
          .h1      = factor(paste(alt, mu)),
          .p.val   = p.value,
          .lower   = lower.CL,
          .upper   = upper.CL,
          .se      = SE,
          .n       = nrow(data),
          .df      = as.integer(df),
          .t.ratio = t.ratio ))
    }
    
    code <- glue("
      stats <- data %>% 
        stats::lm(.resp ~ 1, .) %>% 
        emmeans::emmeans(
          specs = '1', 
          infer = TRUE, 
          level = {level} ) %>% 
        summary(side = '{alt}', null = {mu}) %>% 
        with(tibble(
          .mean    = emmean,
          .h1      = factor('{paste(alt, mu)}'),
          .p.val   = p.value,
          .lower   = lower.CL,
          .upper   = upper.CL,
          .se      = SE,
          .n       = nrow(data),
          .df      = as.integer(df),
          .t.ratio = t.ratio ))" )
    
  }
  
  else if (!is.null(regr) && is.null(stat.by)) {
    
    pairwise <- FALSE # CRAN check requires consistent func signature
    func <- function (data, pair) {
      
      model <- stats_fit_model(data, fit, stat.by = FALSE, regr = TRUE)
      
      gof <- stats_glance(model)
      
      # Find pseudo-minimum p-value
      if (p_slice <- is.null(at)) {
        lo <- min(data$.regr)
        hi <- max(data$.regr)
        at <- labeling::extended(dmin = lo, dmax = hi, m = 100)
        at <- at[at > lo & at < hi]
        remove("lo", "hi")
      }
      
      stats <- model %>% 
        emmeans::emmeans(
          specs = ".regr", 
          at    = list(.regr = at), 
          infer = TRUE, 
          level = level ) %>% 
        summary(side = alt, null = mu) %>% 
        with(tibble(
          !!regr  := .regr, 
          .mean    = emmean,
          .h1      = factor(paste(alt, mu)),
          .p.val   = p.value,
          .lower   = lower.CL,
          .upper   = upper.CL,
          .se      = SE,
          .n       = nrow(data),
          .df      = as.integer(df),
          .t.ratio = t.ratio,
          .r.sqr   = gof$r.squared,
          .adj.r   = gof$adj.r.squared,
          .aic     = gof$AIC,
          .bic     = gof$BIC,
          .loglik  = gof$logLik,
          .fit.p   = gof$p.value ))
      
      if (p_slice)
        stats <- stats %>%
          dplyr::arrange(.p.val, dplyr::desc(.t.ratio), .se) %>% 
          dplyr::slice_head(n = 1)
      
      return (stats)
    }
    
    
    if (is.null(at)) {
      at_str <- "scales::breaks_extended(100)(data$.regr)"
    } else {
      at_str <- as.args(list(at))
    }
    
    
    code <- glue("
      {stats_model_code('data', fit, stat.by = FALSE)}
      
      gof <- broom::glance(model)
      
      stats <- model %>% 
        emmeans::emmeans(
          specs = '.regr', 
          at    = list(.regr = {at_str}), 
          infer = TRUE, 
          level = {level} ) %>% 
        summary(side = '{alt}', null = {mu}) %>% 
        with(tibble(
          {single_quote(regr)} = .regr,
          .mean    = emmean,
          .h1      = factor('{paste(alt, mu)}'),
          .p.val   = p.value,
          .lower   = lower.CL,
          .upper   = upper.CL,
          .se      = SE,
          .n       = nrow(data),
          .df      = as.integer(df),
          .t.ratio = t.ratio,
          .r.sqr   = gof$r.squared,
          .adj.r   = gof$adj.r.squared,
          .aic     = gof$AIC,
          .bic     = gof$BIC,
          .loglik  = gof$logLik,
          .fit.p   = gof$p.value ))" )
    
  }
  
  else if (is.null(regr) && !is.null(stat.by)) {
    
    pairwise <- TRUE # CRAN check requires consistent func signature
    func <- function (data, pair) {
      
      model <- data %>%
        subset(.stat.by %in% pair) %>%
        stats::lm(.resp ~ .stat.by, .)
      
      emm <- model %>% 
        emmeans::emmeans(
          specs = '.stat.by', 
          infer = TRUE )
      
      eff <- emm %>% 
        emmeans::eff_size(
          sigma  = stats::sigma(model), 
          edf    = stats::df.residual(model), 
          simple = '.stat.by' ) %>%
        as.data.frame()
      
      emm %>% 
        .pairs(simple = '.stat.by') %>%
        as.data.frame() %>% 
        with(tibble(
          !!stat.by   := contrast, 
          .mean.diff   = estimate, 
          .h1          = factor('!= 0'), 
          .p.val       = p.value, 
          .effect.size = eff$effect.size, 
          .se          = SE, 
          .n           = sum(data$.stat.by %in% pair),
          .df          = as.integer(df), 
          .t.ratio     = t.ratio ))
    }
    
    code <- glue("
      model <- data %>%
        subset(.stat.by %in% pair) %>%
        stats::lm(.resp ~ .stat.by, .)
      
      emm <- model %>% 
        emmeans::emmeans(
          specs = '.stat.by', 
          infer = TRUE )
      
      eff <- emm %>% 
        emmeans::eff_size(
          sigma  = stats::sigma(model), 
          edf    = stats::df.residual(model), 
          simple = '.stat.by' ) %>%
        as.data.frame()
      
      stats <- emm %>% 
        pairs(simple = '.stat.by') %>%
        as.data.frame() %>% 
        with(tibble(
          {single_quote(stat.by)} = contrast, 
          .mean.diff   = estimate, 
          .h1          = factor('!= 0'), 
          .p.val       = p.value, 
          .effect.size = eff$effect.size, 
          .se          = SE, 
          .n           = sum(data$.stat.by %in% pair), 
          .df          = as.integer(df), 
          .t.ratio     = t.ratio ))")
    
  }
  
  else if (!is.null(regr) && !is.null(stat.by)) {
    
    pairwise <- TRUE # CRAN check requires consistent func signature
    func <- function (data, pair) {
      
      pdata <- subset(data, .stat.by %in% pair)
      model <- stats_fit_model(pdata, fit, TRUE)
      
      gof <- stats_glance(model)
      
      # Find pseudo-minimum p-value
      if (p_slice <- is.null(at)) {
        lo <- max(stats::aggregate(pdata, .regr ~ .stat.by, min)$.regr)
        hi <- min(stats::aggregate(pdata, .regr ~ .stat.by, max)$.regr)
        at <- labeling::extended(dmin = lo, dmax = hi, m = 100)
        at <- at[at > lo & at < hi]
        remove("lo", "hi")
      }
      
      emm <- model %>% 
        emmeans::emmeans(
          infer  = TRUE, 
          specs  = c('.regr', '.stat.by'), 
          at     = list(.regr = at) )
      
      eff <- emm %>% 
        emmeans::eff_size(
          sigma  = stats::sigma(model), 
          edf    = stats::df.residual(model), 
          simple = '.stat.by' ) %>%
        as.data.frame()
      
      stats <- emm %>% 
        .pairs(simple = '.stat.by') %>%
        as.data.frame() %>% 
        with(tibble(
          !!regr      := .regr, 
          !!stat.by   := contrast, 
          .mean.diff   = estimate, 
          .h1          = factor('!= 0'), 
          .p.val       = p.value,
          .effect.size = eff$effect.size, 
          .se          = SE,  
          .n           = nrow(pdata),
          .df          = as.integer(df),
          .t.ratio     = t.ratio, 
          .r.sqr       = gof$r.squared,
          .adj.r       = gof$adj.r.squared,
          .aic         = gof$AIC,
          .bic         = gof$BIC,
          .loglik      = gof$logLik,
          .fit.p       = gof$p.value ))
      
      if (p_slice)
        stats <- stats %>%
          dplyr::arrange(.p.val, dplyr::desc(.t.ratio), .se) %>% 
          dplyr::slice_head(n = 1)
      
      return (stats)
    }
    
    
    if (is.null(at)) {
      at_str <- "at"
      at_def <- glue("
        
        at    <- labeling::extended(
          dmin = max(stats::aggregate(pdata, .regr ~ .stat.by, min)$.regr), 
          dmax = min(stats::aggregate(pdata, .regr ~ .stat.by, max)$.regr), 
          m    = 100 ){ifelse(fit == 'gam', '', '\n        ')}
        ")
    } else {
      at_str <- as.args(list(at))
      at_def <- ""
    }
    
    
    code <- glue("
      pdata <- subset(data, .stat.by %in% pair){at_def}
      {stats_model_code('pdata', fit, stat.by = TRUE)}
      
      gof <- broom::glance(model)
      
      emm <- model %>% 
        emmeans::emmeans(
        specs  = c('.regr', '.stat.by'), 
        at     = list(.regr = {at_str}), 
        infer  = TRUE )
      
      eff <- emm %>% 
        emmeans::eff_size(
          sigma  = stats::sigma(model), 
          edf    = stats::df.residual(model), 
          simple = '.stat.by' ) %>%
        as.data.frame()
      
      stats <- emm %>% 
        pairs(simple = '.stat.by') %>%
        as.data.frame() %>% 
        with(tibble(
          {single_quote(regr)} = .regr, 
          {single_quote(stat.by)} = contrast, 
          .mean.diff   = estimate, 
          .h1          = factor('!= 0'), 
          .p.val       = p.value,
          .effect.size = eff$effect.size, 
          .se          = SE,  
          .n           = nrow(pdata), 
          .df          = as.integer(df),
          .t.ratio     = t.ratio, 
          .r.sqr       = gof$r.squared,
          .adj.r       = gof$adj.r.squared,
          .aic         = gof$AIC,
          .bic         = gof$BIC,
          .loglik      = gof$logLik,
          .fit.p       = gof$p.value ))" )
    
  }
  
  
  
  if (!is.null(regr) && is.null(at))
    code %<>% paste0(' %>% \n  dplyr::arrange(.p.val, desc(.t.ratio), .se) %>% \n  dplyr::slice_head(n = 1)')
  
  stats <- stats_run(df, stat.by, split.by, func, pairwise) %>% aa(code = code)
  stats <- stats_finalize(stats, df, regr, resp, stat.by, split.by, fit, p.adj)
  
  return (stats)
}







stats_emtrends <- function (
    df, regr = NULL, resp = attr(df, 'response'), stat.by = NULL, split.by = NULL, 
    fit = "gam", at = NULL, level = 0.95, alt = "!=", mu = 0, p.adj = 'fdr' ) {
  
  stats_validate()
  
  if (is.null(regr))
    cli_abort('`regr` field is required when `test` = "emtrends".')
  
  
  # for CRAN check only
  .regr <- .regr.trend <- .stat.by <- NULL
  p.value <- lower.CL <- upper.CL <- NULL
  SE <- t.ratio <- contrast <- estimate <- NULL
  
  
  func     <- NULL
  pairwise <- NULL
  .pairs   <- P('graphics::pairs')
  
  if (is.null(stat.by)) {
    
    pairwise <- FALSE
    func <- function (data, pair) {
      
      model <- stats_fit_model(data, fit, stat.by = FALSE)
      
      gof <- stats_glance(model)
      
      # Find pseudo-minimum p-value
      if (p_slice <- is.null(at)) {
        lo <- min(data$.regr)
        hi <- max(data$.regr)
        at <- labeling::extended(dmin = lo, dmax = hi, m = 100)
        at <- at[at > lo & at < hi]
        remove("lo", "hi")
      }
      
      stats <- model %>% 
        emmeans::emtrends(
          specs = NULL, 
          var   = '.regr', 
          at    = list(.regr = at), 
          infer = TRUE, 
          level = level ) %>% 
        summary(side = alt, null = mu) %>%
        with(tibble(
          !!regr  := .regr, 
          .slope   = .regr.trend,
          .h1      = factor(paste(alt, mu)),
          .p.val   = p.value,
          .lower   = lower.CL,
          .upper   = upper.CL,
          .se      = SE,
          .n       = nrow(data),
          .df      = as.integer(df),
          .t.ratio = t.ratio,
          .r.sqr   = gof$r.squared,
          .adj.r   = gof$adj.r.squared,
          .aic     = gof$AIC,
          .bic     = gof$BIC,
          .loglik  = gof$logLik,
          .fit.p   = gof$p.value ))
      
      if (p_slice)
        stats <- stats[which.min(stats$.p.val),,drop=FALSE]
      
      return (stats)
      
    }
    
    
    if (is.null(at)) {
      at_str <- "scales::breaks_extended(100)(data$.regr)"
    } else {
      at_str <- as.args(list(at))
    }
    
    
    code <- glue("
      {stats_model_code('data', fit, stat.by = FALSE)}
      
      gof <- broom::glance(model)
      
      stats <- model %>% 
        emmeans::emtrends(
          specs = NULL, 
          var   = '.regr', 
          at    = list(.regr = {at_str}), 
          infer = TRUE, 
          level = {level} ) %>% 
        summary(side = '{alt}', null = {mu}) %>% 
        with(tibble(
          {single_quote(regr)} = .regr, 
          .slope   = .regr.trend,
          .h1      = factor('{paste(alt, mu)}'),
          .p.val   = p.value,
          .lower   = lower.CL,
          .upper   = upper.CL,
          .se      = SE,
          .n       = nrow(data),
          .df      = as.integer(df),
          .t.ratio = t.ratio,
          .r.sqr   = gof$r.squared,
          .adj.r   = gof$adj.r.squared,
          .aic     = gof$AIC,
          .bic     = gof$BIC,
          .loglik  = gof$logLik,
          .fit.p   = gof$p.value ))" )
    
    
  }
  
  else if (!is.null(stat.by)) {
    
    pairwise <- TRUE
    func <- function (data, pair) {
      
      pdata <- subset(data, .stat.by %in% pair)
      model <- stats_fit_model(pdata, fit, stat.by = TRUE)
      
      gof <- stats_glance(model)
      
      # Find pseudo-minimum p-value
      if (p_slice <- is.null(at)) {
        lo <- max(stats::aggregate(pdata, .regr ~ .stat.by, min)$.regr)
        hi <- min(stats::aggregate(pdata, .regr ~ .stat.by, max)$.regr)
        at <- labeling::extended(dmin = lo, dmax = hi, m = 100)
        at <- at[at > lo & at < hi]
        remove("lo", "hi")
      }
      
      
      emt <- model %>% 
        emmeans::emtrends(
        specs = '.stat.by', 
        var   = '.regr', 
        at    = list(.regr = at), 
        infer = TRUE )
      
      eff <- emt %>% 
        emmeans::eff_size(
        sigma = stats::sigma(model), 
        edf   = stats::df.residual(model) ) %>%
        as.data.frame()
      
      stats <- emt %>% 
        .pairs() %>%
        as.data.frame() %>% 
        with(tibble(
          !!regr      := at, 
          !!stat.by   := contrast, 
          .slope.diff  = estimate, 
          .h1          = factor('!= 0'), 
          .p.val       = p.value,
          .effect.size = eff$effect.size, 
          .se          = SE, 
          .n           = nrow(pdata),
          .df          = as.integer(df), 
          .t.ratio     = t.ratio, 
          .r.sqr       = gof$r.squared,
          .adj.r       = gof$adj.r.squared,
          .aic         = gof$AIC,
          .bic         = gof$BIC,
          .loglik      = gof$logLik,
          .fit.p       = gof$p.value ))
      
      if (p_slice)
        stats <- stats[which.min(stats$.p.val),,drop=FALSE]
      
      return (stats)
    }
    
    
    if (is.null(at)) {
      at_str <- "at"
      at_def <- glue("
        
        at    <- labeling::extended(
          dmin = max(stats::aggregate(pdata, .regr ~ .stat.by, min)$.regr), 
          dmax = min(stats::aggregate(pdata, .regr ~ .stat.by, max)$.regr), 
          m    = 100 ){ifelse(fit == 'gam', '', '\n        ')}
        ")
    } else {
      at_str <- as.args(list(at))
      at_def <- ""
    }
    
    
    code <- glue("
      pdata <- subset(data, .stat.by %in% pair){at_def}
      {stats_model_code('pdata', fit, stat.by = TRUE)}
      
      gof <- broom::glance(model)
      
      emt <- model %>% 
        emmeans::emmeans(
        specs = '.stat.by', 
        var   = '.regr', 
        at    = list(.regr = {at_str}), 
        infer = TRUE )
      
      eff <- emt %>% 
        emmeans::eff_size(
          sigma = stats::sigma(model), 
          edf   = stats::df.residual(model) ) %>%
        as.data.frame()
      
      stats <- emt %>% 
        pairs() %>%
        as.data.frame() %>% 
        with(tibble(
          {single_quote(regr)} = {at_str}, 
          {single_quote(stat.by)} = contrast, 
          .slope.diff  = estimate, 
          .h1          = factor('!= 0'), 
          .p.val       = p.value,
          .effect.size = eff$effect.size, 
          .se          = SE, 
          .n           = nrow(pdata),
          .df          = as.integer(df), 
          .t.ratio     = t.ratio, 
          .r.sqr       = gof$r.squared,
          .adj.r       = gof$adj.r.squared,
          .aic         = gof$AIC,
          .bic         = gof$BIC,
          .loglik      = gof$logLik,
          .fit.p       = gof$p.value ))" )
    
  }
  
  
  
  if (is.null(at))
    code %<>% paste0(' %>% \n  dplyr::slice_min(.p.val)')
  
  stats <- stats_run(df, stat.by, split.by, func, pairwise) %>% aa(code = code)
  stats <- stats_finalize(stats, df, regr, resp, stat.by, split.by, fit, p.adj)
  
  return (stats)
}



# Copied from broom package (https://github.com/tidymodels/broom)
stats_glance <- function (model) {
  
  x <- model
  s <- summary(x)
  
  if (inherits(x, 'gam')) { # broom:::glance.gam
    
    tbl <- tibble(
      df            = sum(x$edf), 
      logLik        = as.numeric(stats::logLik(x)), 
      AIC           = stats::AIC(x), 
      BIC           = stats::BIC(x), 
      deviance      = stats::deviance(x), 
      df.residual   = stats::df.residual(x), 
      nobs          = stats::nobs(x), 
      adj.r.squared = s$r.sq, 
      npar          = s$np )
    
  } else { # broom:::glance.lm
    
    tbl <- tibble(
      r.squared     = s$r.squared, 
      adj.r.squared = s$adj.r.squared, 
      sigma         = s$sigma, 
      statistic     = NA_real_, 
      p.value       = NA_real_, 
      df            = NA_real_, 
      logLik        = as.numeric(stats::logLik(x)), 
      AIC           = stats::AIC(x), 
      BIC           = stats::BIC(x), 
      deviance      = stats::deviance(x), 
      df.residual   = stats::df.residual(x), 
      nobs          = stats::nobs(x) )
    
    
    # check whether the model was fitted with only an intercept, in which
    # case don't compute the fstatistic related columns
    if (!identical(row.names(s$coefficients), "(Intercept)")) {
      
      fs <- s$fstatistic
      
      tbl$statistic <- fs[['value']]
      tbl$p.value   <- stats::pf(fs[['value']], fs[['numdf']], fs[['dendf']], lower.tail = FALSE)
      tbl$df        <- fs[['numdf']]
    }
    
  }
  
  return (tbl)
}

