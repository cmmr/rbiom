
#' Run estimated marginal means (least-squares means).
#' 
#' All arguments should be pre-validated by [stats_table()].
#' 
#' @noRd
#' @keywords internal
#' 

stats_table_num <- function (df, stat.by, regr, resp, test, level, model, split.by) {
  
  
  #________________________________________________________
  # Sanity checks specific to corrplot-like stats.
  #________________________________________________________
  stopifnot(hasName(df, regr) && is.numeric(df[[regr]]))
  stopifnot(!eq(regr, resp))
  stopifnot(!(stat.by %in% split.by))
  
  
  #________________________________________________________
  # Sanity check model function.
  #________________________________________________________
  validate_model()
  
  
  #________________________________________________________
  # Customize the model formula with actual column names.
  #________________________________________________________
  df %<>% rename_cols(setNames(c('.regr', '.resp', '.stat.by'), c(regr, resp, stat.by)))
  model[[2]][['formula']] <- local({
    replacements <- list(x = as.symbol(".regr"), y = as.symbol(".resp"))
    f <- eval(do.call(substitute, list(model[[2]][['formula']], replacements)))
    if (!is.null(stat.by)) f <- as.formula(paste(format(f), "* .stat.by"))
    return (f)
  })
  
  
  
  #________________________________________________________
  # Iterate over `split.by` categories.
  #________________________________________________________
  stats <- plyr::ddply(df, ply_cols(split.by), function (data) {
    
    #________________________________________________________
    # Some facets might only have one level of `stat.by`.
    #________________________________________________________
    if (!is.null(stat.by))
      if (length(unique(data[['.stat.by']])) < 2)
        return (data.frame()[1,])
    
    
    #________________________________________________________
    # Apply the model to the provided data from `df`.
    #________________________________________________________
    m <- run.cmd(
      f     = model[[1]], 
      args  = c(list(data = data), model[[2]]), 
      # hist  = data, 
      # lhs   = "model", 
      envir = baseenv() )
    
    
    #________________________________________________________
    # These values can be called from multiple places below.
    #________________________________________________________
    if (is.null(stat.by)) {
      emm <- function () emmeans(object = m,  specs = '.regr',             level = level, infer = TRUE, .lhs = "emm")
      emt <- function () emtrends(object = m, specs = NULL, var = '.regr', level = level, infer = TRUE)
    } else {
      emm <- function () emmeans(object = m,  specs = c('.regr', '.stat.by'),    level = level, infer = TRUE, .lhs = "emm")
      emt <- function () emtrends(object = m, specs = '.stat.by', var = '.regr', level = level, infer = TRUE)
    }
    
    result <- tryCatch(
      error   = function (...) data.frame()[1,], 
      warning = function (...) data.frame()[1,], 
      expr    = as.data.frame(switch(
        EXPR      = test,
        predict   = m     %>% augment(data = data, se_fit = TRUE, interval = "confidence"),
        terms     = m     %>% tidy(conf.int = TRUE, conf.level = level),
        fit       = m     %>% glance(),
        means     = emm() %>% summary(adjust = 'none'),
        trends    = emt() %>% summary(adjust = 'none'),
        pw_means  = emm() %>% pairs(adjust = 'none', simple = '.stat.by'),
        pw_trends = emt() %>% pairs(adjust = 'none'),
        es_means  = emm() %>% eff_size(sigma = sigma(m), edf = df.residual(m), simple = '.stat.by'),
        es_trends = emt() %>% eff_size(sigma = sigma(m), edf = df.residual(m)),
        data.frame()[1,] )))
    
    
    #________________________________________________________
    # Drop `regr` column from stats tables.
    #________________________________________________________
    if (!test %in% c("predict", "terms", "fit"))
      result <- result[,colnames(result) != '.regr',drop=FALSE]
    
    if (test == "terms") {
      result[['term']] %<>% sub("^.stat.by", "", .)
      result[['term']] %<>% sub("^.regr:.stat.by", ".regr:", .)
      result[['term']] %<>% sub("^.regr", regr, .)
    }
    
    return (result)
    
  }) %>%
    drop_cols('.id') %>%
    rename_cols(list(
      'emmean'        = ".mean",
      'SE'            = ".se",
      'df'            = ".df",
      'lower.CL'      = ".lower",
      'upper.CL'      = ".upper",
      'conf.low'      = ".lower",
      'conf.high'     = ".upper",
      't.ratio'       = ".t.ratio",
      'p.value'       = ".p.val",
      'estimate'      = ".estimate",
      'effect.size'   = ".effect.size",
      'term'          = ".term",
      'std.error'     = ".se",
      'statistic'     = ".t.stat",
      'r.squared'     = ".r.sqr",
      'adj.r.squared' = ".adj.r.sqr",
      'sigma'         = ".sigma",
      'logLik'        = ".log.lik",
      'AIC'           = ".aic",
      'BIC'           = ".bic",
      'deviance'      = ".deviance",
      'df.residual'   = ".df.res",
      'nobs'          = ".n",
      '.regr.trend'   = ".trend",
      'contrast'      = stat.by,
      '.resp'         = resp, 
      '.regr'         = regr, 
      '.stat.by'      = stat.by )) %>%
    as_tibble() %>%
    add_class('rbiom_tbl')
  
  
  
  
  
  
  # #________________________________________________________
  # # Attach stats R commands as attr(,'code') attributes.
  # #________________________________________________________
  # 
  # stats <- local({
  #   
  #   model_fn      <- attr(model_function, 'fn', exact = TRUE)
  #   model_arg_str <- as.args(model_args)
  #   xcol_str      <- as.args(list(xcol))
  #   color.by_str  <- as.args(list(color.by))
  #   facet.by_str  <- as.args(list(facet.by))
  #   specs_str     <- as.args(list(c(xcol, color.by)))
  #   
  #   emm_template <- ifelse(
  #     test = is.null(facet.by), 
  #     yes  = paste0(
  #       sprintf("model <- %s(%s, data = data)\n", model_fn, model_arg_str),
  #       sprintf("emm   <- emmeans(object = model, specs = %s, level = %s, infer = TRUE)\n", specs_str, level),
  #       sprintf("stats <- {cmd}") ),
  #     no   = paste0(
  #       sprintf("stats <- plyr::ddply(data, %s, function (df) {\n", facet.by_str),
  #       sprintf("  model <- %s(%s, data = df)\n", model_fn, model_arg_str),
  #       sprintf("  emm   <- emmeans(object = model, specs = %s, level = %s, infer = TRUE)\n", specs_str, level),
  #       sprintf("  {cmd}\n}") ))
  #   
  #   if (hasName(stats, 'emmeans'))
  #     attr(stats[['emmeans']],  'code') <- emm_template %>%
  #     sub("{cmd}", "summary(object = emm, adjust = 'none')", ., fixed = TRUE) %>%
  #     paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
  #   
  #   if (hasName(stats, 'emm_pairs'))
  #     attr(stats[['emm_pairs']],  'code') <- emm_template %>%
  #     sub("{cmd}", sprintf("pairs(x = emm, simple = %s, adjust = 'none')", color.by_str), ., fixed = TRUE) %>%
  #     paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
  #   
  #   if (hasName(stats, 'emm_eff_size'))
  #     attr(stats[['emm_eff_size']], 'code') <- emm_template %>%
  #     sub("{cmd}", "eff_size(object = emm, sigma = sigma(model), edf = df.residual(model))", ., fixed = TRUE)
  #   
  #   
  #   
  #   emt_template <- ifelse(
  #     test = is.null(facet.by), 
  #     yes  = paste0(
  #       sprintf("model <- %s(%s, data = data)\n", model_fn, model_arg_str),
  #       sprintf("emt   <- emtrends(object = model, specs = %s, var = %s, level = %s, infer = TRUE)\n", color.by_str, xcol_str, level),
  #       sprintf("stats <- {cmd}") ),
  #     no   = paste0(
  #       sprintf("stats <- plyr::ddply(data, %s, function (df) {\n", facet.by_str),
  #       sprintf("  model <- %s(%s, data = df)\n", model_fn, model_arg_str),
  #       sprintf("  emt   <- emtrends(object = model, specs = %s, var = %s, level = %s, infer = TRUE)\n", color.by_str, xcol_str, level),
  #       sprintf("  {cmd}\n}") ))
  #   
  #   if (hasName(stats, 'emtrends'))
  #     attr(stats[['emtrends']],  'code') <- emt_template %>%
  #     sub("{cmd}", "summary(object = emt, adjust = 'none')", ., fixed = TRUE) %>%
  #     paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
  #   
  #   if (hasName(stats, 'emt_pairs'))
  #     attr(stats[['emt_pairs']],  'code') <- emt_template %>%
  #     sub("{cmd}", "pairs(x = emt, adjust = 'none')", ., fixed = TRUE) %>%
  #     paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
  #   
  #   if (hasName(stats, 'emt_eff_size'))
  #     attr(stats[['emt_eff_size']], 'code') <- emt_template %>%
  #     sub("{cmd}", "eff_size(object = emt, sigma = sigma(model), edf = df.residual(model))", ., fixed = TRUE)
  #   
  #   
  #   broom_template <- ifelse(
  #     test = is.null(facet.by), 
  #     yes  = paste0(
  #       sprintf("model <- %s(%s, data = data)\n", model_fn, model_arg_str),
  #       sprintf("stats <- {cmd}") ),
  #     no   = paste0(
  #       sprintf("stats <- plyr::ddply(data, %s, function (df) {\n", facet.by_str),
  #       sprintf("  model <- %s(%s, data = df)\n", model_fn, model_arg_str),
  #       sprintf("  {cmd}\n}") ))
  #   
  #   if (hasName(stats, 'df'))
  #     attr(stats[['df']],  'code') <- broom_template %>%
  #     sub("{cmd}", "broom::augment(x = model, data = {data})", ., fixed = TRUE) %>%
  #     sub("{data}", ifelse(is.null(facet.by), "data", "df"), ., fixed = TRUE)
  #   
  #   if (hasName(stats, 'terms'))
  #     attr(stats[['terms']],  'code') <- broom_template %>%
  #     sub("{cmd}", "broom::tidy(x = model)", ., fixed = TRUE) %>%
  #     paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
  #   
  #   if (hasName(stats, 'fit'))
  #     attr(stats[['fit']], 'code') <- broom_template %>%
  #     sub("{cmd}", "broom::glance(x = model)", ., fixed = TRUE) %>%
  #     paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
  #   
  #   
  #   return (stats)
  # })
  
  
  
  
  
  
  attr(stats, 'tbl_sum') <- c(
    
    'Test'  = switch(
      EXPR      = test,
      fit       = "Model goodness-of-fit.",
      terms     = "Significance of model terms. Pr(estimate != 0)",
      predict   = "Observed versus model-predicted values.",
      means     = "Estimated marginal means (aka least-squares means).",
      trends    = "Estimated marginal means of linear trends.",
      pw_means  = "Estimated marginal means - pairwise.",
      pw_trends = "Estimated marginal means of linear trends - pairwise.",
      es_means  = "Estimated marginal means - effect sizes.",
      es_trends = "Estimated marginal means of linear trends - effect sizes.",
      test ),
    
    'Model' = with(model, { # convert:  .resp ~ .regr * .stat.by  =>  .diversity ~ Age * `Body Site`
      
        replacements <- list(.regr = as.symbol(regr), .resp = as.symbol(resp))
        if (!is.null(stat.by)) replacements %<>% c(list(.stat.by = as.symbol(stat.by)))
        
        fn   <- attr(fun, 'fn', exact = TRUE)
        frm  <- eval(do.call(substitute, list(args[['formula']], replacements)))
        args <- c(list(frm), args[names(args) != 'formula'])
        
        sprintf("%s(%s)", fn, as.args(args))
      }) )
  
  
  if (test == "terms" && !is.null(stat.by))
    attr(stats, 'tbl_sum') %<>% c(
      'Ref. Group' = paste(coan(stat.by), "=", levels(df[['.stat.by']])[[1]]) )
  
  
  return (stats)
}
