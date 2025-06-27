#________________________________________________________
# Computes p-values for categorical differences.
#________________________________________________________
corrplot_stats <- function (params) {
  
  stopifnot(is_bare_environment(params))
  
  
  # Validation of user arguments handled by corrplot_build()
  
  
  #________________________________________________________
  # Omit stats from this plot.
  #________________________________________________________
  if (params$test == 'none')
    return (params)
  
  
  
  #________________________________________________________
  # If anything goes wrong, skip stats.
  #________________________________________________________
  tryCatch(
    error = function (e) { warning(e); return (params) },
    expr  = with(params, {
      .plot_attrs[['stats']] <- stats_table(
        df       = .ggdata, 
        regr     = .xcol, 
        resp     = .ycol, 
        stat.by  = stat.by, 
        split.by = facet.by, 
        test     = test, 
        fit      = fit, 
        at       = at, 
        level    = level, 
        alt      = alt, 
        mu       = mu, 
        p.adj    = p.adj )
    }))
  
  
  
  #________________________________________________________
  # Use p.top to retain only most significant taxa.
  #________________________________________________________
  apply_p.top(params)
  if (plyr::empty(params$.ggdata))
    return (invisible(params))
  
  
  
  #________________________________________________________
  # Show p-values on the plot. One per facet.
  #________________________________________________________
  
  with(params, {
    
    attr(.ggdata, 'stat_labels') <- plyr::ddply(
      .data      = .plot_attrs[['stats']],
      .variables = ply_cols(facet.by),
      .fun       = function (df) {
        
        plabel <- paste("p =", format(df$.adj.p[[1]], digits=2))
        
        if (fit == "lm" && test == "emtrends")
          return (tibble(
            !!.xcol := -Inf, 
            !!.ycol :=  Inf,
            .label   =  plabel, 
            .hjust   = -0.1,
            .vjust   =  1.2 ))
        
        
        xpos   <- df[[.xcol]][[1]]
        xlabel <- paste("x =", xpos)
        
        if (nrow(df) > 1)
          xlabel <- paste(df[[stat.by]][[1]], "\n", xlabel)
        
        tibble(
          !!.xcol := c(xpos, xpos), 
          !!.ycol := c(Inf, -Inf),
          .label   = c(plabel, xlabel), 
          .hjust   = c(0.5,  0.5),
          .vjust   = c(1, 0) )
        
      }) %>% as_tibble()
    
    
    if (fit != "lm" || test != "emtrends")
      attr(.ggdata, 'stat_vline') <- plyr::ddply(
        .data      = .plot_attrs[['stats']],
        .variables = ply_cols(facet.by),
        .fun       = function (df) {
          tibble(!!.xcol := df[[.xcol]][[1]])
        }) %>% as_tibble()
  })
  
  
  set_layer(
    params = params,
    layer  = 'stats_label',
    'mapping|label' = ".label",
    'mapping|hjust' = ".hjust",
    'mapping|vjust' = ".vjust",
    'size'          = 4 )
  
  if ("linewidth" %in% names(ggplot2::GeomLabel$default_aes)) {
    set_layer(params, 'stats_label', 'linewidth' = NA)
  } else {
    set_layer(params, 'stats_label', 'label.size' = NA)
  }
  
  
  set_layer(
    params = params,
    layer  = 'yaxis',
    'expand' = c(0.15, 0, .15, 0) )
  
  if (!plyr::empty(attr(params$.ggdata, 'stat_vline')))
    set_layer(
      params = params,
      layer  = 'stats_vline',
      'mapping|xintercept' = params$.xcol,
      'linetype'           = 'dotdash',
      'alpha'              = 0.2 )
  
  
  
  if (isTRUE(params$caption)) {
    
    .element_text <- P('ggplot2::element_text')
    
    # element_markdown
    set_layer(params, 'theme', plot.caption = .element_text(size = 9, face = "italic"))
    set_layer(params, 'labs', .overwrite = TRUE, caption = local({
      
      curr <- params$layers$labs$caption
      curr <- if (is.null(curr)) '' else paste0(curr, "\n")
      
      interp <- with(params, glue(switch(
        EXPR = paste0(test, "_", is.null(stat.by)),
        emmeans_TRUE   = "trendline mean {alt} {mu}",
        emmeans_FALSE  = "trendline means are different",
        emtrends_TRUE  = "trendline slope {alt} {mu}",
        emtrends_FALSE = "trendline slopes are different" )))
      
      if (nrow(params$.plot_attrs[['stats']]) <= 1) {
        meth <- ""
      } else {
        meth <- switch(
          EXPR = params$p.adj,
          holm       = "Holm",                  # (1979)
          hochberg   = "Hochberg",              # (1988)
          hommel     = "Hommel",                # (1988)
          BH         = "Benjamini & Hochberg",  # (1995)
          fdr        = "Benjamini & Hochberg",  # (1995)
          BY         = "Benjamini & Yekutieli", # (2001)
          bonferroni = "Bonferroni",
          none       = "no" )
        meth <- paste0("\n", meth, " FDR correction in use.")
      }
      
      return(glue("{curr}Low p-value indicates {interp}.{meth}"))
      
    }))
  }
  
  
  
  
  # if (!is.null(stats[['fit']])) {
  #
  #   stats_text <- sprintf(
  #     fmt = "*p* = %s; *R<sup>2</sup>* = %s; *F* = %s",
  #     format(stats[['fit']][['p.value']],   digits=3),
  #     format(stats[['fit']][['r.squared']], digits=3),
  #     format(stats[['fit']][['statistic']], digits=3) )
  #
  #
  #
  #   #________________________________________________________
  #   # Add caption describing the model/formula.
  #   #________________________________________________________
  #   model_cmd <- local({
  #
  #     fun  <- model[[1]]
  #     args <- model[[2]]
  #
  #     fm <- capture.output(args[['formula']])[[1]]
  #     for (i in predictors)
  #       fm %<>% paste(sep = " + ", capture.output(as.symbol(i)))
  #     args[['formula']] <- structure(fm, display = fm)
  #
  #     str <- sprintf("%s(%s)", attr(fun, 'display'), as.args(args, fun = fun))
  #
  #     # Ensure that nothing in the formula is interpreted as markdown syntax.
  #     entities <- c(
  #       setNames(paste0("&#", 33:42,    ";"), strsplit("!\"#$%&'()*", "")[[1]]),
  #       setNames(paste0("&#", c(60,62), ";"), strsplit("<>", "")[[1]]),
  #       setNames(paste0("&#", 91:96,    ";"), strsplit("[\\]^_`", "")[[1]]),
  #       setNames(paste0("&#", 123:126,  ";"), strsplit("{|}~", "")[[1]]) )
  #     for (i in seq_along(entities))
  #       str <- gsub(str, pattern = names(entities)[[i]], replacement = entities[[i]], fixed = TRUE)
  #
  #     return (str)
  #   })
  #
  #   methods_text <- ifelse(
  #     test = isFALSE(ci),
  #     yes  = sprintf("Curve fitted using %s", model_cmd),
  #     no   = sprintf("Curve and %g%% CI fitted using %s", level * 100, model_cmd) )
  #
  #   subtitle <- sprintf("%s<br><span style='font-size:9pt'>%s</span>", stats_text[[1]], methods_text)
  #   set_layer(params, 'labs',  subtitle = subtitle)
  #   set_layer(params, 'theme', plot.subtitle = element_markdown(size = 11, lineheight = 1.2))
  #
  # }
  
  
  
  return (params)
}
