

#________________________________________________________
# Computes p-values for categorical differences.
#________________________________________________________
corrplot_stats <- function (params) {
  
  stopifnot(is_bare_environment(params))
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  with(params, {
    validate_bool('caption')
    validate_var_range('level',   range = c(0, 1))
    validate_var_choices('p.adj', choices = p.adjust.methods)
    validate_var_choices('test',  choices = c(
      'none', 'fit', 'terms', 'means', 'trends', 'pw_means', 'pw_trends' ))
  })
  
  
  
  #________________________________________________________
  # Omit stats from this plot.
  #________________________________________________________
  if (params$test == 'none')
    return (params)
  
  if (!params$test %in% c('fit', 'terms')) {
    group <- names(params$color.by)
    if (is.null(group))                       return (params)
    if (nlevels(params$.ggdata[[group]]) < 2) return (params)
    remove("group")
  }
  
  
  
  #________________________________________________________
  # If anything goes wrong, skip stats.
  #________________________________________________________
  tryCatch(
    error = function (e) { warning(e); return (params) }, 
    expr  = with(params, {
      .plot_attrs[['stats']] <- stats_table(
        test     = test, 
        df       = .ggdata, 
        stat.by  = names(color.by), 
        resp     = .ycol, 
        regr     = .xcol, 
        model    = model, 
        split.by = facet.by, 
        level    = level, 
        p.adj    = p.adj )}))
  
  
  
  #________________________________________________________
  # Use p.top to retain only most significant taxa.
  #________________________________________________________
  if (isTRUE(is.finite(params$p.top))) {
    
    apply_p.top(params)
    
    if (plyr::empty(params$.ggdata))
      return (params)
  }
  
  
  
  #________________________________________________________
  # Show p-values on the plot. One per facet.
  #________________________________________________________
  
  with(params, {
    attr(.ggdata, 'stat_labels') <- plyr::ddply(
      .data      = .plot_attrs[['stats']], 
      .variables = ply_cols(facet.by), 
      .fun       = function (df) {
        p <- min(df[['.adj.p']], Inf, na.rm = TRUE)
        if (is.infinite(p)) return (NULL)
        fmt <- ifelse(nrow(df) == 1, "%s: *p* = %s", "%s: min(*p*) = %s")
        data.frame(.label = sprintf(fmt, test, format(p, digits=2)))
      }) %>% as_tibble()
  })
  
  
  set_layer(
    params = params, 
    layer  = 'stats_ggtext',
    'size'          = 3,
    'fill'          = "#FFFFFFBB",
    'width'         = NULL,
    'box.size'      = 0,
    'box.padding'   = unit(c(3,6,3,6), "pt"),
    'hjust'         = -0.1,
    'vjust'         = 1.2,
    'mapping|x'     = -Inf,
    'mapping|y'     =  Inf,
    'mapping|label' = ".label" )
  
  set_layer(
    params = params, 
    layer  = 'yaxis',
    'expand' = c(0, 0, .1, 0) )
    
  
  
  
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
  #     str <- sprintf("%s(%s)", attr(fun, "fn", exact = TRUE), as.args(args, fun = fun))
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
