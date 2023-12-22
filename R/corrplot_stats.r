

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
    stat.by <- names(params$color.by)
    
    if (is.null(stat.by)) {
      if (startsWith(params$test, 'pw_')) return (params)
    } else {
      if (nlevels(params$.ggdata[[stat.by]]) < 2) return (params)
    }
    
    remove("stat.by")
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
        p.adj    = p.adj )
      
      
      #________________________________________________________
      # Add biom subsetting code before stats code.
      #________________________________________________________
      if (!is.null(attr(.plot_attrs[['stats']], 'code', exact = TRUE)))
        if (exists(".subset_code", inherits = FALSE) && !is.null(.subset_code))
          attr(.plot_attrs[['stats']], 'code') %<>% paste0(.subset_code, "\n\n", .) %>%
        add_class("rbiom_code")
      
    }))
  
  
  
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
        fmt <- ifelse(nrow(df) == 1, "*p* = %s", "min(*p*) = %s")
        data.frame(.label = sprintf(fmt, format(p, digits=2)))
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
  
  
  
  if (isTRUE(params$caption)) {
    
    # element_markdown
    set_layer(params, 'theme', plot.caption = element_text(size = 9, face = "italic"))
    set_layer(params, 'labs',  caption      = local({
      
      test <- switch(
        EXPR = params$test,
        fit       = "Pr(fit is expected by chance)",
        terms     = "Pr(model term's estimate is zero)",
        means     = "Pr(mean is zero)",
        trends    = "Pr(trend slope is flat)",
        pw_means  = "Pr(means are equal)",
        pw_trends = "Pr(trend slopes are equal)" )
      
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
      
      return(glue("{test}, with {meth} FDR correction."))
      
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
