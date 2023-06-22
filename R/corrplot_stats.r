

#________________________________________________________
# Computes p-values for categorical differences.
#________________________________________________________
corrplot_stats <- function (layers) {
  
  params <- attr(layers, 'params', exact = TRUE)
  xcol   <- attr(layers, 'xcol',   exact = TRUE)
  ycol   <- attr(layers, 'ycol',   exact = TRUE)
  ggdata <- attr(layers, 'data',   exact = TRUE)
  
  color.by  <- names(params[['color.by']])
  facet.by  <- params[['facet.by']]
  ci        <- params[['ci']]
  model     <- params[['model']]
  
  color.by   <- unique(setdiff(color.by, facet.by))
  facet.by   <- unique(facet.by)
  emm_specs  <- if (!is.null(color.by)) color.by else facet.by
  emm_by     <- if (!is.null(color.by)) facet.by else NULL
  predictors <- c(emm_specs, emm_by)
  
  #________________________________________________________
  # If anything goes wrong, skip stats.
  #________________________________________________________
  stats <- tryCatch(
    error = function (e) { warning(e); return (NULL) }, 
    expr  = local({
      
      if (is_null(model))             return (NULL)
      if (isFALSE(params[['stats']])) return (NULL)
      
      
      #________________________________________________________
      # Convert from geom_smooth formula to manual stats.
      #________________________________________________________
      m <- local({
        
        model_function <- model[[1]]
        model_args     <- model[[2]]
        
        model_formula <- model_args[['formula']]
        replacements  <- list(x = as.symbol(xcol), y = as.symbol(ycol))
        model_formula <- do.call(substitute, list(model_formula, replacements))
        
        model_formula <- as.formula(paste0(collapse = " + ", c(
          capture.output(model_formula),
          sapply(c(color.by, facet.by), function (x) capture.output(as.symbol(x))),
          "0" )))
        
        attr(ggdata, 'display') <- "data"
        model_args[['data']]    <- ggdata
        model_args[['formula']] <- model_formula
        
        fn <- attr(model_function, 'fn', exact = TRUE)
        m  <- do.call(model_function, model_args, envir = baseenv())
        
        attr(m, 'display') <- "m"
        attr(m, 'cmd')     <- sprintf("%s(%s)", fn, as.args(model_args, fun = model_function))
        attr(m, 'history') <- c(attr(ggdata, 'history'), sprintf("m <- %s", attr(m, 'cmd')))
        
        return (m)
      })
      
      
      #________________________________________________________
      # Basic stats for all models.
      #________________________________________________________
      ggdata <- augment(m, data = ggdata)
      attr(ggdata, 'history') <- c(attr(m, 'history'), "data <- broom::augment(m, data = data)")
      
      stats <- list(
        'fit'     = glance(x = m),
        'terms'   = tidy(x = m) )
      
      
      #________________________________________________________
      # Stats that come into play for color.by / facet.by
      #________________________________________________________
      if (length(emm_specs) > 0) {
        
        emm_options      <- list(ci = ifelse(is.numeric(ci), ci / 100, 0.95), infer = TRUE)
        emm_args         <- list(object = m, specs = emm_specs, by = emm_by, options = emm_options)
        emm              <- do.call(emmeans::emmeans, emm_args) %>% aa(display = "emm")
        attr(emm, 'cmd') <- sprintf("emm <- emmeans::emmeans(%s)", as.args(emm_args, fun = emmeans::emmeans))
        
        stats[['groupwise']] <- summary(object = emm)
        stats[['pairwise']]  <- pairs(x = emm)
        stats[['effect']]    <- eff_size(object = emm, sigma = sigma(m), edf = df.residual(m))
        
        stats <- stats[c(3,4,1,2,5)]
      }
      
      
      stats <- lapply(stats, as.data.frame)
      return (stats)
    }))
  
  
  
  if (!is.null(stats[['fit']])) {
    
    stats_text <- sprintf(
      fmt = "*p* = %s; *R<sup>2</sup>* = %s; *F* = %s",
      format(stats[['fit']][['p.value']],   digits=3), 
      format(stats[['fit']][['r.squared']], digits=3), 
      format(stats[['fit']][['statistic']], digits=3) )
    
    setLayer("labs",  subtitle      = stats_text[[1]])
    setLayer("theme", plot.subtitle = element_markdown())
    
    

    #________________________________________________________
    # Add caption describing the model/formula.
    #________________________________________________________
    model_cmd <- local({
      
      fun  <- model[[1]]
      args <- model[[2]]
      
      fm <- capture.output(args[['formula']])[[1]]
      for (i in predictors)
        fm %<>% paste(sep = " + ", capture.output(as.symbol(i)))
      args[['formula']] <- structure(fm, display = fm)
      
      sprintf("%s(%s)", attr(fun, "fn", exact = TRUE), as.args(args, fun = fun))
    })
    
    caption_text <- ifelse(
      test = isFALSE(ci),
      yes  = sprintf("Curve fitted using %s", model_cmd),
      no   = sprintf("Curve and %g%% CI fitted using %s", ci, model_cmd) )
    setLayer("labs",  caption = caption_text)
    setLayer("theme", plot.caption = element_text(face = "italic"))
    
    
  }
  
  
  params[['stats']] <- NULL
  params[['model']] <- NULL
  
  attr(layers, 'data')   <- ggdata
  attr(layers, "stats")  <- stats
  attr(layers, 'params') <- params
  
  
  return (layers)
}
