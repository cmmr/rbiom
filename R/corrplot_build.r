

corrplot_build <- function (layers) {
  
  params <- attr(layers, 'params', exact = TRUE)
  
  ci       <- params[['ci']]
  color.by <- names(params[['color.by']])
  
  
  #________________________________________________________
  # Theming specific to a numeric x-axis.
  #________________________________________________________
  setLayer("theme", panel.grid.minor.x = element_blank())
  
  
  #________________________________________________________
  # aes() mapping for color/fill
  #________________________________________________________
  if (is_scalar_character(color.by)) {
    
    if (hasName(layers, 'scatter'))
      setLayer("scatter", "mapping|color" = color.by)
    
    if (hasName(layers, 'residual'))
      setLayer("residual", "mapping|color" = color.by)
    
    if (hasName(layers, 'trend')) {
      setLayer("trend", "mapping|color" = color.by)
      
      if (!isFALSE(ci)) {
        setLayer("trend", "mapping|fill" = color.by)
        
        if (!is_null(layers[["color"]][["values"]]))
          setLayer("fill", values = layers[["color"]][["values"]])
      }
    }
  }
  
  
  
  #________________________________________________________
  # Define the curve and confidence interval.
  #________________________________________________________
  if (hasLayer("trend")) {
    
    model <- params[['model']]
    
    
    #________________________________________________________
    # Predefined regression models
    #________________________________________________________
    if (is_scalar_character(model)) {
      
      model <- local({
        
        models <- list(
          'lm'  = list("stats::lm", list(formula = y ~ x      )), 
          'log' = list("stats::lm", list(formula = y ~ log(x) )), 
          'gam' = list("mgcv::gam", list(formula = y ~ s(x, bs = "cs"), method = "REML" )) )
        
        i <- pmatch(tolower(model), names(models), NA)
        if (isTRUE(is.na(i)))
          stop("'", model, "' is not a predefined regression model.")
        
        return (models[[i]])
      })
    }
    
    
    #________________________________________________________
    # Sanity check model function.
    #________________________________________________________
    stopifnot(is_list(model))
    stopifnot(identical(length(model), 2L))
    stopifnot(is_list(model[[2]]))
    stopifnot(is_scalar_character(model[[1]]) && !is_na(model[[1]]))
    stopifnot(is_formula(model[[2]][['formula']]))
    
    model[[1]] <- local({
      fn  <- strsplit(model[[1]], '::', fixed = TRUE)[[1]]
      fun <- if (length(fn) == 1) get(fn[[1]]) else getFromNamespace(fn[[2]], fn[[1]])
      stopifnot(is_function(fun))
      stopifnot(all(c('formula', 'data') %in% formalArgs(fun)))
      attr(fun, 'fn') <- model[[1]]
      return (fun)
    })
    
    
    #________________________________________________________
    # The method for geom_smooth.
    #________________________________________________________
    method <- function (data, ...) do.call(model[[1]], c(model[[2]], list(data = data)))
    attr(method, 'display') <- sprintf(
      fmt = "function (data, ...) %s(%s, data = data)",
      attr(model[[1]], 'fn', exact = TRUE),
      as.args(model[[2]], fun = model[[1]]) )
    setLayer("trend", method = method)
    
    
    if (is.numeric(ci)) setLayer("trend", level = ci / 100)
    params[grep("(^|\\.)(se|level)$", names(params))] <- NULL
    
    
    
    #________________________________________________________
    # Fade out the points when a curve is fitted.
    #________________________________________________________
    if (hasLayer("scatter"))
      setLayer("scatter", size = 0.2, alpha = 0.5)
    
    
    params[['model']] <- model
    
  } else {
    params[['model']] <- NULL
  }
  
  
  
  
  
  
  attr(layers, 'params') <- params
  
  
  #________________________________________________________
  # Build the plot.
  #________________________________________________________
  p <- layers %>% 
    plot_facets() %>% 
    corrplot_stats() %>%
    plot_build()
  
  
  return (p)
}
