

init_corrplot_layers <- function (params = parent.frame()) {
  
  stopifnot(is_bare_environment(params))
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  with(params, {
    stopifnot(is_tibble(.ggdata))
    validate_model()
  })
  
  
  #________________________________________________________
  # How to interpret ggdata.
  #________________________________________________________
  params$.xcol  %<>% if.null(params$x)
  params$.ycol  %<>% if.null(attr(params$.ggdata, "response", exact = TRUE))
  params$.xmode <- "numeric"
  stopifnot(is_scalar_character(params$.xcol))
  stopifnot(is_scalar_character(params$.ycol))
  
  
  #________________________________________________________
  # Clean up the y-axis values.
  #________________________________________________________
  with(params, {
    .ggdata <- .ggdata[is.finite(.ggdata[[.xcol]]),]
    .ggdata <- .ggdata[is.finite(.ggdata[[.ycol]]),]
  })
  
  
  
  #________________________________________________________
  # Convert user's `layers` spec to layer names.
  #________________________________________________________
  
  layer_names <- init_layers(
    params  = params,
    no_init = c("trend", "confidence"),
    choices = c(
      't' = "trend",   'c' = "confidence", 
      's' = "scatter", 'r' = "residual",   'n' = "name") )
  
  
  
  # Merge trend and confidence into a single layer.
  if (any(c('trend', 'confidence') %in% layer_names)) {
    add_layer(params, 'trend')
    if (!'trend'      %in% layer_names) set_layer(params, 'trend', color = NA)
    if (!'confidence' %in% layer_names) set_layer(params, 'trend', se = FALSE)
  }
  
  
  
  
  #________________________________________________________
  # Theming specific to a numeric x-axis.
  #________________________________________________________
  set_layer(params, 'theme', panel.grid.minor.x = element_blank())
  
  
  #________________________________________________________
  # aes() mappings for color/fill
  #________________________________________________________
  if (has_layer(params, 'color')) {
    
    color.by <- names(params$color.by)
    
    if (has_layer(params, 'scatter'))
      set_layer(params, 'scatter', 'mapping|color' = color.by)
    
    if (has_layer(params, 'residual'))
      set_layer(params, 'residual', 'mapping|color' = color.by)
    
    if (has_layer(params, 'trend')) {
      set_layer(params, 'trend', 'mapping|color' = color.by)
      
      if ('confidence' %in% layer_names) {
        set_layer(params, 'trend', 'mapping|fill' = color.by)
        
        if (!is_null(params$layers[['color']][['values']]))
          set_layer(
            params = params, 
            layer  = 'fill', 
            'values' = params$layers[['color']][['values']] )
      }
    }
  }
  
  
  
  #________________________________________________________
  # Define the curve and confidence interval.
  #________________________________________________________
  if (has_layer(params, 'trend')) {
    
    model <- params$model
    
    
    #________________________________________________________
    # The method for geom_smooth.
    #________________________________________________________
    method <- with(params$model, structure(
      .Data   = function (data, ...) do.call(fun, c(args, list(data = data))),
      display = sprintf(
        fmt = "function (data, ...) %s(%s, data = data)",
        attr(fun, 'fn'), as.args(args, fun = fun) )))
    
    set_layer(
      params = params, 
      layer  = 'trend', 
      'method'  = method, 
      'formula' = params$model[['args']][['formula']] )
    
    if (!isTRUE(params$level == 0.95)) # 0.95 is default
      set_layer(params, 'trend', level = params$level)
    
    
    
    #________________________________________________________
    # Fade out the points when a curve is fitted.
    #________________________________________________________
    if (has_layer(params, 'scatter'))
      set_layer(params, 'scatter', size = 0.2, alpha = 0.5)
    
    
    remove("method")
  }
  
  
  
  
  
  #________________________________________________________
  # Connect points to trendline with vertical lines.
  #________________________________________________________
  if (has_layer(params, 'residual')) {
    
    with(params, {
      attr(.ggdata, 'residual') <- stats_table(
          test     = "predict", 
          df       = .ggdata, 
          stat.by  = names(color.by), 
          resp     = .ycol, 
          regr     = .xcol, 
          model    = model, 
          split.by = facet.by ) %>%
        keep_cols('.sample', .xcol, .ycol, '.fitted', names(color.by), facet.by)
    })
    
    set_layer(
      params = params, 
      layer  = 'residual',
      'mapping|x'    = params$.xcol,
      'mapping|xend' = params$.xcol,
      'mapping|y'    = params$.ycol,
      'mapping|yend' = '.fitted' )
  }
  
  
  
  
  
  #________________________________________________________
  # Display sample names at x,y coordinates.
  #________________________________________________________
  if (has_layer(params, 'name'))
    set_layer(params, 'name', 'mapping|label' = ".sample")
  
  
  
  
  
  #________________________________________________________
  # Add layers will have a common `x` and `y` field name.
  #________________________________________________________
  set_layer(
    params = params, 
    layer  = 'ggplot', 
    'mapping|x' = params$.xcol, 
    'mapping|y' = params$.ycol )
  
  
  return (invisible(params))
}
