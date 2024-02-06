

init_corrplot_layers <- function (params = parent.frame()) {
  
  stopifnot(is_bare_environment(params))
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  with(params, {
    stopifnot(is_tibble(.ggdata))
    validate_formula()
    validate_var_range('level', range = c(0.5, 1))
    validate_var_choices('engine', choices = c("lm", "local"))
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
      's' = "scatter", 'n' = "name" ))
  
  
  
  #________________________________________________________
  # Merge trend and confidence into a single layer.
  #________________________________________________________
  if (any(c('trend', 'confidence') %in% layer_names)) {
    
    add_layer(params, 'smooth')
    set_layer(params, 'smooth', formula = params$formula)
    
    if (!'trend'      %in% layer_names) set_layer(params, 'smooth', color = NA)
    if (!'confidence' %in% layer_names) set_layer(params, 'smooth', se = FALSE)
    
    if (params$engine == "lm")
      set_layer(params, 'smooth', method = params$engine)
    
    if ('confidence' %in% layer_names && params$level != 0.95)
      set_layer(params, 'smooth', level = params$level)
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
    
    
    if (has_layer(params, 'smooth')) {
      
      if ('trend' %in% layer_names)
        set_layer(params, 'smooth', 'mapping|color' = color.by)
      
      if ('confidence' %in% layer_names) {
        
        set_layer(params, 'smooth', 'mapping|fill' = color.by)
      
        if (!is_null(params$layers[['color']][['values']]))
          set_layer(
            params = params, 
            layer  = 'fill', 
            'values' = params$layers[['color']][['values']] )
      }
    }
    
  }
  
  
  
  #________________________________________________________
  # Display sample names at x,y coordinates.
  #________________________________________________________
  if (has_layer(params, 'name'))
    set_layer(params, 'name', 'mapping|label' = ".sample")
  
  
  
  
  
  #________________________________________________________
  # All layers will have a common `x` and `y` field name.
  #________________________________________________________
  set_layer(
    params = params, 
    layer  = 'ggplot', 
    'mapping|x' = params$.xcol, 
    'mapping|y' = params$.ycol )
  
  
  return (invisible(params))
}
