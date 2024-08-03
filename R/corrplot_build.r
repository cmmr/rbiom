



corrplot_build <- function (params) {
  
  stopifnot(is_bare_environment(params))
  
  
  with(params, {
    
    
    #________________________________________________________
    # Consistent naming for plot_* functions.
    #________________________________________________________
    .ggdata <- as_tibble(df)
    .xcol   <- x
    .ycol   <- y
    .xmode  <- "numeric"
    remove("df", "x", "y")
    
    
    #________________________________________________________
    # Clean up the x- and y- axis values.
    #________________________________________________________
    .ggdata <- .ggdata[is.finite(.ggdata[[.xcol]]),]
    .ggdata <- .ggdata[is.finite(.ggdata[[.ycol]]),]
    
    
    #________________________________________________________
    # Set lower/upper bounds for confidence interval range.
    #________________________________________________________
    if (is.null(.dots$ci.min)) .dots$ci.min <- -Inf
    if (is.null(.dots$ci.max)) .dots$ci.max <-  Inf
    
  })
  
  
  
  #________________________________________________________
  # Convert user's `layers` spec to layer names.
  #________________________________________________________
  init_layers(
    params  = params,
    choices = c(
      't' = "trend", 'c' = "confidence", 'p' = "point", 
      'n' = "name", 'r' = "residual" ) )
  
  
  
  #________________________________________________________
  # All layers will have a common `x` and `y` field name.
  #________________________________________________________
  set_layer(
    params = params, 
    layer  = 'ggplot', 
    'mapping|x' = params$.xcol, 
    'mapping|y' = params$.ycol )
  
  
  
  
  
  #________________________________________________________
  # Display sample names at x,y coordinates.
  #________________________________________________________
  if (has_layer(params, 'name'))
    set_layer(params, 'name', 'mapping|label' = ".sample")
  
  
  
  #________________________________________________________
  # Do the model fitting just once.
  #________________________________________________________
  if (any(has_layer(params, c('confidence', 'trend', 'residual'))) || isTRUE(params$check))
    params$.models <- local({
      
      models <- with(params, {
        .ggdata %>% 
          plyr::dlply(ply_cols(c(stat.by, facet.by)), function (d) {
            
            data.frame(x = d[[.xcol]], y = d[[.ycol]]) %>%
              stats_fit_model(fit = fit, regr = 'x', resp = 'y') %>% 
              suppressWarnings() %>% 
              tryCatch(error = function (e) NULL)
            
          })
      })
      
      # Drop models that returned NULL
      attrs  <- attributes(models)
      keep   <- !sapply(models, is.null)
      models <- models[keep] %>%
        aa(split_type   = attrs$split_type) %>%
        aa(split_labels = attrs$split_labels[keep,,drop=FALSE])
      
      return (models)
    })
  
  
  
  #________________________________________________________
  # Residual segments from trendline to point.
  #________________________________________________________
  
  if (has_layer(params, 'residual') || isTRUE(params$check)) {
    
    if (has_layer(params, 'residual'))
      set_layer(params, 'residual', 'mapping|yend' = ".yend")
    
    
    with(params, {
      
      attr(.ggdata, 'residual') <- plyr::ldply(.models, function (m) {
        
        mf <- stats::model.frame(m)
        mf$x %<>% signif(digits = 10)
        at <- list(x = unique(mf$x))
        
        emmeans::emmeans(m, specs = 'x', at = at) %>% 
          summary() %>% 
          with(tibble(x = x, .se = SE, .fitted = emmean)) %>%
          dplyr::right_join(mf, by = 'x') %>% 
          dplyr::mutate(
            .residual = y - .fitted,
            .yend     = .fitted %>% pmax(.dots$ci.min) %>% pmin(.dots$ci.max), 
            .sd       = sqrt(nrow(mf)) * .se,
            .std.res  = .residual / .sd ) %>% 
          dplyr::select(!!.xcol := x, !!.ycol := y, dplyr::starts_with('.')) %>% 
          dplyr::relocate(.yend, .se, .sd, .after = ".std.res") %>% 
          suppressWarnings() %>% 
          tryCatch(error = function (e) tibble())
        
      }) %>% 
        as_tibble()
      
    })
    
    
    if (plyr::empty(attr(params$.ggdata, 'residual')))
      del_layer(params, 'residual')
    
  }
  
  
  
  
  #________________________________________________________
  # Means (and optionally CI) at 100 x-axis locations.
  #________________________________________________________
  
  if (any(has_layer(params, c('confidence', 'trend')))) {
    
    
    if (has_layer(params, 'confidence'))
      set_layer(
        params = params, 
        layer  = 'confidence', 
        'mapping|ymin' = ".ymin", 
        'mapping|ymax' = ".ymax",  
        'alpha'        = 0.25,
        'linewidth'    = 0.2 )
    
    
    with(params, {
      
      # For confidence, set infer = TRUE
      if (hasName(layers, 'confidence')) {
        
        attr(.ggdata, 'fit') <- as_tibble(plyr::ldply(.models, function (m) {
          
          dr <- range(stats::model.frame(m)$x %||% exp(stats::model.frame(m)$`log(x)`))
          at <- list(x = seq(from = dr[[1]], to = dr[[2]], length.out = 100))
          
          emmeans::emmeans(m, specs = 'x', at = at, infer = TRUE, level = level) %>% 
            summary() %>% 
            with(tibble(
              !!.xcol := x, 
              !!.ycol := emmean   %>% pmax(.dots$ci.min) %>% pmin(.dots$ci.max), 
              .ymin    = lower.CL %>% pmax(.dots$ci.min), 
              .ymax    = upper.CL %>% pmin(.dots$ci.max) )) %>%
            suppressWarnings() %>% 
            tryCatch(error = function (e) tibble() )
          
        }))
        
        
      } else {
        
        attr(.ggdata, 'fit') <- as_tibble(plyr::ldply(.models, function (m) {
          
          dr <- range(stats::model.frame(m)$x %||% exp(stats::model.frame(m)$`log(x)`))
          at <- list(x = seq(from = dr[[1]], to = dr[[2]], length.out = 100))
          
          if (fit == 'lm') at %<>% within(x %<>% range()) # list(x = c(min,max))
          
          emmeans::emmeans(m, specs = 'x', at = at) %>% 
            summary() %>% 
            with(tibble(
              !!.xcol := x, 
              !!.ycol := emmean %>% pmax(.dots$ci.min) %>% pmin(.dots$ci.max) )) %>%
            suppressWarnings() %>% 
            tryCatch(error = function (e) tibble())
          
        }))
        
      }
      
    })
    
    
    if (plyr::empty(attr(params$.ggdata, 'fit'))) {
      del_layer('confidence')
      del_layer('trend')
    }
    
  }
  
  
  
  #________________________________________________________
  # Rarefaction curve's vertical line.
  #________________________________________________________
  if (!is_null(params$rline))
    set_layer(
      params     = params, 
      layer      = 'vline', 
      xintercept = params$rline, 
      color      = "red", 
      linetype   = "dashed" )
  
  
  
  #________________________________________________________
  # Theming specific to a numeric x-axis.
  #________________________________________________________
  set_layer(params, 'theme', panel.grid.minor.x = element_blank())
  
  
  
  
  #________________________________________________________
  # Build the plot.
  #________________________________________________________
  p <- params %>% 
    corrplot_stats() %>% 
    plot_facets() %>% 
    plot_build()
  
  
  
  attr(p, 'models') <- params$.models
  
  
  return (p)
}




