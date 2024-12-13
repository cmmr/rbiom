
#________________________________________________________
# Assemble faceting formula and attach nrow/ncol/etc attributes
#________________________________________________________

plot_facets <- function (params) {
  
  stopifnot(is_bare_environment(params))
  
  ggdata   <- params$.ggdata
  facet.by <- params$facet.by
  
  
  #________________________________________________________
  # No facets
  #________________________________________________________
  if (is_null(facet.by) || plyr::empty(ggdata)) {
    
    del_layer(params, 'facet')
    
    params$.plot_attrs$facet.nrow  <- 1
    params$.plot_attrs$facet.ncol  <- 1
    params$.plot_attrs$facet.count <- 1
    
    return (invisible(params))
  }
  
  
  
  
  #________________________________________________________
  # Use user's or auto-selected free_x/y scales.
  #________________________________________________________
  set_layer(params, 'facet', scales = local({
    
    
    #________________________________________________________
    # User has explicitly defined scales.
    #________________________________________________________
    scales <- params$layers[['facet']][['scales']]
    
    if (!is.null(scales)) {
      
      stopifnot(is_string(scales, c('fixed', 'free', 'free_x', 'free_y')))
      
      params$.free_x <- scales %in% c('free', 'free_x')
      params$.free_y <- scales %in% c('free', 'free_y')
      
      return (scales)
    }
    
    
    
    #________________________________________________________
    # Examine the variance in min/max per facet.
    #________________________________________________________
    xmode <- params$.xmode
    xcol  <- params$.xcol
    ycol  <- params$.ycol
    stopifnot(is_scalar_character(xmode))
    stopifnot(is_scalar_character(xcol) || is.null(xcol))
    stopifnot(is_scalar_character(ycol))
    
    df <- if (xmode == "numeric") {
      
      coords <- ggdata
      
      # For corrplots with only trendlines.
      if (!any(has_layer(params, c('point', 'residual', 'name')))) { 
        if (has_layer(params, 'confidence')) {
          ycol   <- ".ymax"
          coords <- attr(ggdata, 'fit')
        } else if (has_layer(params, 'trend')) {
          coords <- attr(ggdata, 'fit')
        }
      }
      
      plyr::ddply(coords, ply_cols(facet.by), function (d) {
        
        if (is.null(xcol)) d[[xcol <- '.x']] <- 0
        
        data.frame(
          x_min = min(c(d[[xcol]], 0), na.rm = TRUE),
          x_max = max(c(d[[xcol]], 0), na.rm = TRUE),
          x_pct = 1,
          y_min = min(c(d[[ycol]], 0), na.rm = TRUE),
          y_max = max(c(d[[ycol]], 0), na.rm = TRUE) )
      })
      
    } else {
      
      plyr::ddply(ggdata, ply_cols(facet.by), function (d) {
        
        if (is.null(xcol)) d[[xcol <- '.x']] <- 0
        
        data.frame(
          x_min = 0,
          x_max = 1,
          x_pct = length(unique(d[[xcol]])) / nlevels(d[[xcol]]),
          y_min = min(c(d[[ycol]], 0), na.rm = TRUE),
          y_max = max(c(d[[ycol]], 0), na.rm = TRUE) )
      })
      
    }
    
    
    
    #________________________________________________________
    # Automatically free the scales if ggdata needs it.
    #________________________________________________________
    
    params$.free_x <- with(df, any(
      isTRUE(params$.free_x),
      min(x_pct) < 0.9,
      diff(range(x_max)) > abs(max(x_max) / 2),
      diff(range(x_min)) > abs(min(x_min) / 2) ))
    
    params$.free_y <- with(df, any(
      isTRUE(params$.free_y),
      diff(range(y_max)) > abs(max(y_max) / 2),
      diff(range(y_min)) > abs(min(y_min) / 2) ))
    
    return (switch(
      paste(params$.free_x, params$.free_y),
      'TRUE TRUE'   = "free", 
      'TRUE FALSE'  = "free_x", 
      'FALSE TRUE'  = "free_y",
      'FALSE FALSE' = "fixed" ))
    
  }))
  
  
  
  
  
  #________________________________________________________
  # Facet arguments
  #________________________________________________________
  args <- list()
  
  # if (length(facet.by) == 1) {
  #   
  #   args[['facets']] <- sprintf(
  #     fmt = "~ %s", 
  #     capture.output(as.name(facet.by)) ) %>%
  #     as.formula()
  #   
  #   nfacets <- length(unique(data[[facet.by]]))
  #   autodim <- ggplot2::wrap_dims(
  #     n    = nfacets,
  #     nrow = layers[['facet']][['nrow']],
  #     ncol = layers[['facet']][['ncol']] )
  #   
  #   attr(layers, 'facet.nrow')  <- autodim[[1]]
  #   attr(layers, 'facet.ncol')  <- autodim[[2]]
  #   attr(layers, 'facet.count') <- nfacets
  #   
  # } else 
    
  if (length(facet.by) == 2) {
    
    args[['rows']] <- sprintf(
      fmt = "%s ~ %s",
      capture.output(as.name(facet.by[[2]])),
      capture.output(as.name(facet.by[[1]])) ) %>%
      as.formula()
    
    f_rows <- length(unique(ggdata[[facet.by[[2]]]]))
    f_cols <- length(unique(ggdata[[facet.by[[1]]]]))
    
    params$.plot_attrs$facet.nrow  <- f_rows
    params$.plot_attrs$facet.ncol  <- f_cols
    params$.plot_attrs$facet.count <- f_rows * f_cols
    
  } else {
    
    nfacets <- ggdata[,facet.by,drop=FALSE] %>%
      apply(1L, paste, collapse="|@#|") %>%
      unique() %>%
      length()
    
    autodim <- ggplot2::wrap_dims(
      n    = nfacets,
      nrow = params$layers[['facet']][['nrow']],
      ncol = params$layers[['facet']][['ncol']] )
    
    args[['facets']] <- do.call(ggplot2::vars, lapply(facet.by, as.name))
    
    attr(args[['facets']], 'display') <- facet.by %>%
        sapply(as.name) %>%
        sapply(capture.output) %>%
        paste(collapse = ", ") %>%
        sprintf(fmt = "vars(%s)")
    
    params$.plot_attrs$facet.nrow  <- autodim[[1]]
    params$.plot_attrs$facet.ncol  <- autodim[[2]]
    params$.plot_attrs$facet.count <- nfacets
    
  }
  
  set_layer(params, 'facet', args)
  
  # To enable %>% chaining
  return (invisible(params))
}



