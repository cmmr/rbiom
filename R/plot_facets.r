
#--------------------------------------------------------------
# Assemble faceting formula and attach nrow/ncol/etc attributes
#--------------------------------------------------------------

plot_facets <- function (layers) {
  
  data   <- attr(layers, "data",   exact = TRUE)
  params <- attr(layers, 'params', exact = TRUE)
  
  facet.by <- params[['facet.by']]
  
  
  #--------------------------------------------------------------
  # No facets
  #--------------------------------------------------------------
  if (is.null(facet.by)) {
    
    layers[['facet']] <- NULL
    
    attr(layers, 'facet.nrow')  <- 1
    attr(layers, 'facet.ncol')  <- 1
    attr(layers, 'facet.count') <- 1
    
    return (layers)
  }
  
  
  #--------------------------------------------------------------
  # Facet arguments
  #--------------------------------------------------------------
  args <- list()
  
  if (length(facet.by) == 1) {
    
    args[['facets']] <- sprintf(
      fmt = "~ %s", 
      capture.output(as.name(facet.by)) ) %>%
      as.formula()
    
    nfacets <- length(unique(data[[facet.by]]))
    autodim <- ggplot2::wrap_dims(
      n    = nfacets,
      nrow = layers[['facet']][['nrow']],
      ncol = layers[['facet']][['ncol']] )
    
    attr(layers, 'facet.nrow')  <- autodim[[1]]
    attr(layers, 'facet.ncol')  <- autodim[[2]]
    attr(layers, 'facet.count') <- nfacets
    
  } else if (length(facet.by) > 1) {
    
    args[['rows']] <- sprintf(
      fmt = "%s ~ %s",
      capture.output(as.name(facet.by[[2]])),
      capture.output(as.name(facet.by[[1]])) ) %>%
      as.formula()
    
    f_rows <- length(unique(data[[facet.by[[2]]]]))
    f_cols <- length(unique(data[[facet.by[[1]]]]))
    
    attr(layers, 'facet.nrow')  <- f_rows
    attr(layers, 'facet.ncol')  <- f_cols
    attr(layers, 'facet.count') <- f_rows * f_cols
    
  }
  
  setLayer(layer = "facet", args)
  
  return (layers)
}



#--------------------------------------------------------------
# Specify individual y-axis limits and breaks for each facet.
#--------------------------------------------------------------

boxplot_facets <- function (layers) {
  
  ggdata   <- attr(layers, 'data', exact = TRUE)
  params   <- attr(layers, 'params', exact = TRUE)
  facet.by <- params[['facet.by']]
  
  
  #--------------------------------------------------------------
  # When facet.by=NULL, just pass layers on to plot_facets()
  #--------------------------------------------------------------
  if (!is.null(facet.by)) {
    
    
    #--------------------------------------------------------------
    # Pull facet's scales parameter from user args
    #--------------------------------------------------------------
    if (!hasLayer('facet'))
      initLayer("facet")
    
    
    #--------------------------------------------------------------
    # Automatically pick scales (free, fixed, free_x, or free_y)
    #--------------------------------------------------------------
    if (is.null(layers[['facet']][['scales']]) && nrow(ggdata) > 0) {
      
      xcol <- attr(layers, 'xcol', exact = TRUE)
      ycol <- attr(layers, 'ycol', exact = TRUE)
      
      plyby <- ply_cols(facet.by)
      
      
      # Don't get hung up on NA values in the plyby columns
      #--------------------------------------------------------------
      clean_df <- ggdata[,c(facet.by, xcol, ycol),drop=FALSE]
      clean_df <- clean_df[complete.cases(clean_df),,drop=FALSE]
      
      
      # At most 50% whitespace on average at the top of all plots
      #--------------------------------------------------------------
      y_maxs <- plyr::daply(clean_df, plyby, function (df) { max(df[[ycol]], na.rm = TRUE) })
      free_y <- isTRUE(mean(y_maxs / max(y_maxs, na.rm = TRUE)) < 0.5)
      
      
      # All x categories need to be present >50% of the time.
      #--------------------------------------------------------------
      x_cats <- sapply(levels(ggdata[[xcol]]), function (xval) {
        mean(plyr::daply(clean_df, plyby, function (df) { xval %in% df[[xcol]] }))
      })
      free_x <- isTRUE(any(x_cats <= 0.5))
      
      remove("clean_df")
      
      
      # Override free_x and free_y with layer attribute settings
      #--------------------------------------------------------------
      for (i in c('free_x', 'free_y'))
        if (!is.null(attr(layers, i, exact = TRUE)))
          assign(i, attr(layers, i, exact = TRUE))
      
      
      scales <- if (free_x && free_y) { "free" 
      } else    if (free_x)           { "free_x"
      } else    if (free_y)           { "free_y"
      } else                          { "fixed" }
      
      setLayer(layer = "facet", scales = scales)
    }
    
  }
  
  
  #--------------------------------------------------------------
  # Usual faceting configurations
  #--------------------------------------------------------------
  layers <- plot_facets(layers)
  
  
  #--------------------------------------------------------------
  # Rotate the x-axis tick mark labels to avoid overlap
  #--------------------------------------------------------------
  layer <- "theme"
  xcol  <- attr(layers, 'xcol', exact = TRUE)
  
  if (identical(xcol, '.all')) {
    
    if (isTRUE(params[['flip']])) {
      setLayer(
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank() )
      
    } else {
      setLayer(
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank() )
    }
    
  } else {
    
    if (identical(substr(xcol, 1, 1), ".")) {
      if (isTRUE(params[['flip']])) { setLayer(axis.title.y = element_blank())
      } else                        { setLayer(axis.title.x = element_blank()) }
    }
      
    
    ggdata     <- attr(layers, 'data',       exact = TRUE)
    facet.cols <- attr(layers, 'facet.ncol', exact = TRUE)
    color.by   <- params[['color.by']]
    xlab.angle <- params[['xlab.angle']] %||% "auto"
    
    # x-axis Text Angle
    if (tolower(xlab.angle) == 'auto') {
      if (!is.null(xcol) && !isTRUE(params[['flip']])) {
        charCount <- sum(nchar(unique(as.character(ggdata[[xcol]]))), na.rm = TRUE)
        charCount <- charCount * facet.cols
        if (charCount > 40)
          xlab.angle <- 30
        remove("charCount")
      }
    }
    
    if (xlab.angle == 90 || tolower(xlab.angle) == "vertical") {
      setLayer(axis.text.x = element_text(angle=-90, vjust=0.3, hjust=0))
      
    } else if (xlab.angle == 30 || tolower(xlab.angle) == "angled") {
      setLayer(axis.text.x = element_text(angle=-30, vjust=1, hjust=0))
      
      # Ensure long x-axis labels don't get truncated at the figure's edge
      rpad <- strwidth(tail(levels(ggdata[[xcol]]), 1), units="inches")
      rpad <- rpad * 0.8660254 # sin((90-30) * pi / 180) / sin(90 * pi / 180)
      if (!is.null(color.by))
        rpad <- rpad - max(c(
          strwidth(color.by, units="inches", cex=1.2) + .20,
          strwidth(levels(ggdata[[color.by]]), units="inches") + .52 ))
      rpad <- max(.1, signif(rpad, digits = 3))
      setLayer(plot.margin = as.cmd(unit(x=c(.1, rpad, .1, .1), units='inches'), list(rpad=rpad)))
    }
    
    remove("ggdata", "facet.cols", "color.by", "xlab.angle")
  }
  remove("layer", "xcol")
  
  
  
  return (layers)
}

