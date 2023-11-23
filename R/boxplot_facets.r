



#________________________________________________________
# Specify individual y-axis limits and breaks for each facet.
#________________________________________________________

boxplot_facets <- function (params) {
  
  stopifnot(is_bare_environment(params))
  
  ggdata   <- params$.ggdata
  facet.by <- params$facet.by
  
  
  #________________________________________________________
  # Usual faceting configurations
  #________________________________________________________
  plot_facets(params)
  
  
  
  #________________________________________________________
  # Flipping the axes requires flipping the scales too.
  #________________________________________________________
  if (has_layer(params, 'facet') && isTRUE(params$flip))
    params$layers[['facet']] %<>% within({
      if (scales %in% c("free_x", "free_y"))
        scales %<>% switch(free_x = "free_y", free_y = "free_x")
    })
  
  
  
  #________________________________________________________
  # Re-stripe to account for facets missing x-categories.
  #________________________________________________________
  if (has_layer(params, 'stripe') && isTRUE(params$.free_x)) {
    
    attr(params$.ggdata, 'stripe') <- with(params, {
      plyr::ddply(.ggdata, ply_cols(facet.by), function (z) {
        data.frame(x = seq(2, length(unique(z[[.xcol]])), 2)) }) })
    
    set_layer(params, 'stripe', 'data' = ~ attr(., "stripe"), .overwrite = TRUE)
  }
  
  
  
  
  #________________________________________________________
  # Rotate the x-axis tick mark labels to avoid overlap
  #________________________________________________________
  xcol <- params$.xcol
  
  if (eq(xcol, '.all')) {
    
    if (isTRUE(params$flip)) {
      set_layer(
        params = params, 
        layer  = 'theme',
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank() )
      
    } else {
      set_layer(
        params = params, 
        layer  = 'theme',
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank() )
    }
    
  } else {
    
    if (eq(substr(xcol, 1, 1), ".")) {
      if (isTRUE(params$flip)) { set_layer(params, 'theme', axis.title.y = element_blank())
      } else                   { set_layer(params, 'theme', axis.title.x = element_blank()) }
    }
    
    
    ggdata     <- params$.ggdata
    facet.cols <- params$.plot_attrs$facet.ncol
    color.by   <- names(params$color.by)
    xlab.angle <- params$xlab.angle
    
    # x-axis Text Angle
    if (tolower(xlab.angle) == 'auto') {
      if (!is_null(xcol) && !isTRUE(params$flip)) {
        charCount <- sum(nchar(unique(as.character(ggdata[[xcol]]))), na.rm = TRUE)
        charCount <- charCount * facet.cols
        if (charCount > 40)
          xlab.angle <- 30
        remove("charCount")
      }
    }
    
    if (xlab.angle == 90 || tolower(xlab.angle) == "vertical") {
      set_layer(params, 'theme', axis.text.x = element_text(angle=-90, vjust=0.3, hjust=0))
      
    } else if (xlab.angle == 30 || tolower(xlab.angle) == "angled") {
      set_layer(params, 'theme', axis.text.x = element_text(angle=-30, vjust=1, hjust=0))
      
      # Ensure long x-axis labels don't get truncated at the figure's edge
      rpad <- strwidth(tail(levels(ggdata[[xcol]]), 1), units="inches")
      rpad <- rpad * 0.8660254 # sin((90-30) * pi / 180) / sin(90 * pi / 180)
      if (!is_null(color.by))
        rpad <- rpad - max(c(
          strwidth(color.by, units="inches", cex=1.2) + .20,
          strwidth(levels(ggdata[[color.by]]), units="inches") + .52 ))
      rpad <- max(.1, signif(rpad, digits = 3))
      set_layer(
        params = params, 
        layer  = 'theme', 
        'plot.margin' = as.cmd(unit(x=c(.1, rpad, .1, .1), units='inches'), list(rpad=rpad)) )
    }
    
    remove("ggdata", "facet.cols", "color.by", "xlab.angle")
  }
  
  
  
  # To enable %>% chaining
  return (invisible(params))
}
