

boxplot_build <- function (params, data_func, layers_func) {
  
  
  #________________________________________________________
  # Set up the x-axis variables 
  #________________________________________________________
  if (identical(params[['x']], ".")) {
    params[['xval.by']] <- NA
  } else {
    params[['xval.by']] <- as.vector(params[['x']])
    
    for (i in params[['xval.by']]) {
      
      if (identical(substr(i, 1, 1), ".")) next
      i <- sub("^[!=]=", "", i)
      
      vals <- params[['biom']][['metadata']][[i]]
      if (is.character(vals)) {
        params[['biom']][['metadata']][[i]] <- as.factor(vals)
      } else if (!is.factor(vals)) {
        stop("Non-categorical '", i, "' cannot be a boxplot x-axis value.")
      }
      remove("vals")
    }
  }
  
  
  
  #________________________________________________________
  # Convert user's `layers` spec to layer names.
  #________________________________________________________
  layer_names <- local({
    layerlist <- c( 'b' = "box", 'v' = "violin",   'e' = "errorbar", 
                    'r' = "bar", 's' = "strip",    'l' = "linerange",
                    'd' = "dot", 'c' = "crossbar", 'p' = "pointrange" )
    
    layer_match(params[['layers']], choices = layerlist, default = "rls") %>%
      c('ggplot', ., 'xaxis', 'yaxis', 'labs', 'theme', 'theme_bw')
  })
  
  
  
  #________________________________________________________
  # Ignore shapes/etc without applicable layers
  #________________________________________________________
  if (!any(c('dot', 'strip')         %in% layer_names)) params[['shape.by']]   <- NULL
  if (!any(c('box', 'bar', 'violin') %in% layer_names)) params[['pattern.by']] <- NULL
  
  
  #________________________________________________________
  # Create the plot's data frame
  #________________________________________________________
  ggdata <- data_func(params)        # Plot-specific
  ggdata <- metadata_filters(ggdata) # General
  remove("params")
    
  
  
  #________________________________________________________
  # Avoid warning messages that come from this:
  # ggplot(data.frame(x="A", y=c(1, 1 - 1e-15))) 
  # + ggbeeswarm::geom_quasirandom(aes(x, y))
  #________________________________________________________
  ggdata[['.y']] %<>% round(12)
  
  
  #________________________________________________________
  # Initialize the `layers` object
  #________________________________________________________
  layers <- list()
  
  # Move attributes from `ggdata` to `layers`
  for (i in c('xcol', 'ycol', 'params', 'stats')) {
    attr(layers, i) <- attr(ggdata, i, exact = TRUE)
    attr(ggdata, i) <- NULL
  }
  attr(layers, 'data')  <- ggdata
  attr(layers, 'xmode') <- "factor"
  
  initLayer(layer_names)
  layers <- layers_func(layers)     # Plot-specific
  layers <- metadata_layers(layers) # General
  
  remove("i", "ggdata", "layer_names")
  
  
  
  
  
  #________________________________________________________
  # Help multiple plots overlay well together
  #________________________________________________________
  
  dodged    <- hasLayer(c("color", "shape", "pattern")) %>% any()
  patterned <- hasLayer("pattern")
  dotted    <- hasLayer(c("dot", "strip")) %>% any()
  
  dodge  <- as.cmd(position_dodge(width = 0.8))
  jitter <- as.cmd(position_jitter(width = 0.25, height = 0, seed = 0))
  jdodge <- as.cmd(position_jitterdodge(
    dodge.width = 0.8, jitter.width = 0.05, jitter.height = 0, seed = 0) )
  
  
  
  #________________________________________________________
  # Violin arguments
  #________________________________________________________
  layer <- "violin"
  if (hasLayer()) {
    setLayer(color = "black")
    
    if (dodged)
      setLayer(position = dodge)
    
    if (hasLayer("bar"))
      setLayer(alpha = 0.5)
    
    if (patterned) {
      setLayer(pattern="magick", fill="white", color="black", pattern_res=120)
      
      if (!hasLayer("color"))
        setLayer(pattern_fill="black")
    }
  }
  
  
  #________________________________________________________
  # Bar plot arguments
  #________________________________________________________
  layer <- "bar"
  if (hasLayer()) {
    setLayer(stat="summary", fun="mean", alpha=0.6)
    
    if (dodged)
      setLayer(position=dodge, width=0.7)
    
    if (patterned) {
      setLayer(pattern="magick", fill="white", color="black", pattern_res=120)
      
      if (!hasLayer("color"))
        setLayer(pattern_fill = "black")
    }
    
    # Bar charts need extra help on non-linear axes
    if (!is.null(layers[['yaxis']][['trans']])) {
      trans <- layers[['yaxis']][['trans']]
      fun <- if (trans == 'sqrt')  { function (y) sqrt(mean(y * y))
      } else if (trans == 'log1p') { function (y) log1p(mean(exp(y) - 1))
      } else { stop('Bar charts can only be re-scaled using sqrt or log1p.') }
      setLayer(fun=fun)
      remove("fun", "trans")
    }
  }
  
  
  #________________________________________________________
  # Boxplot arguments
  #________________________________________________________
  layer <- "box"
  if (hasLayer()) {
    setLayer(color="black", width=0.7)
    
    if (hasLayer("violin")) {
      setLayer(fill="white", width=0.1)
      
    } else if (patterned)     {
      setLayer(pattern="magick", fill="white", color="black", pattern_res=120)
      
      if (!hasLayer("color"))
        setLayer(pattern_fill="black")
      
    } else {
      setLayer(alpha=0.6)
    }
    
    if (dodged)
      setLayer(position=dodge, outlier.shape=NA)
  }
  
  
  #________________________________________________________
  # Adjust default point size based on number of points.
  #________________________________________________________
  ptsize <- local({
    n <- nrow(attr(layers, 'data'))
    round(log(1002 - min(1000, max(10, n))) / 6, digits = 2)
  })
  
  
  #________________________________________________________
  # Dotplot arguments for ggbeeswarm::geom_beeswarm
  #________________________________________________________
  layer <- "dot"
  if (hasLayer()) {
    
    setLayer(cex = 0.5, size = ptsize)
    
    if (packageVersion("ggbeeswarm") >= "0.7.0.9000") {
      setLayer(corral = "random", corral.width = 0.9)
    } else {
      setLayer(groupOnX = TRUE)
    }
    
    
    if (dodged)
      setLayer(dodge.width = 0.8)
    
    if (any(hasLayer(c("violin", "box", "bar"))))
      setLayer(color = "black", fill = "black", 'mapping|group' = ".group")
    
    if (any(hasLayer(c("crossbar", "errorbar", "linerange", "pointrange"))))
      setLayer(stroke = 0, alpha = 0.4)
  }
  
  
  #________________________________________________________
  # Stripchart arguments for ggbeeswarm::geom_quasirandom
  #________________________________________________________
  layer <- "strip"
  if (hasLayer()) {
    
    setLayer(size = ptsize)
    
    if (packageVersion("ggbeeswarm") < "0.7.0.9000")
      setLayer(groupOnX = TRUE)
    
    if (dodged)
      setLayer(dodge.width = 0.8)
    
    if (any(hasLayer(c("violin", "box", "bar"))))
      setLayer(color = "black", 'mapping|group' = ".group")
    
    if (any(hasLayer(c("crossbar", "errorbar", "linerange", "pointrange"))))
      setLayer(stroke = 0, alpha = 0.4)
  }
  
  
  #________________________________________________________
  # Vertical line arguments
  #________________________________________________________
  if (any(hasLayer(c("crossbar", "errorbar", "linerange", "pointrange")))) {
    
    ci <- attr(layers, 'params', exact = TRUE)[['ci']]
    
    if (is.numeric(ci)) {
      
      if (ci < 75 || ci > 100)
        stop("`ci` must be between 75 and 100. Value given = ", ci)
      
      vlineFn <- function (vals) {
        tt <- try(t.test(vals, conf.level = ci / 100), silent = TRUE)
        if (!is(tt, "htest")) tt <- list(estimate = NA, conf.int = c(NA, NA))
        data.frame(
          .y    = ifelse(isTRUE(unname(tt$estimate) >= 0), unname(tt$estimate), 0),
          .ymin = ifelse(isTRUE(tt$conf.int[1]      >= 0), tt$conf.int[1],      0),
          .ymax = ifelse(isTRUE(tt$conf.int[2]      >= 0), tt$conf.int[2],      0)
        )
      }
    } else if (ci == "range") {
      vlineFn <- function (vals) {
        data.frame(.y = mean(vals), .ymin = min(vals), .ymax = max(vals))
      }
    } else if (ci == "mad") {
      vlineFn <- function (vals) {
        med <- median(vals); dev <- mad(vals, med)
        data.frame(.y = med, .ymin = med - dev, .ymax = med + dev)
      }
    } else if (ci == "sd") {
      vlineFn <- function (vals) {
        avg <- mean(vals); dev <- sd(vals)
        data.frame(.y = avg, .ymin = avg - dev, .ymax = avg + dev)
      }
    } else if (ci == "se") {
      vlineFn <- function (vals) {
        avg <- mean(vals); dev <- sqrt(var(vals)/length(vals))
        data.frame(.y = avg, .ymin = avg - dev, .ymax = avg + dev)
      }
    } else {
      stop("`ci` must be either NULL, a number between 75 and 100, 'mad', 'sd', 'se', or 'range'. Value given = ", ci)
    }
    
    
    #________________________________________________________
    # Append rows to ggdata for the vline function
    #________________________________________________________
    ggdata <- attr(layers, 'data', exact = TRUE)
    
    group_cols <- as.vector(unlist(c(
      attr(layers, 'xcol', exact = TRUE),
      ".group",
      attr(layers, 'params', exact = TRUE)[c('color.by', 'shape.by', 'pattern.by', 'facet.by')] )))
    
    plyby    <- ply_cols(intersect(group_cols, colnames(ggdata)))
    vline_df <- plyr::ddply(ggdata, plyby, function (v) { vlineFn(v[[".y"]]) })
    
    if (nrow(vline_df) > 0) {
      attr(ggdata, "vline") <- vline_df
      attr(layers, 'data')  <- ggdata
    }
    
    remove("ggdata", "plyby", "group_cols", "vline_df")
    
    
    #________________________________________________________
    # Additional arguments for the vline layer. `data` and
    # `mapping` parameters are handled by plot_build().
    #________________________________________________________
    args <- list()
    
    if (dodged)                            args[['position']] <- dodge
    if (patterned)                         args[['color']]    <- "black"
    if (any(hasLayer(c("violin", "box")))) args[['color']]    <- "black"
    
    if (length(args) > 0)
      for (layer in intersect(c("errorbar", "linerange", "pointrange"), names(layers)))
        setLayer(layer, args)
    
    if (hasLayer("errorbar")) setLayer("errorbar", list('width' = 0.25))
    if (hasLayer("crossbar")) setLayer("crossbar", list('width' = 0.25))
    
    if (hasLayer("crossbar")) {
      if (any(hasLayer(c("violin", "box", "bar")))) { args %<>% c(fill="transparent")
      } else if (patterned)                         { args %<>% c(pattern="magick", pattern_res=120)
      } else                                        { args %<>% c(fill="transparent") }
      setLayer("crossbar", args)
    }
    
    remove("ci", "args")
  }
  
  
  #________________________________________________________
  # Theme arguments
  #________________________________________________________
  layer <- "theme"
  setLayer(text=element_text('size' = 14))
  
  
  
  #________________________________________________________
  # aes() defaults
  #________________________________________________________
  specs <- list(
    "box"        = c('x', 'y', 'group', 'color', 'fill', 'pattern_type', 'pattern_fill'),
    "violin"     = c('x', 'y', 'group', 'color', 'fill', 'pattern_type', 'pattern_fill'),
    "bar"        = c('x', 'y', 'group', 'color', 'fill', 'pattern_type', 'pattern_fill'),
    "dot"        = c('x', 'y', 'group', 'color',         'shape'),
    "strip"      = c('x', 'y', 'group', 'color',         'shape'),
    "crossbar"   = c('x', 'y', 'group', 'color', 'fill', 'ymin', 'ymax'),
    "pointrange" = c('x', 'y', 'group', 'color', 'fill', 'ymin', 'ymax'),
    "linerange"  = c('x', 'y', 'group', 'color',         'ymin', 'ymax'),
    "errorbar"   = c('x', 'y', 'group', 'color',         'ymin', 'ymax') )
  
  ggdata <- attr(layers, 'data',   exact = TRUE)
  params <- attr(layers, 'params', exact = TRUE)
  
  xcol <- attr(layers, 'xcol', exact = TRUE)
  args <- c(list('x' = xcol), .qw(y, xmin, xmax, ymin, ymax, xend, yend, label))
  if (hasName(ggdata, ".group"))             args[['group']]        <- ".group"
  if (hasLayer("color"))                     args[['color']]        <- params[['color.by']]
  if (hasLayer("color"))                     args[['fill']]         <- params[['color.by']]
  if (hasLayer("shape"))                     args[['shape']]        <- params[['shape.by']]
  if (hasLayer("pattern"))                   args[['pattern_type']] <- params[['pattern.by']]
  if (all(hasLayer(c("color", "pattern"))))  args[['pattern_fill']] <- params[['color.by']]
  
  for (layer in intersect(names(layers), names(specs))) {
    layerArgs <- args[intersect(specs[[layer]], names(args))]
    layerArgs <- args[setdiff(names(layerArgs), names(layers[[layer]]))]
    
    if (length(layerArgs) == 0) next
    
    names(layerArgs) <- paste0("mapping|", names(layerArgs))
    setLayer(layer=layer, layerArgs)
  }
  
  remove(list = c("specs", "ggdata", "xcol", "args", "layer", "layerArgs", "params") %>% intersect(ls()))
  
  
  p <- layers %>%
    boxplot_facets() %>%
    boxplot_stats() %>%
    plot_build()
  
  
  return (p)
}
