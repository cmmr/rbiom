

init_boxplot_layers <- function (params = parent.frame()) {
  
  stopifnot(is_bare_environment(params))
  
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  with(params, {
     
    stopifnot(is_tibble(.ggdata))
    
    validate_bool('flip')
    validate_bool('caption')
    validate_bool('stripe',   default = flip)
    validate_bool('outliers', null_ok = TRUE)
    
    validate_var_range('p.label', range = c(0, 1), max = Inf)
    validate_var_range('level',   range = c(0, 1))
    
    validate_var_choices('ci',         choices = c('ci', 'range', 'sd', 'se', 'mad'))
    validate_var_choices('p.adj',      choices = p.adjust.methods)
    validate_var_choices('xlab.angle', choices = c('auto', '0', '30', '90'))
  })
  
  
  
  #________________________________________________________
  # Clean up the y-axis values.
  #________________________________________________________
  with(params, {
    
    .ycol   <- attr(.ggdata, 'response', exact = TRUE)
    stopifnot(hasName(.ggdata, .ycol))
    .ggdata <- .ggdata[is.finite(.ggdata[[.ycol]]),]
    
    # Avoid warning messages that come from this:
    # ggplot(data.frame(x="A", y=c(1, 1 - 1e-15))) 
    # + ggbeeswarm::geom_quasirandom(aes(x, y))
    #________________________________________________________
    .ggdata[[.ycol]] %<>% round(12)
  })
  
  
  
  #________________________________________________________
  # Establish how to group samples along the x-axis/facets.
  #________________________________________________________
  with(params, {
    
    if (is.null(x) || x %in% facet.by) {
      x <- ".all"
      .ggdata[[x]] <- factor("all")
    }
    
    
    .group_by <- intersect(colnames(.ggdata), unique(c(
      x, facet.by, names(color.by), names(shape.by), names(pattern.by) )))
    
    .ggdata[['.group']] <- .ggdata[,.group_by] %>%
      apply(1L, paste, collapse="|") %>%
      unname()
  })
  
  
  
  #________________________________________________________
  # Convert user's `layers` spec to layers environment.
  #________________________________________________________
  
  params$.xcol  <- params$x
  params$.xmode <- "factor"
  
  init_layers(
    params  = params,
    choices = c( 'x' = "box", 'v' = "violin",   'e' = "errorbar", 
                 'b' = "bar", 's' = "strip",    'l' = "linerange",
                 'd' = "dot", 'c' = "crossbar", 'p' = "pointrange" ) )
  
  layers <- params$layers
  
  
  
  #________________________________________________________
  # Help multiple plots overlay well together
  #________________________________________________________
  
  dodged    <- has_layer(params, c('color', 'shape', 'pattern')) %>% any()
  patterned <- has_layer(params, 'pattern')
  dotted    <- has_layer(params, c('dot', 'strip')) %>% any()
  
  dodge  <- as.cmd(position_dodge(width = 0.8))
  jitter <- as.cmd(position_jitter(width = 0.25, height = 0, seed = 0))
  jdodge <- as.cmd(position_jitterdodge(
    dodge.width = 0.8, jitter.width = 0.05, jitter.height = 0, seed = 0) )
  
  
  
  #________________________________________________________
  # Violin arguments
  #________________________________________________________
  if (has_layer(params, 'violin')) {
    set_layer(params, 'violin', color = "black")
    
    if (dodged)
      set_layer(params, 'violin', position = dodge)
    
    if (has_layer(params, 'bar'))
      set_layer(params, 'violin', alpha = 0.5)
    
    if (patterned) {
      set_layer(params, 'violin', pattern="magick", fill="white", color="black", pattern_res=120)
      
      if (!has_layer(params, 'color'))
        set_layer(params, 'violin', pattern_fill="black")
    }
  }
  
  
  #________________________________________________________
  # Bar plot arguments
  #________________________________________________________
  if (has_layer(params, 'bar')) {
    set_layer(params, 'bar', stat="summary", alpha=0.6)
    
    if (dodged)
      set_layer(params, 'bar', position=dodge, width=0.7)
    
    if (patterned) {
      set_layer(params, 'bar', pattern="magick", fill="white", color="black", pattern_res=120)
      
      if (!has_layer(params, 'color'))
        set_layer(params, 'bar', pattern_fill = "black")
    }
    
    # Bar charts need extra help on non-linear axes
    if (is_null(layers[['yaxis']][['trans']])) {
      set_layer(params, 'bar', fun="mean")
      
    } else {
      trans <- layers[['yaxis']][['trans']]
      fun <- if (trans == 'sqrt')  { as.cmd(function (y) sqrt(mean(y * y)))
      } else if (trans == 'log1p') { as.cmd(function (y) log1p(mean(exp(y) - 1)))
      } else { stop('Bar charts can only be re-scaled using sqrt or log1p.') }
      set_layer(params, 'bar', fun=fun)
      remove("fun", "trans")
      
    }
  }
  
  
  #________________________________________________________
  # Boxplot arguments
  #________________________________________________________
  if (has_layer(params, 'box')) {
    set_layer(params, 'box', color="black", width=0.7)
    
    if (has_layer(params, 'violin')) {
      set_layer(params, 'box', fill="white", width=0.1)
      
    } else if (patterned)     {
      set_layer(params, 'box', pattern="magick", fill="white", color="black", pattern_res=120)
      
      if (!has_layer(params, 'color'))
        set_layer(params, 'box', pattern_fill="black")
      
    } else {
      set_layer(params, 'box', alpha=0.6)
    }
    
    if (dodged)
      set_layer(params, 'box', position=dodge) #, outlier.shape=NA)
    
    
    
    # Show/hide/customize boxplot outliers.
    #________________________________________________________
    with(params, {
      
      if (is.null(outliers))
        outliers <- !any(hasName(layers, c('dot', 'strip')))
      
      if (isFALSE(outliers)) {
        layers$box$outlier.shape <- NA
      } else {
        layers$box$outlier.size  <- params$pt.size
        layers$box$outlier.alpha <- params$pt.alpha
      }
    })
    
  }
  
  
  #________________________________________________________
  # Adjust default point size based on number of points.
  #________________________________________________________
  ptsize <- local({
    n <- nrow(params$.ggdata)
    round(log(1002 - min(1000, max(10, n))) / 6, digits = 2)
  })
  
  
  #________________________________________________________
  # Dotplot arguments for ggbeeswarm::geom_beeswarm
  #________________________________________________________
  if (has_layer(params, 'dot')) {
    
    set_layer(params, 'dot', cex = 0.5, size = ptsize, method="center")
    
    if (packageVersion("ggbeeswarm") >= "0.7.0.9000") {
      set_layer(params, 'dot', corral = "random", corral.width = 0.9)
    } else {
      set_layer(params, 'dot', groupOnX = TRUE)
    }
    
    
    if (dodged)
      set_layer(params, 'dot', dodge.width = 0.8)
    
    if (any(has_layer(params, c("violin", "box", "bar"))))
      set_layer(params, 'dot', color = "black", fill = "black", 'mapping|group' = ".group")
    
    if (any(has_layer(params, c("crossbar", "errorbar", "linerange", "pointrange"))))
      set_layer(params, 'dot', stroke = 0, alpha = 0.4)
  }
  
  
  #________________________________________________________
  # Stripchart arguments for ggbeeswarm::geom_quasirandom
  #________________________________________________________
  if (has_layer(params, 'strip')) {
    
    set_layer(params, 'strip', size = ptsize)
    
    if (packageVersion("ggbeeswarm") < "0.7.0.9000")
      set_layer(params, 'strip', groupOnX = TRUE)
    
    if (dodged)
      set_layer(params, 'strip', dodge.width = 0.8)
    
    if (any(has_layer(params, c("violin", "box", "bar"))))
      set_layer(params, 'strip', color = "black", 'mapping|group' = ".group")
    
    if (any(has_layer(params, c("crossbar", "errorbar", "linerange", "pointrange")))) {
      
      set_layer(params, 'strip', alpha = 0.4)
      
      if (is_null(params$shape.by))
        set_layer(params, 'strip', stroke = 0)
    }
  }
  
  
  #________________________________________________________
  # Vertical line arguments
  #________________________________________________________
  if (any(has_layer(params, c("crossbar", "errorbar", "linerange", "pointrange")))) {
    
    ci <- params$ci
    
    if (ci == "ci") {
      
      params$.vlineFun <- function (vals) {
        level <- params$level
        tt    <- try(t.test(vals, conf.level = level), silent = TRUE)
        if (!is(tt, "htest")) tt <- list(estimate = NA, conf.int = c(NA, NA))
        data.frame(
          .y    = ifelse(isTRUE(unname(tt$estimate) >= 0), unname(tt$estimate), 0),
          .ymin = ifelse(isTRUE(tt$conf.int[1]      >= 0), tt$conf.int[1],      0),
          .ymax = ifelse(isTRUE(tt$conf.int[2]      >= 0), tt$conf.int[2],      0)
        )
      }
      
    } else if (ci == "range") {
      params$.vlineFun <- function (vals) {
        data.frame(.y = mean(vals), .ymin = min(vals), .ymax = max(vals))
      }
    } else if (ci == "mad") {
      params$.vlineFun <- function (vals) {
        med <- median(vals); dev <- mad(vals, med)
        data.frame(.y = med, .ymin = med - dev, .ymax = med + dev)
      }
    } else if (ci == "sd") {
      params$.vlineFun <- function (vals) {
        avg <- mean(vals); dev <- sd(vals)
        data.frame(.y = avg, .ymin = avg - dev, .ymax = avg + dev)
      }
    } else if (ci == "se") {
      params$.vlineFun <- function (vals) {
        avg <- mean(vals); dev <- sqrt(var(vals)/length(vals))
        data.frame(.y = avg, .ymin = avg - dev, .ymax = avg + dev)
      }
    }
    
    
    #________________________________________________________
    # Generate the vline function's data frame.
    #________________________________________________________
    with(params, {
      
      attr(.ggdata, "vline") <- plyr::ddply(
        .data      = .ggdata, 
        .variables = ply_cols(intersect(colnames(.ggdata), c('.group', .group_by))), 
        .fun       = function (v) { .vlineFun(v[[.ycol]]) })
      
      remove('.vlineFun')
    })
    
    for (layer in intersect(c("crossbar", "errorbar", "linerange", "pointrange"), env_names(layers)))
      set_layer(params, layer, 'mapping|y' = ".y")
    
    
    #________________________________________________________
    # Additional arguments for the vline layer. `data` and
    # `mapping` parameters are handled by plot_build().
    #________________________________________________________
    args <- list()
    
    if (dodged)                            args[['position']] <- dodge
    if (patterned)                         args[['color']]    <- "black"
    if (any(has_layer(params, c("violin", "box")))) args[['color']]    <- "black"
    
    if (length(args) > 0)
      for (layer in intersect(c("errorbar", "linerange", "pointrange"), env_names(layers)))
        set_layer(params, layer, args)
    
    if (has_layer(params, 'errorbar')) set_layer(params, 'errorbar', list('width' = 0.25))
    if (has_layer(params, 'crossbar')) set_layer(params, 'crossbar', list('width' = 0.25))
    
    if (has_layer(params, 'crossbar')) {
      if (any(has_layer(params, c("violin", "box", "bar")))) { args %<>% c(fill="transparent")
      } else if (patterned)                         { args %<>% c(pattern="magick", pattern_res=120)
      } else                                        { args %<>% c(fill="transparent") }
      set_layer(params, 'crossbar', args)
    }
    
    remove("ci", "args")
  }
  
  
  #________________________________________________________
  # Shade background of x-axis positions
  #________________________________________________________
  if (isTRUE(params$stripe) && nlevels(params$.ggdata[[params$.xcol]]) > 1) {
    
    stripe_x <- seq(2, nlevels(params$.ggdata[[params$.xcol]]), 2)
    
    set_layer(
      params = params, 
      layer  = 'stripe',
      'data'    = as.cmd(data.frame(x = stripe_x), list(stripe_x = stripe_x)),
      'mapping' = aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -Inf, ymax = Inf),
      'fill'    = 'black', color = NA, alpha = 0.05 )
    
    remove("stripe_x")
  }
  
  
  
  
  
  
  #________________________________________________________
  # Theme arguments
  #________________________________________________________
  set_layer(params, 'theme', text=element_text('size' = 14))
  
  
  
  #________________________________________________________
  # aes() defaults
  #________________________________________________________
  specs <- list(
    "box"        = c('x', 'y', 'group', 'color', 'fill', 'pattern_type', 'pattern_fill'),
    "violin"     = c('x', 'y', 'group', 'color', 'fill', 'pattern_type', 'pattern_fill'),
    "bar"        = c('x', 'y', 'group', 'color', 'fill', 'pattern_type', 'pattern_fill'),
    "dot"        = c('x', 'y', 'group', 'color',         'shape'),
    "strip"      = c('x', 'y', 'group', 'color',         'shape'),
    "pointrange" = c('x', 'y', 'group', 'color', 'fill', 'shape', 'ymin', 'ymax'),
    "crossbar"   = c('x', 'y', 'group', 'color', 'fill',          'ymin', 'ymax'),
    "linerange"  = c('x', 'y', 'group', 'color',                  'ymin', 'ymax'),
    "errorbar"   = c('x', 'y', 'group', 'color',                  'ymin', 'ymax') )
  
  
  args <- c(list('x' = params$.xcol, 'y' = params$.ycol), .qw(xmin, xmax, ymin, ymax, xend, yend, label))
  if (hasName(params$.ggdata, '.group'))             args[['group']]        <- ".group"
  if (has_layer(params, 'color'))                    args[['color']]        <- names(params$color.by)
  if (has_layer(params, 'color'))                    args[['fill']]         <- names(params$color.by)
  if (has_layer(params, 'shape'))                    args[['shape']]        <- names(params$shape.by)
  if (has_layer(params, 'pattern'))                  args[['pattern_type']] <- names(params$pattern.by)
  if (all(has_layer(params, c('color', 'pattern')))) args[['pattern_fill']] <- names(params$color.by)
  if (!is_null(layers[['color']][['values']])) set_layer(params, 'fill', values = layers[['color']][['values']])
  
  for (layer in intersect(names(layers), names(specs))) {
    layerArgs <- args[intersect(specs[[layer]], names(args))]
    layerArgs <- args[setdiff(names(layerArgs), names(layers[[layer]]))]
    
    if (length(layerArgs) == 0) next
    
    names(layerArgs) <- paste0("mapping|", names(layerArgs))
    set_layer(params, layer, layerArgs)
  }
  
  
  
  return (invisible(params))
}
