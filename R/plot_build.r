
#______________________________________________________________
# Create the plot and add each layer with its arguments.
# Also attach a human-readable version of the plot command.
#______________________________________________________________
plot_build <- function (params) {
  
  
  layers <- params$layers
  ggdata <- params$.ggdata
  
  stopifnot(has_layer(params, 'ggplot'))
  stopifnot(inherits(ggdata, 'data.frame'))
  
  
  #______________________________________________________________
  # Return a placeholder plot if no data is present
  #______________________________________________________________
  if (nrow(ggdata) == 0) {
    p <- ggplot() + 
      geom_text(mapping = aes(x=1,y=1), label="No Data to Display") + 
      theme_void()
    p$plot_env <- emptyenv()
    attr(p, 'facet.nrow')  <- 1
    attr(p, 'facet.ncol')  <- 1
    attr(p, 'facet.count') <- 1
    return (p)
  }
  
  
  
  
  #________________________________________________________
  # Add in coord_flip when flip=TRUE
  #________________________________________________________
  if (isTRUE(params$flip)) add_layer(params, 'flip')
  
  
  #________________________________________________________
  # Theme arguments
  #________________________________________________________
  set_layer(params, 'theme', text=element_text('size' = 14))
  
  
  #______________________________________________________________
  # Suppress vertical gridlines (horizontal if flipped).
  #______________________________________________________________
  if (isTRUE(params$flip)) {
    set_layer(params, 'theme', panel.grid.major.y = element_blank())
  } else {
    set_layer(params, 'theme', panel.grid.major.x = element_blank())
  }
  
  
  
  #______________________________________________________________
  # x-axis and y-axis transformations
  #______________________________________________________________
  
  for (axis in c('x', 'y')) {
    
    var   <- sprintf(".%scol", axis)
    layer <- sprintf("%saxis", axis)
    label <- params$layers[['labs']][[axis]] %||% params[[var]] %||% ''
    xy_transform <- params$layers[[layer]][['transform']]
    
    
    if (is.null(xy_transform)) next
    
    
    # Set axis label to, e.g., `.ylab <- "Abundance (log scale)"`
    #______________________________________________________________
    if (!has_layer(params, 'labs')) add_layer(params, 'labs')
    params$layers[['labs']][[axis]] <- sprintf("%s\n(%s scale)", label, xy_transform)
    
    
    
    # Change 10000 to 10k for log scales.
    #______________________________________________________________
    if (isTRUE(xy_transform %in% c("log", "log10", "log1p")))
      set_layer(params, layer, labels = label_number(scale_cut = cut_si("")))
    
    
    
    # Force sqrt scale to display zero tick mark.
    # Handle sqrt- and log1p- transforming of Inf and -Inf
    #______________________________________________________________
    if (eq(xy_transform, "log1p") && isTRUE(params$stripe)) {
      
      params$layers[[layer]][['transform']] <- trans_new(
        name      = paste0(axis, "_log1p"), 
        transform = function (y) { y[is.finite(y)] <- base::log1p(y[is.finite(y)]); return (y); }, 
        inverse   = as.cmd(base::expm1) )
    
    } else if (eq(xy_transform, "sqrt")) {
      
      if (isTRUE(params$stripe)) {
        
        params$layers[[layer]][['transform']] <- trans_new(
          name      = paste0(axis, "_sqrt"), 
          transform = function (y) { y[is.finite(y)] <- base::sqrt(y[is.finite(y)]); return (y); }, 
          inverse   = function (y) ifelse(y<0, 0, y^2) )
        
      } else {
        
        params$layers[[layer]][['transform']] <- trans_new(
          name      = paste0(axis, "_sqrt"), 
          transform = as.cmd(base::sqrt), 
          inverse   = function (y) ifelse(y<0, 0, y^2) )
      }  
      
    }
  }
  
  
  
  #________________________________________________________
  # aes() defaults
  #________________________________________________________
  specs <- list(
    "trend"      = c('color'),
    "name"       = c('color'),
    "residual"   = c('color'),
    "linerange"  = c('color'),
    "errorbar"   = c('color'),
    "ellipse"    = c('color'),
    "spider"     = c('color'),
    "name"       = c('color'),
    "confidence" = c('color', 'fill'),
    "violin"     = c('color', 'fill'),
    "bar"        = c('color', 'fill'),
    "crossbar"   = c('color', 'fill'),
    "box"        = c('color', 'fill'),
    "stack"      = c('color', 'fill'),
    "pointrange" = c('color', 'fill', 'shape'),
    "point"      = c('color',         'shape'),
    "dot"        = c('color', 'fill', 'shape'),
    "strip"      = c('color', 'fill', 'shape') )
  
  args <- list()
  
  colors <- if (!is.null(params$colors))   list(values   = params$colors)   else NULL
  shapes <- if (!is.null(params$shapes))   list(values   = params$shapes)   else NULL
  fills  <- if (!is.null(params$patterns)) list(patterns = params$patterns) else colors
  
  if (!is.null(colors)) args[['color']] <- if (hasName(params, 'color.by')) params$color.by else params$stat.by
  if (!is.null(shapes)) args[['shape']] <- if (hasName(params, 'shape.by')) params$shape.by else params$stat.by
  if (!is.null(fills))  args[['fill']]  <- if (hasName(params, 'color.by')) params$color.by else params$stat.by
  
  
  
  for (layer in intersect(names(layers), names(specs))) {
    
    layerArgs <- args[intersect(specs[[layer]], names(args))]
    layerArgs <- args[setdiff(names(layerArgs), names(layers[[layer]]))]
    
    
    #________________________________________________________
    # `stat.by` is set, by colors/shapes/patterns = FALSE.
    #________________________________________________________
    if (length(layerArgs) == 0) {
      
      if (!is.null(params$stat.by))
        set_layer(params, layer, 'mapping|group' = params$stat.by)
      
      next
    }
    
    
    #________________________________________________________
    # Initialize colors/shapes/pattern scales.
    #________________________________________________________
    if (hasName(layerArgs, 'color')) set_layer(params, 'color', colors)
    if (hasName(layerArgs, 'shape')) set_layer(params, 'shape', shapes)
    if (hasName(layerArgs, 'fill'))  set_layer(params, 'fill',  fills)
    
    
    names(layerArgs) <- paste0("mapping|", names(layerArgs))
    set_layer(params, layer, layerArgs)
  }
  
  
  
  
  
  #______________________________________________________________
  # Clean up args
  #______________________________________________________________
  
  for (layer in names(layers)) {
    
    args <- layers[[layer]]
    fun  <- attr(args, 'function', exact = TRUE)
    src  <- attr(args, 'src',      exact = TRUE)
    
    
    if (is_null(src)) { data_cols <- colnames(ggdata)
    } else            { data_cols <- colnames(attr(ggdata, src, exact = TRUE)) }
    
    
    # Create the aes object for `mapping`=
    #______________________________________________________________
    if (hasName(args, 'mapping') && !inherits(args[['mapping']], "uneval")) {
    
      aes_args <- args[['mapping']]
      aes_args <- aes_args[!sapply(aes_args, is.null)]
      aes_args <- lapply(aes_args, as.vector)
      
      for (arg in names(aes_args))
        if (is_string(aes_args[[arg]], data_cols)) 
          aes_args[[arg]] <- as.name(aes_args[[arg]])
      
      args[['mapping']] <- NULL
      if (isTRUE(length(aes_args) > 0))
        args[['mapping']] <- do.call(aes, aes_args)
      
      remove("aes_args")
    }
    
    
    # Don't specify an argument if it's already the default
    #______________________________________________________________
    defaults <- formals(fun)
    for (i in intersect(names(defaults), names(args)))
      if (eq(defaults[[i]], args[[i]]))
        args[[i]] <- NULL
    
    
    layers[[layer]] <- list(fun = fun, args = args)
    
    remove("args", "fun", "src", "defaults")
  }
  
  
  #______________________________________________________________
  # Standardize the list order: ggplot() first, theme() last, etc
  #______________________________________________________________
  layer_order <- c(
    'ggplot', 'confidence', 'stats_bg', 'stripe', 'density', 'violin', 
    'residual', 'point', 'trend', 'bar', 'box', 'spider', 'dot', 'ellipse', 'strip', 
    'name', 'crossbar', 'linerange', 'rect', 'errorbar', 'pointrange', 'mean', 
    'arrow', 'taxon', 'brackets', 'stats_vline', 'stats_label', 'stats_text', 'stack', 'hline', 'vline', 
    'labs', 'color', 'fill', 'shape', 'pattern', 'size', 'continuous_scale', 
    'scale_size', 'facet', # 'free_y', 
    'xaxis', 'yaxis', 'flip', 'theme_bw', 'theme' )
  layer_order <- c(
    intersect(layer_order, env_names(layers)),
    setdiff(env_names(layers), layer_order) )
  
  
  
  
  #______________________________________________________________
  # Assemble the ggplot object (p)
  #______________________________________________________________
  p    <- NULL
  cmds <- NULL
  
  for (layer in layer_order) {
    
    fun  <- layers[[layer]][['fun']]
    args <- layers[[layer]][['args']]
    fn   <- attr(fun, 'fn', exact = TRUE)
    
    
    # Skip theme() unless it has arguments
    #______________________________________________________________
    if (layer == 'theme' && length(args) == 0)
      next
    
    
    
    # Show ggplot() layer as "ggplot(data)", rest more verbosely
    #______________________________________________________________
    if (fn == "ggplot") {
      
      args[['data']] <- structure(ggdata, display = "data")
      p    <- do.call(fun, args) 
      cmds <- attr(p, 'display')
      
    } else {
      args[['.indent']] <- 4
      p    <- p + (q <- do.call(fun, args))
      cmds <- c(cmds, attr(q, 'display'))
    }
    
  }
  
  
  
  #______________________________________________________________
  # Append QQ plots with patchwork
  #______________________________________________________________
  if (isTRUE(params$check)) {
    
    p <- local({
      
      .fitted <- .residual <- .std.res <- NULL # for CRAN check only
      
      df <- with(
        data = attr(params$.ggdata, 'residual'),
        expr = tibble(
          .y        = get(params$.ycol), 
          .fitted   = .fitted, 
          .residual = .residual, 
          .std.res  = .std.res ))
      
      qq <- ggplot2::ggplot(df, aes(sample = .std.res)) + 
        ggplot2::geom_qq_line(alpha = 0.5, fullrange = TRUE) +
        ggplot2::geom_qq(size = 0.2, alpha = 0.5) +
        ggplot2::labs(
          title    = "Normal Q-Q",
          x        = "Theoretical Quantiles",
          y        = "Standardized Residual",
          subtitle = "Points should fall on the line." )
      
      rf <- ggplot2::ggplot(df, aes(x = .fitted, y = .residual)) + 
        ggplot2::geom_hline(yintercept = 0, alpha = 0.5, linetype = "dashed") +
        ggplot2::geom_point(size = 0.2, alpha = 0.5) +
        ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
        ggplot2::labs(
          title    = "Residuals vs Fitted",
          x        = "Fitted Value",
          y        = "Residual",
          subtitle = "Flat and horizontal trendline." )
      
      sl <- ggplot2::ggplot(df, aes(x = .fitted, y = sqrt(abs(.std.res)))) + 
        ggplot2::geom_point(size = 0.2, alpha = 0.5) +
        ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
        ggplot2::labs(
          title    = "Scale-Location",
          x        = "Fitted Value",
          y        = expression(sqrt(abs(` Standardized Residual `))),
          subtitle = "Flat and horizontal trendline." )
      
      diag <- lapply(list(qq, rf, sl), function (x) {
        x +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            plot.title      = ggplot2::element_text(face = "bold", size = 10), 
            plot.subtitle   = ggplot2::element_text(face = "italic", size = 9),
            axis.title      = ggplot2::element_text(size = 9),
            plot.background = ggplot2::element_rect(linewidth = 1, color = "gray") )
      })
      
      patchwork::wrap_plots(
        ncol    = 1,
        heights = c(3, 1),
        patchwork::free(patchwork::wrap_plots(p)), 
        patchwork::wrap_plots(diag, nrow = 1) )
    })
  }
  
  
  
  #______________________________________________________________
  # Attach code, facet info, etc
  #______________________________________________________________
  for (i in names(params$.plot_attrs))
    attr(p, i) <- params$.plot_attrs[[i]]
  
  attr(p, 'code') <- cmds %>% 
    paste(collapse=" +\n  ") %>% 
    add_class('rbiom_code')
  
  
  #________________________________________________________
  # Enable accessing attributes with `$`.
  #________________________________________________________
  p %<>% add_class('rbiom_plot')
  
  
  #______________________________________________________________
  # Prevent p from including all of this environment's variables.
  #______________________________________________________________
  p$plot_env <- emptyenv()
  
  
  return (p)
}
