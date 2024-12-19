


boxplot_build <- function (params) {
  
  stopifnot(is_bare_environment(params))
  
  
  with(params, {
    
    
    #________________________________________________________
    # Validate and pre-process user's arguments.
    #________________________________________________________
    if (inherits(df, 'rbiom')) df <- df$metadata
    if (!inherits(df, 'data.frame'))
      cli_abort("`df` must be a data.frame or rbiom object, not {.type {df}}.")
    
    # Enforce unnamed vectors.
    for (.i in colnames(df)) df[[.i]] %<>% unname()
    remove(".i")
    
    validate_df_field('x',        col_type = "cat", null_ok = TRUE)
    validate_df_field('y',        col_type = "num")
    validate_df_field('stat.by',  col_type = "cat", null_ok = TRUE)
    validate_df_field('facet.by', col_type = "cat", null_ok = TRUE, max = Inf)
    
    validate_var_choices('test',       choices = c('auto', 'none'))
    validate_var_choices('p.adj',      choices = p.adjust.methods)
    validate_var_choices('ci',         choices = c('ci', 'range', 'sd', 'se', 'mad'))
    validate_var_choices('xlab.angle', choices = c('auto', '0', '30', '90'))
    
    validate_var_range('level',   range = c(0.5, 1), n = 1)
    validate_var_range('p.label', range = c(0, 1))
    
    validate_bool('flip')
    validate_bool('stripe',   default = flip)
    validate_bool('outliers', null_ok = TRUE)
    validate_bool('caption')
    
    
    #________________________________________________________
    # Consistent naming for plot_* functions.
    #________________________________________________________
    .ggdata <- as_tibble(df)
    .xcol  <- x
    .ycol  <- y
    .xmode <- "factor"
    remove("df", "x", "y")
    
    
    #________________________________________________________
    # Check for uniqueness among named metadata fields.
    #________________________________________________________
    if (any(.xcol %in% facet.by)) {
      cli_warn("{.val {(.xcol)}} is assigned to `x`; removing from `facet.by`.", .frequency = "always")
      facet.by %<>% setdiff(.xcol)
    }
    if (any(stat.by %in% facet.by)) {
      cli_warn("{.val {stat.by}} is assigned to `stat.by`; removing from `facet.by`.", .frequency = "always")
      facet.by %<>% setdiff(stat.by)
    }
    
    
    
    #________________________________________________________
    # Clean up the y-axis values.
    #________________________________________________________
    .ggdata <- .ggdata[is.finite(.ggdata[[.ycol]]),]
    
    # Avoid warning messages that come from this:
    # ggplot(data.frame(x="A", y=c(1, 1 - 1e-15))) 
    # + ggbeeswarm::geom_quasirandom(aes(x, y))
    #________________________________________________________
    .ggdata[[.ycol]] %<>% round(12)
    
    
    #________________________________________________________
    # Set lower/upper bounds for confidence interval range.
    #________________________________________________________
    if (is.null(.dots$ci.min)) .dots$ci.min <- -Inf
    if (is.null(.dots$ci.max)) .dots$ci.max <-  Inf
    
  })
  
  
  
  #________________________________________________________
  # Convert user's `layers` spec to layers environment.
  #________________________________________________________
  
  init_layers(
    params  = params,
    choices = c( 'x' = "box", 'v' = "violin",   'e' = "errorbar", 
                 'b' = "bar", 's' = "strip",    'l' = "linerange",
                 'd' = "dot", 'c' = "crossbar", 'p' = "pointrange" ) )
  
  layers <- params$layers
  
  
  #________________________________________________________
  # All layers share `x`, and `y` metadata fields.
  #________________________________________________________
  set_layer(
    params = params, 
    layer  = 'ggplot', 
    'mapping|x' = if.null(params$.xcol, 0), 
    'mapping|y' = params$.ycol )
  
  
  
  
  
  #________________________________________________________
  # Help multiple layers overlay well together
  #________________________________________________________
  
  dodged    <- !eq(params$.xcol, params$stat.by)
  patterned <- !is.null(params$patterns)
  dotted    <- has_layer(params, c('dot', 'strip')) %>% any()
  
  dodge  <- as.cmd(position_dodge(width = 0.8))
  jitter <- as.cmd(position_jitter(width = 0.25, height = 0, seed = 0))
  jdodge <- as.cmd(position_jitterdodge(
    dodge.width = 0.8, jitter.width = 0.05, jitter.height = 0, seed = 0) )
  
  
  
  #________________________________________________________
  # Violin ====
  #________________________________________________________
  if (has_layer(params, 'violin')) {
    
    if (!patterned)
      set_layer(params, 'violin', color = "black")
    
    if (dodged)
      set_layer(params, 'violin', position = dodge)
    
    if (has_layer(params, 'bar'))
      set_layer(params, 'violin', alpha = 0.5)
  }
  
  
  
  
  
  #________________________________________________________
  # Bar ====
  #________________________________________________________
  if (has_layer(params, 'bar')) {
    set_layer(params, 'bar', stat="summary", alpha=0.6)
    
    if (dodged)
      set_layer(params, 'bar', position=dodge, width=0.7)
    
    
    # Bar charts need extra help on non-linear axes
    if (is_null(layers[['yaxis']][['transform']])) {
      set_layer(params, 'bar', fun="mean")
      
    } else {
      transform <- layers[['yaxis']][['transform']]
      fun <- if (transform == 'sqrt')  { as.cmd(function (y) sqrt(mean(y * y)))
      } else if (transform == 'log1p') { as.cmd(function (y) log1p(mean(exp(y) - 1)))
      } else { stop('Bar charts can only be re-scaled using sqrt or log1p.') }
      set_layer(params, 'bar', fun=fun)
      remove("fun", "transform")
      
    }
  }
  
  
  #________________________________________________________
  # Box ====
  #________________________________________________________
  if (has_layer(params, 'box')) {
    
    if (has_layer(params, 'violin')) {
      set_layer(params, 'box', fill="white", width=0.1)
      
    } else {
      set_layer(params, 'box', width=0.7, alpha=0.6)
    }
    
    if (dodged)
      set_layer(params, 'box', position=dodge) #, outlier.shape=NA)
    
    
    
    # Show/hide/customize boxplot outliers.
    #________________________________________________________
    if (is.null(params$outliers))
      params$outliers <- !any(hasName(params$layers, c('dot', 'strip')))
    
    if (isFALSE(params$outliers)) {
      params$layers$box$outlier.shape <- NA
    } else {
      params$layers$box$outlier.size  <- params$pt.size
      params$layers$box$outlier.alpha <- params$pt.alpha
    }
    
  }
  
  
  #________________________________________________________
  # Adjust default point size based on number of points.
  #________________________________________________________
  ptsize <- local({
    n <- nrow(params$.ggdata)
    round(log(1002 - min(1000, max(10, n))) / 6, digits = 2)
  })
  
  
  #________________________________________________________
  # Dot ====
  # ggbeeswarm::geom_beeswarm
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
    
    if (patterned)
      set_layer(params, 'dot', fill = "black")
    
    if (any(has_layer(params, c("violin", "box", "bar"))))
      set_layer(params, 'dot', color = "black", fill = "black")
    
    if (any(has_layer(params, c("crossbar", "errorbar", "linerange", "pointrange"))))
      set_layer(params, 'dot', alpha = 0.4)
  }
  
  
  #________________________________________________________
  # Strip ====
  # ggbeeswarm::geom_quasirandom
  #________________________________________________________
  if (has_layer(params, 'strip')) {
    
    set_layer(params, 'strip', size = ptsize)
    
    if (packageVersion("ggbeeswarm") < "0.7.0.9000")
      set_layer(params, 'strip', groupOnX = TRUE)
    
    if (dodged)
      set_layer(params, 'strip', dodge.width = 0.8)
    
    if (patterned)
      set_layer(params, 'strip', fill = "black")
    
    if (any(has_layer(params, c("violin", "box", "bar"))))
      set_layer(params, 'strip', color = "black", fill = "black")
    
    if (any(has_layer(params, c("crossbar", "errorbar", "linerange", "pointrange")))) {
      
      set_layer(params, 'strip', alpha = 0.4)
    }
  }
  
  
  #________________________________________________________
  # CI ====
  #________________________________________________________
  if (any(has_layer(params, c("crossbar", "errorbar", "linerange", "pointrange")))) {
    
    ci <- params$ci
    
    if (ci == "ci") {
      
      params$.vlineFun <- function (vals) {
        level <- params$level
        tt    <- try(t.test(vals, conf.level = level), silent = TRUE)
        if (!inherits(tt, "htest")) tt <- list(estimate = NA, conf.int = c(NA, NA))
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
        med <- median(vals); dev <- stats::mad(vals, med)
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
          .variables = ply_cols(c(.xcol, stat.by, facet.by)), 
          .fun       = function (v) { .vlineFun(v[[.ycol]]) }) %>%
        dplyr::rename(!!.ycol := .y)
      
      remove('.vlineFun')
    })
    
    for (layer in intersect(c("crossbar", "errorbar", "linerange", "pointrange"), env_names(layers)))
      set_layer(
        params = params, 
        layer  = layer, 
        'mapping|ymin' = ".ymin", 
        'mapping|ymax' = ".ymax" )
    
    
    #________________________________________________________
    # Additional arguments for the vline layer. `data` and
    # `mapping` parameters are handled by plot_build().
    #________________________________________________________
    args <- list()
    
    if (dodged)                                     args[['position']] <- dodge
    if (patterned)                                  args[['color']]    <- "black"
    if (any(has_layer(params, c("violin", "box")))) args[['color']]    <- "black"
    
    if (length(args) > 0)
      for (layer in intersect(c("errorbar", "linerange", "pointrange"), env_names(layers)))
        set_layer(params, layer, args)
    
    if (has_layer(params, 'errorbar')) set_layer(params, 'errorbar', list('width' = 0.25))
    if (has_layer(params, 'crossbar')) set_layer(params, 'crossbar', list('width' = 0.25))
    
    if (has_layer(params, 'crossbar')) {
      if (any(has_layer(params, c("violin", "box", "bar")))) { args %<>% c(fill="transparent")
      } else if (!patterned)                                 { args %<>% c(fill="transparent") }
      set_layer(params, 'crossbar', args)
    }
    
    remove("ci", "args")
  }
  
  
  #________________________________________________________
  # Shade background of x-axis positions
  #________________________________________________________
  if (isTRUE(params$stripe) && nlevels(params$.ggdata[[params$.xcol]]) > 1) {
    
    stripe_x <- seq(2, nlevels(params$.ggdata[[params$.xcol]]), 2)
    
    x <- NULL # only for CRAN check
    
    set_layer(
      params = params, 
      layer  = 'stripe',
      'data'    = as.cmd(data.frame(x = stripe_x), list(stripe_x = stripe_x)),
      'mapping' = aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -Inf, ymax = Inf, x = NULL, y = NULL),
      'fill'    = 'black', color = NA, alpha = 0.05 )
    
    remove("stripe_x", "x")
  }
  
  
  
  
  
  
  
  #________________________________________________________
  # Build the plot.
  #________________________________________________________
  fig <- params %>% 
    plot_facets() %>% 
    boxplot_stats() %>% 
    plot_build()
  
  
  return (fig)
}

