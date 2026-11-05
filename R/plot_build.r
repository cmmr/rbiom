
#______________________________________________________________
# Create the plot and add each layer with its arguments.
# Also attach the code that rebuilds the plot from `data`.
#
# Inputs read from `params`:
#   layers        - the layer environment from init_layers().
#   .ggdata       - the data.frame given to ggplot(). Derived tables that
#                   some layers draw from (the statistics annotations, CI
#                   intervals, trendlines, ...) are attributes of it; see
#                   layer_src() for which layer reads which attribute.
#   .data_code    - optional character vector: code lines that turn the
#                   `data` produced by the table function into .ggdata
#                   (NA filtering etc.). Written in terms of `.data_var`.
#   .derived_code - optional named list: code that computes a derived
#                   table (named by its attribute) from `data`. Tables
#                   without an entry are written out as literal tibbles.
#   .data_var     - the name of the data variable in the generated code
#                   ("data", or "d1", "d2", ... for multi-panel plots).
#   .plot_var     - when set, the plot is assigned to this variable
#                   instead of being the code's final expression.
#   .plot_attrs   - attributes to attach to the plot (stats, facet.*).
#   check         - append model diagnostic panels (corrplots).
#
# The generated code is a record of the computation; it is never
# evaluated here. Every string that originates from the data (column
# names, factor levels, sample and taxa names) reaches the code only
# through coan() (backticked names), double_quote(), or deparse_vec().
#______________________________________________________________

DERIVED_VARS <- c(
  stat_labels   = "labels",
  stat_brackets = "brackets",
  stat_vline    = "stat_vline",
  vline         = "intervals",
  fit           = "fit",
  residual      = "residuals",
  spider        = "spider",
  taxa_coords   = "taxa_coords",
  ellipse       = "ellipse_data" )


plot_build <- function (params) {


  layers   <- params$layers
  ggdata   <- params$.ggdata
  data_var <- params$.data_var %||% "data"
  plot_var <- params$.plot_var
  panel    <- !is.null(params$.data_var)

  stopifnot(has_layer(params, 'ggplot'))
  stopifnot(inherits(ggdata, 'data.frame'))

  if (isTRUE(params$check)) plot_var <- plot_var %||% "p"


  #________________________________________________________
  # Provenance-tracked ggplot2 functions.
  #________________________________________________________
  .ggplot        <- P('ggplot2::ggplot')
  .geom_text     <- P('ggplot2::geom_text')
  .theme_void    <- P('ggplot2::theme_void')
  .element_blank <- P('ggplot2::element_blank')
  .element_text  <- P('ggplot2::element_text')
  .trans_new     <- P('scales::trans_new')
  .label_number  <- P('scales::label_number')
  .cut_si        <- P('scales::cut_si')


  #______________________________________________________________
  # Return a placeholder plot if no data is present
  #______________________________________________________________
  if (nrow(ggdata) == 0) {

    p <- .ggplot() +
      .geom_text(mapping = aes(x = 1, y = 1), label = "No Data to Display") +
      .theme_void()

    code <- "ggplot() + geom_text(mapping = aes(x = 1, y = 1), label = \"No Data to Display\") + theme_void()"
    if (!is.null(plot_var)) code <- sprintf("%s <- %s", plot_var, code)

    return (plot_finalize(p, params, code, ggdata))
  }


  #________________________________________________________
  # Add in coord_flip when flip=TRUE
  #________________________________________________________
  if (isTRUE(params$flip)) add_layer(params, 'flip')


  #________________________________________________________
  # Theme arguments
  #________________________________________________________
  set_layer(params, 'theme', text = .element_text('size' = 14))


  #______________________________________________________________
  # Suppress vertical gridlines (horizontal if flipped).
  #______________________________________________________________
  if (isTRUE(params$flip)) {
    set_layer(params, 'theme', panel.grid.major.y = .element_blank())
  } else {
    set_layer(params, 'theme', panel.grid.major.x = .element_blank())
  }


  #______________________________________________________________
  # Angle of the x-axis labels (categorical x only).
  #______________________________________________________________
  if (eq(params$.xmode, "factor") && !is.null(params$xlab.angle) && !isTRUE(params$flip)) {

    angle <- as.character(params$xlab.angle)

    if (eq(angle, "auto")) {
      lvls  <- if (is.null(params$.xcol)) "" else levels(as.factor(ggdata[[params$.xcol]]))
      width <- (sum(nchar(lvls)) + 2 * length(lvls)) * (params$.plot_attrs$facet.ncol %||% 1)
      angle <- if (width <= 50) "0" else if (width <= 100) "30" else "90"
    }

    if (eq(angle, "30")) {
      set_layer(params, 'theme', axis.text.x = .element_text(angle = -30, vjust = 1, hjust = 0))
    } else if (eq(angle, "90")) {
      set_layer(params, 'theme', axis.text.x = .element_text(angle = -90, vjust = 0.3, hjust = 0))
    }

    remove("angle")
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
    if (is_scalar_character(xy_transform)) {
      if (!has_layer(params, 'labs')) add_layer(params, 'labs')
      params$layers[['labs']][[axis]] <- sprintf("%s\n(%s scale)", label, xy_transform)
    }


    # Change 10000 to 10k for log scales.
    #______________________________________________________________
    if (isTRUE(xy_transform %in% c("log", "log10", "log1p")))
      set_layer(params, layer, labels = .label_number(scale_cut = .cut_si("")))


    # Force sqrt scale to display zero tick mark.
    # Handle sqrt- and log1p- transforming of Inf and -Inf (from the
    # stripe/background rectangles), which the stock transforms can't.
    #______________________________________________________________
    if (eq(xy_transform, "log1p") && isTRUE(params$stripe)) {

      params$layers[[layer]][['transform']] <- .trans_new(
        name      = paste0(axis, "_log1p"),
        transform = aa(function (y) { y[is.finite(y)] <- log1p(y[is.finite(y)]); y }, display = "function (y) { y[is.finite(y)] <- log1p(y[is.finite(y)]); y }"),
        inverse   = aa(function (y) expm1(y), display = "function (y) expm1(y)") )

    } else if (eq(xy_transform, "sqrt")) {

      if (isTRUE(params$stripe)) {

        params$layers[[layer]][['transform']] <- .trans_new(
          name      = paste0(axis, "_sqrt"),
          transform = aa(function (y) { y[is.finite(y)] <- sqrt(y[is.finite(y)]); y }, display = "function (y) { y[is.finite(y)] <- sqrt(y[is.finite(y)]); y }"),
          inverse   = aa(function (y) ifelse(y < 0, 0, y^2), display = "function (y) ifelse(y < 0, 0, y^2)") )

      } else {

        params$layers[[layer]][['transform']] <- .trans_new(
          name      = paste0(axis, "_sqrt"),
          transform = aa(function (y) sqrt(y), display = "function (y) sqrt(y)"),
          inverse   = aa(function (y) ifelse(y < 0, 0, y^2), display = "function (y) ifelse(y < 0, 0, y^2)") )
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




  #________________________________________________________
  # Every `...` argument must have reached a layer.
  #________________________________________________________
  check_unused_dots(params)


  #______________________________________________________________
  # Clean up args; collect the derived tables that layers draw from.
  #______________________________________________________________
  derived <- list()

  if (isTRUE(params$check) && !is.null(attr(ggdata, 'residual', exact = TRUE)))
    derived[['residual']] <- attr(ggdata, 'residual', exact = TRUE)

  for (layer in env_names(layers)) {

    args <- layers[[layer]]
    fun  <- attr(args, 'function', exact = TRUE)
    src  <- attr(args, 'src',      exact = TRUE)


    # Layers that draw a derived table get it as their `data`.
    #______________________________________________________________
    if (is_null(src)) {
      data_cols <- colnames(ggdata)

    } else {

      df <- attr(ggdata, src, exact = TRUE)

      if (is_empty_df(df)) {
        del_layer(params, layer)
        next
      }

      derived[[src]] <- df
      args[['data']] <- structure(df, display = DERIVED_VARS[[src]])
      data_cols      <- colnames(df)
    }


    # Create the aes object for `mapping`=
    #______________________________________________________________
    if (hasName(args, 'mapping') && !inherits(args[['mapping']], "uneval")) {

      aes_args <- args[['mapping']]
      aes_args <- aes_args[!sapply(aes_args, is.null)]
      aes_args <- lapply(aes_args, as.vector)
      aes_args <- aes_args[aes_order(names(aes_args))]

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
    defaults <- attr(fun, 'formals', exact = TRUE) %||% formals(fun)
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
    'scale_size', 'facet',
    'xaxis', 'yaxis', 'flip', 'theme_bw', 'theme' )
  layer_order <- c(
    intersect(layer_order, env_names(layers)),
    setdiff(env_names(layers), layer_order) )




  #______________________________________________________________
  # Assemble the ggplot object (p)
  #______________________________________________________________
  plot_data <- if (panel) strip_prov(ggdata) else ggdata

  p    <- NULL
  cmds <- NULL

  for (layer in layer_order) {

    fun  <- layers[[layer]][['fun']]
    args <- layers[[layer]][['args']]
    fn   <- attr(fun, 'display')


    # Skip theme() unless it has arguments
    #______________________________________________________________
    if (layer == 'theme' && length(args) == 0)
      next



    # Show ggplot() layer as "ggplot(data)", rest more verbosely
    #______________________________________________________________
    if (identical(fn, "ggplot")) {

      args[['data']] <- structure(plot_data, display = data_var)
      p    <- do.call(fun, args)
      cmds <- attr(p, 'display')

    } else {
      args[['.indent']] <- 4
      p    <- p + (q <- do.call(fun, args))
      cmds <- c(cmds, attr(q, 'display'))
    }

  }

  p$data <- plot_data



  #______________________________________________________________
  # The code: data prep, derived tables, then the ggplot expression.
  #______________________________________________________________
  code <- c()

  if (length(params$.data_code) > 0)
    code <- c(code, params$.data_code, "")

  for (src in names(derived)) {

    src_code <- params$.derived_code[[src]]

    # "" marks a table that the data-prep code already defines.
    if (identical(src_code, "")) next

    if (is.null(src_code))
      src_code <- derived_literal_code(DERIVED_VARS[[src]], derived[[src]])

    code <- c(code, src_code, "")
  }

  plot_code <- paste(cmds, collapse = " +\n  ")
  if (!is.null(plot_var)) plot_code <- sprintf("%s <- %s", plot_var, plot_code)
  code <- c(code, plot_code)



  #______________________________________________________________
  # Append model diagnostic panels with patchwork.
  #______________________________________________________________
  if (isTRUE(params$check)) {

    residuals <- derived[['residual']]

    if (is_empty_df(residuals)) {
      cli_warn("`check = TRUE` needs a fitted model; no diagnostic panels were added.")

    } else {

      .fitted <- .residual <- .std.res <- NULL # for CRAN check only

      qq <- ggplot2::ggplot(residuals, aes(sample = .std.res)) +
        ggplot2::geom_qq_line(alpha = 0.5, fullrange = TRUE) +
        ggplot2::geom_qq(size = 0.2, alpha = 0.5) +
        ggplot2::labs(
          title    = "Normal Q-Q",
          x        = "Theoretical Quantiles",
          y        = "Standardized Residual",
          subtitle = "Points should fall on the line." )

      rf <- ggplot2::ggplot(residuals, aes(x = .fitted, y = .residual)) +
        ggplot2::geom_hline(yintercept = 0, alpha = 0.5, linetype = "dashed") +
        ggplot2::geom_point(size = 0.2, alpha = 0.5) +
        ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
        ggplot2::labs(
          title    = "Residuals vs Fitted",
          x        = "Fitted Value",
          y        = "Residual",
          subtitle = "Flat and horizontal trendline." )

      sl <- ggplot2::ggplot(residuals, aes(x = .fitted, y = sqrt(abs(.std.res)))) +
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

      # The main plot goes last so that it is patchwork's outer plot (and
      # `p$data` is its data); `design` places it above the diagnostics.
      p <- patchwork::wrap_plots(
        A       = patchwork::wrap_plots(diag, nrow = 1),
        B       = patchwork::free(p),
        design  = "B\nA",
        heights = c(3, 1) )

      code <- c(code, "", plot_check_code(plot_var, assign = panel))
    }
  }


  return (plot_finalize(p, params, code, ggdata))
}



#______________________________________________________________
# Attach code, facet info, etc. and hide the build environment.
#______________________________________________________________
plot_finalize <- function (p, params, code, ggdata) {

  for (i in names(params$.plot_attrs))
    attr(p, i) <- params$.plot_attrs[[i]]

  attr(p, 'code') <- add_class(paste(code, collapse = "\n"), 'rbiom_code')

  # The command that produces `data`, for plot_wrap() and assemble_code().
  if (!is.null(attr(ggdata, 'cmd', exact = TRUE)))
    attr(p, 'data_cmd') <- attr(ggdata, 'cmd', exact = TRUE)

  p %<>% add_class('rbiom_plot')

  # Prevent p from including all of this environment's variables.
  if (!inherits(p, 'patchwork')) p$plot_env <- emptyenv()

  return (p)
}



#______________________________________________________________
# A derived table written out as a literal tibble.
#______________________________________________________________
derived_literal_code <- function (var, df, comment = NULL) {

  comment <- switch(
    EXPR = var,
    labels      = "# p-values from the statistics table (see `p$stats$code`), placed above the boxes.",
    brackets    = "# Bracket segments connecting the compared groups.",
    stat_vline  = "# Position(s) at which the statistics were evaluated.",
    taxa_coords = "# Taxa biplot coordinates (see `p$taxa_stats$code`).",
    comment )

  cols <- lapply(as.list(df), deparse_vec, indent = 4)

  c(comment, sprintf("%s <- %s", var, tibble_code(cols)))
}



#______________________________________________________________
# The code for the diagnostic panels added by `check = TRUE`. Mirrors
# the construction in plot_build() line for line.
#______________________________________________________________
plot_check_code <- function (plot_var, assign = FALSE) {
  c(
    "# Model diagnostics: residuals from the fitted trendline(s).",
    "qq <- ggplot(residuals, aes(sample = .std.res)) +",
    "  geom_qq_line(alpha = 0.5, fullrange = TRUE) +",
    "  geom_qq(size = 0.2, alpha = 0.5) +",
    "  labs(title = \"Normal Q-Q\", x = \"Theoretical Quantiles\", y = \"Standardized Residual\", subtitle = \"Points should fall on the line.\")",
    "",
    "rf <- ggplot(residuals, aes(x = .fitted, y = .residual)) +",
    "  geom_hline(yintercept = 0, alpha = 0.5, linetype = \"dashed\") +",
    "  geom_point(size = 0.2, alpha = 0.5) +",
    "  geom_smooth(method = \"gam\", formula = y ~ s(x, bs = \"cs\")) +",
    "  labs(title = \"Residuals vs Fitted\", x = \"Fitted Value\", y = \"Residual\", subtitle = \"Flat and horizontal trendline.\")",
    "",
    "sl <- ggplot(residuals, aes(x = .fitted, y = sqrt(abs(.std.res)))) +",
    "  geom_point(size = 0.2, alpha = 0.5) +",
    "  geom_smooth(method = \"gam\", formula = y ~ s(x, bs = \"cs\")) +",
    "  labs(title = \"Scale-Location\", x = \"Fitted Value\", y = expression(sqrt(abs(` Standardized Residual `))), subtitle = \"Flat and horizontal trendline.\")",
    "",
    "diag <- lapply(list(qq, rf, sl), function (x) {",
    "  x + theme_minimal() + theme(",
    "    plot.title      = element_text(face = \"bold\", size = 10),",
    "    plot.subtitle   = element_text(face = \"italic\", size = 9),",
    "    axis.title      = element_text(size = 9),",
    "    plot.background = element_rect(linewidth = 1, color = \"gray\") )",
    "})",
    "",
    sprintf("%spatchwork::wrap_plots(A = patchwork::wrap_plots(diag, nrow = 1), B = patchwork::free(%s), design = \"B\\nA\", heights = c(3, 1))", if (assign) paste(plot_var, "<- ") else "", plot_var) )
}



#______________________________________________________________
# Cache key for a plot built from a data.frame (`df`): provenance
# attributes are excluded so that identical data from different
# pipelines share an entry. rbiom objects are hashed by content.
#______________________________________________________________
plot_cache_file <- function (fn, params) {
  key <- as.list(params, all.names = TRUE)
  if (is.data.frame(key$df)) key$df <- strip_prov(key$df)
  get_cache_file(fn, key)
}



#______________________________________________________________
# Attach provenance to a plot after any cache lookup: the command that
# made it, and the history of the rbiom object (or of the data.frame)
# it was built from. The statistics tables attached to the plot get
# the same history, since theirs was fixed when the plot was cached.
#______________________________________________________________
plot_prov <- function (p, cmd, biom = NULL, data = NULL) {
  
  if (is.null(p)) return (NULL)
  
  p <- prov_attach(p, cmd, biom, data)
  
  for (a in c('stats', 'taxa_stats')) {
    s <- attr(p, a, exact = TRUE)
    if (is.null(s)) next
    attr(s, 'history') <- attr(p, 'history', exact = TRUE)
    attr(p, a) <- s
  }
  
  return (p)
}



#______________________________________________________________
# Code that applies to `data` the same cleaning that the validators
# applied to `df`: rows with NA (or non-finite x/y) removed, and the
# categorical columns converted to factors without unused levels
# (validate_df_field() does this with refactor()).
#______________________________________________________________
plot_data_prep_code <- function (before, after, cat_cols, num_cols, dv) {
  
  lines <- character()
  cols  <- unique(c(cat_cols, num_cols))
  
  if (nrow(after) < nrow(before))
    lines <- c(lines, sprintf(
      fmt = "%s <- %s[stats::complete.cases(%s[%s])%s, ]",
      dv, dv, dv, deparse_vec(cols),
      paste(collapse = "", sprintf(" & is.finite(%s[[%s]])", dv, double_quote(num_cols))) ))
  
  for (col in cat_cols)
    if (!is.factor(before[[col]]) || !identical(levels(before[[col]]), levels(after[[col]])))
      lines <- c(lines, sprintf("%s[[%s]] <- droplevels(as.factor(%s[[%s]]))", dv, double_quote(col), dv, double_quote(col)))
  
  return (lines)
}



#______________________________________________________________
# apply_p.top() (taxa_table.r), plus the code that applies the same
# taxa filter to `data`.
#______________________________________________________________
apply_p_top <- function (params) {
  
  before <- params$.ggdata
  apply_p.top(params)
  after  <- params$.ggdata
  
  if (!hasName(after, '.taxa') || !hasName(after, '.rank') || nrow(after) == nrow(before))
    return (invisible(params))
  
  dv   <- params$.data_var %||% "data"
  keep <- unique(data.frame(
    .rank = as.character(after[['.rank']]),
    .taxa = as.character(after[['.taxa']]),
    stringsAsFactors = FALSE ))
  
  params$.data_code <- c(
    params$.data_code,
    sprintf("# p.top = %s: keep only the taxa with the most significant p-values (see `p$stats`).", deparse_vec(params$p.top)),
    sprintf("keep <- %s", tibble_code(list(.rank = deparse_vec(keep$.rank), .taxa = deparse_vec(keep$.taxa)))),
    sprintf("%s <- dplyr::semi_join(dplyr::mutate(%s, .rank = as.character(.rank), .taxa = as.character(.taxa)), keep, by = c(\".rank\", \".taxa\"))", dv, dv),
    sprintf("%s[[\".rank\"]] <- factor(%s[[\".rank\"]], levels = %s)", dv, dv, deparse_vec(levels(after[['.rank']]))),
    sprintf("%s[[\".taxa\"]] <- factor(%s[[\".taxa\"]], levels = %s)", dv, dv, deparse_vec(levels(after[['.taxa']]))) )
  
  return (invisible(params))
}



#______________________________________________________________
# The order in which aesthetics are written in generated aes() calls
# (the same order utils.r's aes_toString() uses), so that the object
# and its code agree.
#______________________________________________________________
aes_order <- function (keys) {
  std  <- sub("^color$", "colour", keys)
  sort <- c("x", "xend", "xmin", "xmax", "y", "yend", "ymin", "ymax", "colour", "fill", "shape", "group", "label")
  keys[order(match(std, sort, nomatch = length(sort) + seq_along(std)))]
}



#______________________________________________________________
# For a plot built from a data.frame that no rbiom function produced
# (no `cmd` attribute), the code can't recreate `data`; say so.
#______________________________________________________________
data_origin_code <- function (df) {
  if (inherits(df, 'rbiom') || !is.null(attr(df, 'cmd', exact = TRUE))) return (character())
  sprintf("# data <- %s   # not reproducible from code", obj_repr(df))
}
