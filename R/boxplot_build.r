
#________________________________________________________
# Build a boxplot-family figure from a params environment.
#
# Inputs read from `params` (see stats_boxplot() for their meaning):
#   df, x, y, layers, stat.by, facet.by, colors, shapes, patterns, test,
#   flip, stripe, ci, level, p.adj, p.top, outliers, xlab.angle, p.label,
#   strata, caption, .dots, and optionally y.transform (taxa_boxplot),
#   .labs (default axis labels / caption), and for beta diversity
#   boxplots biom, bdiv, tree, alpha, cpus, seed, permutations.
#
# The data-cleaning steps applied to `df` are recorded in `.data_code`
# so that the generated code reproduces `.ggdata` from `data`.
#________________________________________________________

boxplot_build <- function (params) {

  stopifnot(is_bare_environment(params))

  dv <- params$.data_var %||% "data"
  params$.dv <- dv
  params$.data_code %<>% if.null(character())
  default('p.top', Inf, env = params)
  default('strata', NULL, env = params)


  with(params, {


    #________________________________________________________
    # Validate and pre-process user's arguments.
    #________________________________________________________
    if (inherits(df, 'rbiom')) {
      df <- df$metadata
      .data_code <- c(.data_code, sprintf("%s <- biom$metadata", .dv))
    }

    if (!inherits(df, 'data.frame'))
      cli_abort("`df` must be a data.frame or rbiom object, not {.type {df}}.")

    # Enforce unnamed vectors.
    for (.i in colnames(df)) df[[.i]] %<>% unname()
    remove(".i")

    .df0 <- df

    validate_df_field('x',        col_type = "cat", null_ok = TRUE)
    validate_df_field('y',        col_type = "num")
    validate_df_field('stat.by',  col_type = "cat", null_ok = TRUE)
    validate_df_field('facet.by', col_type = "cat", null_ok = TRUE, max = Inf)

    # `strata` is a sample field: it's in `df` except for beta diversity
    # boxplots, where `df` holds sample pairs and the tests use `biom`.
    if (eq(y, '.distance') && exists('biom', inherits = FALSE)) {
      validate_biom_field('strata', col_type = "cat", null_ok = TRUE)
    } else {
      validate_df_field('strata', col_type = "cat", null_ok = TRUE)
    }

    validate_var_choices('test',       choices = c('auto', 'none', 'wilcox', 'emmeans', 'kruskal', 'adonis2', 'mrpp', 'betadisper'))
    validate_var_choices('p.adj',      choices = p.adjust.methods)
    validate_var_choices('ci',         choices = c('ci', 'range', 'sd', 'se', 'mad'))
    xlab.angle <- as.character(xlab.angle)
    validate_var_choices('xlab.angle', choices = c('auto', '0', '30', '90'))

    validate_var_range('level',   range = c(0.5, 1), n = 1)
    validate_var_range('p.label', range = c(0, 1))

    validate_bool('flip')
    validate_bool('stripe',   default = flip)
    validate_bool('outliers', null_ok = TRUE)
    validate_bool('caption')

    if (exists('y.transform', inherits = FALSE))
      validate_var_choices('y.transform', c('none', 'log1p', 'sqrt'))


    #________________________________________________________
    # Consistent naming for plot_* functions.
    #________________________________________________________
    .ggdata <- as_tibble(df)
    .xcol   <- x
    .ycol   <- y
    .xmode  <- "factor"


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
    # Drop rows without a finite y value. Record the cleaning (NA
    # removal and factor conversion by the validators above, plus
    # this) in the code.
    #________________________________________________________
    .ggdata <- .ggdata[is.finite(.ggdata[[.ycol]]),]

    .data_code <- c(.data_code, plot_data_prep_code(
      before   = .df0,
      after    = .ggdata,
      cat_cols = unique(c(.xcol, stat.by, facet.by, if (!eq(.ycol, '.distance')) strata)),
      num_cols = .ycol,
      dv       = .dv ))

    remove("df", "x", "y", ".df0")

  })


  #________________________________________________________
  # Provenance-tracked ggplot2 functions.
  #________________________________________________________
  .position_dodge       <- P('ggplot2::position_dodge')
  .position_jitter      <- P('ggplot2::position_jitter')
  .position_jitterdodge <- P('ggplot2::position_jitterdodge')



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
  # Default labels: the response label from the table function, plus
  # any set by the calling function. The user's `labs.*` win.
  #________________________________________________________
  if (eq(params$.ycol, attr(params$.ggdata, 'response', exact = TRUE)))
    set_layer(params, 'labs', y = attr(params$.ggdata, 'resp_label', exact = TRUE))

  if (length(params$.labs) > 0)
    set_layer(params, 'labs', params$.labs)

  if (is.null(params$.xcol))
    set_layer(params, 'labs', x = NULL)


  #________________________________________________________
  # Transform the y-axis (taxa_boxplot's `y.transform`).
  #________________________________________________________
  if (!is.null(params$y.transform) && !eq(params$y.transform, 'none'))
    set_layer(params, 'yaxis', transform = params$y.transform)


  #________________________________________________________
  # All layers share `x`, and `y` metadata fields.
  #________________________________________________________
  set_layer(
    params = params,
    layer  = 'ggplot',
    'mapping|x' = if.null(params$.xcol, ""),
    'mapping|y' = params$.ycol )



  #________________________________________________________
  # Help multiple layers overlay well together
  #________________________________________________________

  dodged    <- !eq(params$.xcol, params$stat.by)
  patterned <- !is.null(params$patterns)
  dotted    <- any(has_layer(params, c('dot', 'strip')))
  beeswarm  <- pkg_available('ggbeeswarm')
  dodge     <- .position_dodge(width = 0.8)


  #________________________________________________________
  # Adjust default point size based on number of points.
  #________________________________________________________
  ptsize <- local({
    n <- nrow(params$.ggdata)
    round(log(1002 - min(1000, max(10, n))) / 6, digits = 2)
  })


  #________________________________________________________
  # ggbeeswarm's layers warn about values that differ only in the
  # last few digits. Round for them, and record it in the code.
  #________________________________________________________
  if (dotted && beeswarm) {
    params$.ggdata[[params$.ycol]] %<>% round(12)
    params$.data_code <- c(params$.data_code, sprintf(
      fmt = "%s[[%s]] <- round(%s[[%s]], 12)",
      dv, double_quote(params$.ycol), dv, double_quote(params$.ycol) ))
  }



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
    transform <- layers[['yaxis']][['transform']]

    if (is_null(transform)) {
      set_layer(params, 'bar', fun="mean")

    } else if (eq(transform, 'sqrt')) {
      set_layer(params, 'bar', fun = as.cmd(function (y) sqrt(mean(y * y))))

    } else if (eq(transform, 'log1p')) {
      set_layer(params, 'bar', fun = as.cmd(function (y) log1p(mean(exp(y) - 1))))

    } else {
      cli_abort("Bar charts can only be re-scaled using 'sqrt' or 'log1p', not {.val {transform}}.")
    }

    remove("transform")
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
      set_layer(params, 'box', position=dodge)


    # Show/hide/customize boxplot outliers.
    #________________________________________________________
    if (is.null(params$outliers))
      params$outliers <- !dotted

    if (isFALSE(params$outliers)) {
      set_layer(params, 'box', outlier.shape = NA)
    } else {
      set_layer(params, 'box', outlier.size = ptsize)
    }

  }



  #________________________________________________________
  # Dot ====
  # ggbeeswarm::geom_beeswarm, or geom_point() with jitter when
  # ggbeeswarm is not installed.
  #________________________________________________________
  if (has_layer(params, 'dot')) {

    if (beeswarm) {

      set_layer(params, 'dot', cex = 0.5, size = ptsize, method = "center", corral = "wrap", corral.width = 0.9)

      # Keep the exact y values (and avoid ggbeeswarm's warning about it).
      if ('preserve.data.axis' %in% attr(attr(layers[['dot']], 'function', exact = TRUE), 'formalArgs', exact = TRUE))
        set_layer(params, 'dot', preserve.data.axis = TRUE)

      if (dodged)
        set_layer(params, 'dot', dodge.width = 0.8)

    } else {

      set_layer(params, 'dot', size = ptsize, position = boxplot_jitter(params, dodged, .position_jitter, .position_jitterdodge))
    }

    if (patterned)
      set_layer(params, 'dot', fill = "black")

    if (any(has_layer(params, c("violin", "box", "bar"))))
      set_layer(params, 'dot', color = "black", fill = "black")

    if (any(has_layer(params, c("crossbar", "errorbar", "linerange", "pointrange"))))
      set_layer(params, 'dot', alpha = 0.4)
  }


  #________________________________________________________
  # Strip ====
  # ggbeeswarm::geom_quasirandom, or geom_point() with jitter when
  # ggbeeswarm is not installed.
  #________________________________________________________
  if (has_layer(params, 'strip')) {

    set_layer(params, 'strip', size = ptsize)

    if (beeswarm) {

      if (dodged)
        set_layer(params, 'strip', dodge.width = 0.8)

    } else {

      set_layer(params, 'strip', position = boxplot_jitter(params, dodged, .position_jitter, .position_jitterdodge))
    }

    if (patterned)
      set_layer(params, 'strip', fill = "black")

    if (any(has_layer(params, c("violin", "box", "bar"))))
      set_layer(params, 'strip', color = "black", fill = "black")

    if (any(has_layer(params, c("crossbar", "errorbar", "linerange", "pointrange"))))
      set_layer(params, 'strip', alpha = 0.4)
  }


  #________________________________________________________
  # CI ====
  # The intervals themselves are computed after the statistics, in
  # boxplot_intervals(), so that `p.top` has already dropped taxa.
  #________________________________________________________
  if (any(has_layer(params, c("crossbar", "errorbar", "linerange", "pointrange")))) {

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

    remove("args")
  }


  #________________________________________________________
  # Shade background of every other x-axis position.
  #________________________________________________________
  if (isTRUE(params$stripe) && !is.null(params$.xcol) && nlevels(params$.ggdata[[params$.xcol]]) > 1) {

    stripe_x <- seq(2, nlevels(params$.ggdata[[params$.xcol]]), 2)

    set_layer(
      params = params,
      layer  = 'stripe',
      'geom'  = "rect",
      'xmin'  = stripe_x - 0.5,
      'xmax'  = stripe_x + 0.5,
      'ymin'  = -Inf,
      'ymax'  = Inf,
      'fill'  = "black",
      'color' = NA,
      'alpha' = 0.05 )

    remove("stripe_x")
  }



  #________________________________________________________
  # Build the plot.
  #________________________________________________________
  fig <- params %>%
    plot_facets() %>%
    boxplot_stats() %>%
    boxplot_intervals() %>%
    plot_build()


  return (fig)
}



#________________________________________________________
# The jitter position used in place of ggbeeswarm's layers.
#________________________________________________________
boxplot_jitter <- function (params, dodged, .position_jitter, .position_jitterdodge) {

  # position_jitterdodge() needs an aesthetic to dodge by.
  dodge_aes <- !is.null(params$stat.by) && (!is.null(params$colors) || !is.null(params$shapes) || !is.null(params$patterns))

  if (dodged && dodge_aes) {
    .position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.8, seed = 0)
  } else {
    .position_jitter(width = 0.25, height = 0, seed = 0)
  }
}



#________________________________________________________
# The center/min/max table drawn by the crossbar, errorbar, linerange,
# and pointrange layers: one row per x / stat.by / facet group, computed
# with dplyr so that the code text below is the same computation.
#________________________________________________________
boxplot_intervals <- function (params) {

  if (!any(has_layer(params, c("crossbar", "errorbar", "linerange", "pointrange"))))
    return (invisible(params))

  ggdata <- params$.ggdata
  ycol   <- params$.ycol
  ci     <- params$ci
  level  <- params$level
  dv     <- params$.data_var %||% "data"
  by     <- unique(c(params$.xcol, params$stat.by, params$facet.by))

  if (is_empty_df(ggdata)) return (invisible(params))


  #________________________________________________________
  # One function per `ci` option, and the same function as code.
  #________________________________________________________
  .y <- .ymin <- .ymax <- NULL # for CRAN check only

  fn <- switch(
    EXPR = ci,

    'ci' = list(
      fun = function (v) {
        tt <- tryCatch(stats::t.test(v, conf.level = level), error = function (e) list(estimate = NA_real_, conf.int = c(NA_real_, NA_real_)))
        tibble(.y = unname(tt$estimate), .ymin = tt$conf.int[[1]], .ymax = tt$conf.int[[2]]) },
      code = c(
        sprintf("tt <- tryCatch(stats::t.test(v, conf.level = %s), error = function (e) list(estimate = NA_real_, conf.int = c(NA_real_, NA_real_)))", deparse_vec(level)),
        "tibble::tibble(%s = unname(tt$estimate), .ymin = tt$conf.int[[1]], .ymax = tt$conf.int[[2]])" )),

    'range' = list(
      fun  = function (v) tibble(.y = mean(v), .ymin = min(v), .ymax = max(v)),
      code = "tibble::tibble(%s = mean(v), .ymin = min(v), .ymax = max(v))" ),

    'mad' = list(
      fun = function (v) {
        med <- stats::median(v); dev <- stats::mad(v, med)
        tibble(.y = med, .ymin = med - dev, .ymax = med + dev) },
      code = c(
        "med <- stats::median(v); dev <- stats::mad(v, med)",
        "tibble::tibble(%s = med, .ymin = med - dev, .ymax = med + dev)" )),

    'sd' = list(
      fun = function (v) {
        avg <- mean(v); dev <- stats::sd(v)
        tibble(.y = avg, .ymin = avg - dev, .ymax = avg + dev) },
      code = c(
        "avg <- mean(v); dev <- stats::sd(v)",
        "tibble::tibble(%s = avg, .ymin = avg - dev, .ymax = avg + dev)" )),

    'se' = list(
      fun = function (v) {
        avg <- mean(v); dev <- sqrt(stats::var(v) / length(v))
        tibble(.y = avg, .ymin = avg - dev, .ymax = avg + dev) },
      code = c(
        "avg <- mean(v); dev <- sqrt(stats::var(v) / length(v))",
        "tibble::tibble(%s = avg, .ymin = avg - dev, .ymax = avg + dev)" )) )


  #________________________________________________________
  # Compute.
  #________________________________________________________
  vline <- ggdata
  if (length(by) > 0) vline <- dplyr::group_by(vline, !!!rlang::syms(by))
  vline <- vline %>%
    dplyr::group_modify(function (d, keys) {
      v <- d[[ycol]]
      dplyr::rename(fn$fun(v), !!ycol := .y)
    }) %>%
    dplyr::ungroup()

  # A y-axis transform with a restricted domain can't draw bounds outside
  # it (sqrt of a negative CI bound is NaN): squish them to its edge.
  floor <- interval_floor(params$layers[['yaxis']][['transform']], ggdata[[ycol]])

  if (!is.null(floor))
    for (col in c(ycol, '.ymin', '.ymax'))
      vline[[col]] <- pmax(vline[[col]], floor)

  attr(ggdata, 'vline') <- as_tibble(vline)
  params$.ggdata        <- ggdata


  #________________________________________________________
  # Record.
  #________________________________________________________
  body <- sub("%s", coan(ycol), fn$code, fixed = TRUE)

  params$.derived_code[['vline']] <- paste(collapse = "\n", c(
    "# Center and interval of each group, for the crossbar / errorbar / linerange / pointrange layers.",
    sprintf("intervals <- %s %%>%%", dv),
    if (length(by) > 0) sprintf("  dplyr::group_by(%s) %%>%%", paste(coan(by), collapse = ", ")),
    "  dplyr::group_modify(function (d, keys) {",
    sprintf("    v <- d[[%s]]", double_quote(ycol)),
    paste0("    ", body),
    "  }) %>%",
    if (is.null(floor)) "  dplyr::ungroup()" else c(
      "  dplyr::ungroup() %>%",
      sprintf("  dplyr::mutate(dplyr::across(c(%s, .ymin, .ymax), function (v) pmax(v, %s)))   # the y-axis transform's domain",
        coan(ycol), deparse_vec(floor)) ) ))


  return (invisible(params))
}



#________________________________________________________
# The smallest value a y-axis transform can display, or NULL when it
# can show any value. sqrt and log1p get 0; log axes get the smallest
# positive y value in the data.
#________________________________________________________
interval_floor <- function (transform, y) {

  name <- if (is.character(transform)) transform[[1]]
  else if (inherits(transform, c('transform', 'trans'))) transform$name
  else NULL

  if (!is_scalar_character(name)) return (NULL)
  name <- sub("^[xy]_", "", name)

  if (name %in% c("sqrt", "log1p")) return (0)

  if (name %in% c("log", "log10", "log2") || startsWith(name, "log-")) {
    pos <- y[is.finite(y) & y > 0]
    return (if (length(pos) > 0) min(pos) else NULL)
  }

  return (NULL)
}
