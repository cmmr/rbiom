#________________________________________________________
# A plot is assembled from "layers": a named environment (`params$layers`)
# in which each entry holds the arguments for one ggplot2 call. The
# plot_* functions populate it with init_layers() / set_layer(), and
# plot_build() turns it into a ggplot object plus the code that rebuilds
# that object from `data`.
#
# Arguments given in `...` by the user are routed to layers by prefix
# (`b.alpha = 0.5` goes to the bar layer) or, unprefixed, to every layer
# whose function has an argument of that name (`scales = "free"` goes to
# facet_wrap()). add_layer() records which ones were used so that
# plot_build() can reject the rest with a helpful message.
#
# Optional packages (ggbeeswarm, ggrepel, ggdensity, fillpattern) are used
# when installed; the dot/strip/taxon layers fall back to ggplot2 geoms.
#________________________________________________________



#________________________________________________________
# Is an optional package available? The `rbiom.disable_pkgs` option lets
# the tests exercise the fallbacks without uninstalling anything.
#________________________________________________________
pkg_available <- function (pkg) {
  if (pkg %in% getOption('rbiom.disable_pkgs', character())) return (FALSE)
  nzchar(system.file(package = pkg))
}



#________________________________________________________
# Converts user's `layers` spec to `layers` environment.
#________________________________________________________
init_layers <- function (
    params = parent.frame(), choices = NULL, var = "layers", do_init = NULL ) {


  layer_names <- do_init
  if (!is_null(params$facet.by)) layer_names %<>% c('facet')
  if (!is.null(choices))         layer_names %<>% c(validate_layers(var, choices, env = params))



  #________________________________________________________
  # Ignore shapes/patterns without applicable layers.
  #________________________________________________________
  if (!any(c('dot', 'strip', 'pointrange', 'point') %in% layer_names)) params$shapes   <- NULL
  if (!any(c('box', 'bar', 'violin', 'stack')       %in% layer_names)) params$patterns <- NULL


  #________________________________________________________
  # Resolve `colors`, `patterns`, and `shapes`.
  #________________________________________________________

  if (is_palette(params$colors))
    params$colors <- color_palette(params$colors)

  for (i in c('colors', 'shapes', 'patterns')) {

    x <- params[[i]]

    if (is.null(x)) next
    if (isFALSE(x)) { params[[i]] <- NULL; next }


    # Use `color.by` instead of `stat.by`
    by <- sub("s$", ".by", i)
    by <- if (hasName(params, by)) params[[by]] else params$stat.by
    if (is.null(by)) { params[[i]] <- NULL; next }
    by_vals <- params$.ggdata[[by]]
    if (is.null(by_vals)) { params[[i]] <- NULL; next }

    # Gradient colors finish up here.
    if (!is.factor(by_vals)) next


    stat_lvls  <- levels(by_vals)
    stat_nlvls <- length(stat_lvls)


    if (isTRUE(x))
      x <- switch(
        EXPR = i,
        'colors'   = get_n_colors(n = stat_nlvls),
        'shapes'   = get_n_shapes(n = stat_nlvls),
        'patterns' = get_n_patterns(n = stat_nlvls) )


    # Re-order named aesthetic values to match levels() order.
    if (!is.null(names(x))) {

      xn <- names(x)
      validate_var_choices('xn', evar = i, choices = stat_lvls, max = Inf)
      names(x) <- xn
      remove("xn")

      if (length(missing <- setdiff(stat_lvls, names(x))) > 0)
        cli_abort("Missing `{i}` for {.field {by}} {qty(length(missing))} level{?s} {.val {missing}}.")

      if (length(dups <- unique(names(x)[duplicated(names(x))])) > 0)
        cli_abort("Duplicated {qty(length(dups))} name{?s} in `{i}`: {.val {dups}}.")

      x <- as.vector(x[stat_lvls])
    }


    params[[i]] <- rep_len(x, stat_nlvls)
  }

  remove(list = c("i", "x", "by", "by_vals", "stat_lvls", "stat_nlvls") %>% intersect(ls()))


  if (!is.null(params$patterns) && !pkg_available('fillpattern'))
    package_missing('fillpattern', 'to fill with patterns (`patterns`)')


  if (!is.null(names(params$.dots)))
    if (any(startsWith(names(params$.dots), 'labs.')))
      layer_names %<>% c('labs')


  init_layer_names <- layer_names %>%
    c('ggplot', ., 'xaxis', 'yaxis', 'theme_rbiom', 'theme') %>%
    unique()


  #________________________________________________________
  # Create and populate the layers environment.
  #________________________________________________________
  params$.plot_attrs %<>% if.null(list())
  params$.dots_used  %<>% if.null(character())
  params$layers      <- rlang::new_environment()

  for (layer_name in init_layer_names)
    add_layer(params, layer_name)


  return (invisible(layer_names))
}




#________________________________________________________
# Finds parent's `layers` variable for a particular layer
#________________________________________________________
has_layer <- function (params, layer) {
  sapply(layer, hasName, x = params$layers)
}



#________________________________________________________
# Deletes an entire layer.
#________________________________________________________
del_layer <- function (params, layer) {

  for (i in layer)
    if (hasName(params$layers, i))
      remove(list = i, pos = params$layers)

  return (invisible(NULL))
}



#________________________________________________________
# Finds `layers` variable in parent env, then updates it.
#________________________________________________________
set_layer <- function (params, layer, ..., .fn = NULL, .overwrite = FALSE) {

  layers <- params$layers


  if (!hasName(layers, layer))
    add_layer(params, layer, fn = .fn)

  dots <- list(...)


  # Expand lists. For cases of e.g.
  # set_layer(params, 'box', fill="white", list(linewidth=2))
  keyvals <- list()
  for (i in seq_along(dots)) {
    if (isTRUE(nzchar(key <- names(dots)[[i]]))) {
      if (is_null(dots[[i]])) {
        keyvals[key] <- list(NULL)
      } else {
        keyvals[[key]] <- dots[[i]]
      }

    } else {
      keyvals <- c(keyvals, dots[[i]])
    }
  }


  for (key in names(keyvals)) {
    val <- keyvals[[key]]
    key <- strsplit(key, "|", fixed = TRUE)[[1]]

    # set_layer(params, 'bar', fill = "green")
    if (length(key) == 1) {
      if (!hasName(layers[[layer]], key) || .overwrite)
        if (is_null(val)) {
          layers[[layer]][key] <- list(NULL)
        } else {
          layers[[layer]][[key]] <- val
        }

    # set_layer(params, 'bar', "mapping|fill" = "Body Site")
    } else if (length(key) == 2) {
      key1 <- key[[1]]
      key2 <- key[[2]]

      if (!hasName(layers[[layer]], key1))
        layers[[layer]][[key1]] <- list()

      if (!hasName(layers[[layer]][[key1]], key2) || .overwrite)
        if (is_null(val)) {
          layers[[layer]][[key1]][key2] <- list(NULL)
        } else {
          layers[[layer]][[key1]][[key2]] <- val
        }

    } else {
      cli_abort("Key vector is too long. Key {.val {key}} split on '|'.")
    }

  }

  return (invisible(NULL))
}



#________________________________________________________
# The ggplot2 (or optional package) function behind each layer.
#________________________________________________________
layer_function <- function (layer, params, fn = NULL) {

  xmode     <- params$.xmode
  patterned <- !is_null(params$patterns)
  facetDims <- length(params$facet.by)

  switch(
    EXPR = layer,
    'ggplot'       = P('ggplot2::ggplot'),
    'arrow'        = P('ggplot2::geom_segment'),
    'bar'          = P('ggplot2::geom_bar'),
    'box'          = P('ggplot2::geom_boxplot'),
    'brackets'     = P('ggplot2::geom_segment'),
    'ci'           = P('ggplot2::geom_errorbar'),
    'color'        = P('ggplot2::scale_color_manual'),
    'confidence'   = P('ggplot2::geom_ribbon'),
    'crossbar'     = P('ggplot2::geom_crossbar'),
    'density'      = P('ggdensity::geom_hdr'),
    'dot'          = if (pkg_available('ggbeeswarm')) P('ggbeeswarm::geom_beeswarm') else P('ggplot2::geom_point'),
    'ellipse'      = P('ggplot2::stat_ellipse'),
    'errorbar'     = P('ggplot2::geom_errorbar'),
    'facet'        = if (facetDims == 2) P('ggplot2::facet_grid') else P('ggplot2::facet_wrap'),
    'fill'         = if (patterned) P('fillpattern::scale_fill_pattern') else P('ggplot2::scale_fill_manual'),
    'flip'         = P('ggplot2::coord_flip'),
    'hline'        = P('ggplot2::geom_hline'),
    'labs'         = P('ggplot2::labs'),
    'linerange'    = P('ggplot2::geom_linerange'),
    'mean'         = P('ggplot2::geom_point'),
    'name'         = P('ggplot2::geom_text'),
    'point'        = P('ggplot2::geom_point'),
    'pointrange'   = P('ggplot2::geom_pointrange'),
    'rect'         = P('ggplot2::geom_rect'),
    'residual'     = P('ggplot2::geom_segment'),
    'shape'        = P('ggplot2::scale_shape_manual'),
    'spider'       = P('ggplot2::geom_segment'),
    'stack'        = P('ggplot2::geom_col'),
    'stats_bg'     = P('ggplot2::annotate'),
    'stats_label'  = P('ggplot2::geom_label'),
    'stats_text'   = P('ggplot2::geom_text'),
    'stats_vline'  = P('ggplot2::geom_vline'),
    'strip'        = if (pkg_available('ggbeeswarm')) P('ggbeeswarm::geom_quasirandom') else P('ggplot2::geom_point'),
    'stripe'       = P('ggplot2::annotate'),
    'taxon'        = if (pkg_available('ggrepel')) P('ggrepel::geom_label_repel') else P('ggplot2::geom_label'),
    'theme'        = P('ggplot2::theme'),
    'theme_rbiom'  = P('rbiom::theme_rbiom'),
    'trend'        = P('ggplot2::geom_line'),
    'violin'       = P('ggplot2::geom_violin'),
    'vline'        = P('ggplot2::geom_vline'),
    'xaxis'        = if (eq(xmode, "factor")) P('ggplot2::scale_x_discrete') else P('ggplot2::scale_x_continuous'),
    'yaxis'        = if (eq(params$.ymode, "factor")) P('ggplot2::scale_y_discrete') else P('ggplot2::scale_y_continuous'),
    'zero'         = P('ggplot2::geom_vline'),
    if (is.null(fn)) P(paste0('ggplot2::', layer)) else get(fn) )
}



#________________________________________________________
# Layers whose geoms draw lines or outlines, where `size` must be given
# as `linewidth`. (Point and text layers keep `size`; pointrange uses
# both, so it's excluded.)
#________________________________________________________
LINE_LAYERS <- c(
  'arrow', 'bar', 'box', 'brackets', 'ci', 'confidence', 'crossbar', 'ellipse',
  'errorbar', 'hline', 'linerange', 'rect', 'residual', 'spider', 'stack',
  'stats_vline', 'trend', 'violin', 'vline', 'zero' )



#________________________________________________________
# Which derived table (an attribute of `.ggdata`) a layer draws from.
# See plot_build() for how these become `labels <- ...` etc. in the code.
#________________________________________________________
layer_src <- function (layer) {
  switch(
    EXPR = layer,
    'arrow'        = "taxa_coords",
    'brackets'     = "stat_brackets",
    'confidence'   = "fit",
    'crossbar'     = "vline",
    'errorbar'     = "vline",
    'linerange'    = "vline",
    'mean'         = "taxa_coords",
    'pointrange'   = "vline",
    'residual'     = "residual",
    'spider'       = "spider",
    'stats_label'  = "stat_labels",
    'stats_text'   = "stat_labels",
    'stats_vline'  = "stat_vline",
    'taxon'        = "taxa_coords",
    'trend'        = "fit",
    NULL )
}



#________________________________________________________
# The `...` prefix(es) that route an argument to a layer.
#________________________________________________________
layer_regex <- function (layer) {
  switch(
    EXPR = layer,
    'ggplot'       = "^g(|gplot)\\.",
    'arrow'        = "^a(|rrow)\\.",
    'bar'          = "^(b|bar)\\.",
    'box'          = "^box\\.",
    'brackets'     = "^brackets\\.",
    'ci'           = "^ci\\.",
    'color'        = "^color\\.",
    'confidence'   = "^c(|onfidence)\\.",
    'crossbar'     = "^c(|rossbar)\\.",
    'density'      = "^d(|ensity)\\.",
    'dot'          = "^(pt|d|dot)\\.",
    'errorbar'     = "^e(|rrorbar)\\.",
    'ellipse'      = "^e(|llipse)\\.",
    'facet'        = "^f(|acet)\\.",
    'fill'         = "^(fill|pattern|patterns)\\.",
    'flip'         = "^flip\\.",
    'hline'        = "^[rh](|line)\\.",
    'labs'         = "^labs\\.",
    'linerange'    = "^l(|inerange)\\.",
    'mean'         = "^m(|ean)\\.",
    'name'         = "^n(|ame)\\.",
    'pointrange'   = "^(pt|p|pointrange)\\.",
    'point'        = "^(pt|p|point)\\.",
    'residual'     = "^r(|esidual)\\.",
    'shape'        = "^shape\\.",
    'size'         = "^size\\.",
    'spider'       = "^s(|pider)\\.",
    'stack'        = "^s(|tack)\\.",
    'stats_label'  = "^stats\\.",
    'stats_text'   = "^stats\\.",
    'strip'        = "^(pt|s|strip)\\.",
    'taxon'        = "^taxon\\.",
    'theme'        = "^theme\\.",
    'theme_rbiom'  = "^theme_(rbiom|bw)\\.",
    'trend'        = "^t(|rend)\\.",
    'vline'        = "^[rv](|line)\\.",
    'violin'       = "^v(|iolin)\\.",
    'xaxis'        = "^x(|axis)\\.",
    'yaxis'        = "^y(|axis)\\.",
    'zero'         = "^z(|ero)\\.",
    sprintf("^%s\\.", layer) )
}



#________________________________________________________
# Set up a layer with its function, function name, etc.
# Uses attr(layers, 'params') to add in user's params.
# Can set the underlying function using 'fn' parameter.
#________________________________________________________
add_layer <- function (params, layer, fn = NULL) {

  stopifnot(is_scalar_character(layer))
  stopifnot(is_bare_environment(params))

  layers <- params$layers

  if (hasName(layers, layer))
    cli_abort("Cannot add layer {.val {layer}} twice.")


  # Get the plotting function's `...` arguments.
  #________________________________________________________
  layer_params <- params$.dots


  # Ensure NULLs can be inserted into a list
  #________________________________________________________
  for (i in seq_along(layer_params))
    if (is_null(layer_params[[i]]))
      layer_params[[i]] <- list(NULL)


  # Arrays and functions need list wrappers too.
  # plot(yaxis.expand=1:4, yaxis.transform=sqrt)
  #________________________________________________________
  for (i in seq_along(layer_params))
    if (length(layer_params[[i]]) > 1 || is.function(layer_params[[i]]))
      layer_params[[i]] <- list(layer_params[[i]])


  result <- list()
  attr(result, 'function') <- layer_function(layer, params, fn)
  attr(result, 'src')      <- layer_src(layer)

  regex <- layer_regex(layer)
  used  <- character()


  # Unprefixed dot arguments, e.g. 'scales'="free_x"
  #________________________________________________________
  layer_fun <- attr(result, 'function', exact = TRUE)
  f_args    <- intersect(
    x = names(layer_params),
    y = attr(layer_fun, 'formalArgs') %||% formalArgs(layer_fun) )
  for (i in f_args)
    if (!isFALSE(attr(layer_params[[i]], 'display'))) {
      result[i] <- layer_params[[i]]
      used      <- c(used, i)
    }


  # Prefixed dot arguments, e.g. 'facet.scales'="free_x"
  #________________________________________________________
  for (i in grep(regex, names(layer_params), value = TRUE))
    if (!isFALSE(attr(layer_params[[i]], 'display')))
      if (!i %in% c('p.top', 'p.label', 'p.adj')) {
        result[sub(regex, "", i, perl = TRUE)] <- layer_params[[i]]
        used <- c(used, i)
      }


  # ggplot2 >= 3.4 sizes lines with `linewidth`; `size` is deprecated for
  # them. Translate it, e.g. `t.size = 2` for a thicker trendline.
  if (layer %in% LINE_LAYERS && hasName(result, 'size')) {
    if (!hasName(result, 'linewidth')) result['linewidth'] <- result['size']
    result[['size']] <- NULL
  }

  params$.dots_used <- union(params$.dots_used, used)
  layers[[layer]]   <- result


  return (invisible(NULL))
}



#________________________________________________________
# Reject `...` arguments that no layer consumed, naming the prefixes
# that would have worked. Called by plot_build() once every layer exists.
#________________________________________________________
check_unused_dots <- function (params, exempt = NULL) {

  dots <- names(params$.dots)
  if (length(dots) == 0) return (invisible(NULL))

  unused <- setdiff(dots, c(params$.dots_used, exempt, 'ci.min', 'ci.max'))
  unused <- unused[nzchar(unused)]
  if (length(unused) == 0) return (invisible(NULL))

  layers   <- setdiff(env_names(params$layers), c('ggplot', 'stats_bg', 'stripe', 'stats_text', 'stats_label', 'stats_vline', 'brackets'))
  prefixes <- paste0(sort(c(layers, if (!'labs' %in% layers) 'labs')), ".")

  fn <- params$.fn %||% "this function"

  cli_abort(c(
    'x' = "{.val {unused}} {?is/are} not {?an argument/arguments} of {fn} or of any of its layers.",
    'i' = "Prefix an argument with a layer name to pass it to that layer's ggplot2 function, e.g. {.code {prefixes[[1]]}size = 2}. Available prefixes: {.val {prefixes}}.",
    'i' = "Groups are colored by `stat.by`; see also `colors`, `shapes`, and `facet.by`." ))
}
