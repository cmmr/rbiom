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
    by <- params$.ggdata[[by]]
    if (is.null(by)) { params[[i]] <- NULL; next }
    
    # Gradient colors finish up here.
    if (!is.factor(by)) next
    
    
    stat_lvls  <- levels(by)
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
        cli_abort("Missing `{i}` for {stat.by} {qty(missing)} level{?s} {.val {missing}}.")
      
      if (length(dups <- unique(names(x)[duplicated(names(x))])) > 0)
        cli_abort("Duplicated {qty(dups)} name{?s} in `{i}`: {.val {dups}}.")
      
      x <- as.vector(x[stat_lvls])
    }
    
    
    params[[i]] <- rep_len(x, stat_nlvls)
  }
  
  remove(list = c("i", "x", "by", "stat_lvls", "stat_nlvls") %>% intersect(ls()))
  
  
  
  
  if (!is.null(names(params$.dots)))
    if (any(startsWith(names(params$.dots), 'labs.')))
      layer_names %<>% c('labs')
  
  
  init_layer_names <- layer_names %>%
    c('ggplot', ., 'xaxis', 'yaxis', 'theme', 'theme_bw') %>%
    unique()
  
  
  #________________________________________________________
  # Create and populate the layers environment.
  #________________________________________________________
  params$.plot_attrs %<>% if.null(list())
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
# Set up a layer with its function, function name, etc.
# Uses attr(layers, 'params') to add in user's params.
# Can set the underlying function using 'fn' parameter.
#________________________________________________________
add_layer <- function (params, layer, fn = NULL) {
  
  stopifnot(is_scalar_character(layer))
  stopifnot(is_bare_environment(params))
  
  layers <- params$layers
  stopifnot(is_bare_environment(params))
  
  if (hasName(layers, layer))
    cli_abort("Cannot add layer {.val {layer}} twice.")
  
  
  xmode     <- params$.xmode
  patterned <- !is_null(params$patterns)
  facetDims <- length(params$facet.by)
  
  
  
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
  
  
  
  attr(result, 'function') <- switch(
    EXPR = layer,
    'ggplot'       = ggplot,
    'arrow'        = geom_segment,
    'bar'          = geom_bar,
    'box'          = geom_boxplot,
    'brackets'     = geom_segment,
    'color'        = scale_color_manual,
    'confidence'   = geom_ribbon,
    'crossbar'     = geom_crossbar,
    'density'      = geom_hdr,
    'dot'          = geom_beeswarm,
    'errorbar'     = geom_errorbar,
    'ellipse'      = stat_ellipse,
    'facet'        = if (facetDims == 2) facet_grid else facet_wrap,
    'fill'         = if (patterned) scale_fill_pattern else scale_fill_manual,
    'flip'         = coord_flip,
    # 'free_y'     = facetted_pos_scales,
    'hline'        = geom_hline,
    'label'        = geom_label,
    'labs'         = labs,
    'linerange'    = geom_linerange,
    'mean'         = geom_point,
    'name'         = geom_text,
    'pointrange'   = geom_pointrange,
    'point'        = geom_point,
    'rect'         = geom_rect,
    'residual'     = geom_segment,
    'shape'        = scale_shape_manual,
    'smooth'       = stat_smooth,
    'spider'       = geom_segment,
    'stack'        = geom_col,
    'stats_bg'     = geom_rect,
    'stats_ggtext' = geom_textbox,
    'stats_label'  = geom_label,
    'stats_text'   = geom_text,
    'stats_vline'  = geom_vline,
    'strip'        = geom_quasirandom,
    'stripe'       = geom_rect,
    'table'        = annotation_custom,
    'taxon'        = geom_label_repel,
    'theme'        = theme,
    'theme_bw'     = theme_bw,
    'topo'         = geom_hdr_lines,
    'trend'        = geom_line,
    'violin'       = geom_violin,
    'vline'        = geom_vline,
    'xaxis'        = if (xmode == "factor") scale_x_discrete else scale_x_continuous,
    'yaxis'        = scale_y_continuous,
    if (is.null(fn)) get(layer) else get(fn) )
  
  
  attr(result, 'src') <- switch(
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
    'stats_ggtext' = "stat_labels",
    'stats_label'  = "stat_labels",
    'stats_text'   = "stat_labels",
    'stats_vline'  = "stat_vline",
    'taxon'        = "taxa_coords",
    'trend'        = "fit" )
  
  
  result[['data']] <- switch(
    EXPR = layer,
    'arrow'        = ~ attr(., "taxa_coords"),
    'brackets'     = ~ attr(., "stat_brackets"),
    'confidence'   = ~ attr(., "fit"),
    'crossbar'     = ~ attr(., "vline"),
    'errorbar'     = ~ attr(., "vline"),
    'linerange'    = ~ attr(., "vline"),
    'mean'         = ~ attr(., "taxa_coords"),
    'pointrange'   = ~ attr(., "vline"),
    'residual'     = ~ attr(., "residual"),
    'spider'       = ~ attr(., "spider"),
    'stats_ggtext' = ~ attr(., "stat_labels"),
    'stats_label'  = ~ attr(., "stat_labels"),
    'stats_text'   = ~ attr(., "stat_labels"),
    'stats_vline'  = ~ attr(., "stat_vline"),
    'taxon'        = ~ attr(., "taxa_coords"),
    'trend'        = ~ attr(., "fit") )
  
  
  regex <- switch(
    EXPR = layer,
    'ggplot'       = "^g(|gplot)\\.",
    'arrow'        = "^a(|rrow)\\.",
    'bar'          = "^(r|bar)\\.",
    'box'          = "^b(|ox)\\.",
    'brackets'     = "^h(|line)\\.",
    'cladelab'     = "^c(|lade|ladelab)\\.",
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
    'heat'         = "^h(|eat|eatmap)\\.",
    'hexpand'      = "^h(|expand)\\.",
    'hline'        = "^[rh](|line)\\.",
    'label'        = "^l(|abel)\\.",
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
    'stats_ggtext' = "^stats\\.",
    'stats_label'  = "^stats\\.",
    'stats_text'   = "^stats\\.",
    'strip'        = "^(pt|s|strip)\\.",
    'taxon'        = "^taxon\\.",
    'theme'        = "^theme\\.",
    'theme_bw'     = "^theme_bw\\.",
    'tiplab'       = "^tip(|lab)\\.",
    'trend'        = "^t(|rend)\\.",
    'topo'         = "^topo\\.",
    'vline'        = "^[rv](|line)\\.",
    'violin'       = "^v(|iolin)\\.",
    'xaxis'        = "^x(|axis)\\.",
    'yaxis'        = "^y(|axis)\\.",
    sprintf("^%s\\.", layer) )
  
  
  
  
  # Unprefixed dot arguments, e.g. 'scales'="free_x"
  #________________________________________________________
  layer_fun <- attr(result, 'function', exact = TRUE)
  f_args    <- intersect(
    x = names(layer_params), 
    y = attr(layer_fun, 'formalArgs') %||% formalArgs(layer_fun) )
  for (i in f_args)
    if (!isFALSE(attr(layer_params[[i]], 'display')))
      result[i] <- layer_params[[i]]
  
  
  # Prefixed dot arguments, e.g. 'facet.scales'="free_x"
  #________________________________________________________
  for (i in grep(regex, names(layer_params), value = TRUE))
    if (!isFALSE(attr(layer_params[[i]], 'display')))
      if (!i %in% c('p.top', 'p.label', 'p.adj'))
        result[sub(regex, "", i, perl = TRUE)] <- layer_params[[i]]
  
  
  
  # #________________________________________________________
  # # Synchronize user's facet scales with params$.free_x/y
  # #________________________________________________________
  # if (layer == "facet") {
  #   
  #   scales <- if.null(layers[['facet']][['scales']], 'fixed')
  #   
  #   params$.free_x %<>% if.null(scales %in% c('free', 'free_x'))
  #   params$.free_y %<>% if.null(scales %in% c('free', 'free_y'))
  #   
  #   if (params$.free_x || params$.free_y)
  #     layers[['facet']][['scales']] <- switch(
  #       paste(params$.free_x, params$.free_y),
  #       'TRUE TRUE'  = "free", 
  #       'TRUE FALSE' = "free_x", 
  #       'FALSE TRUE' = "free_y" )
  # }
  
  
  
  
  layers[[layer]] <- result
  
  
  return (invisible(NULL))
}



# #________________________________________________________
# #' Subset `biom` by `*.by` args. Assign colors, etc.
# #' 
# #' @noRd
# #' @keywords internal
# sync_metadata <- function (md, params = parent.frame()) {
#   
#   
#   by_params <- c(
#     "x", "y", "color.by", "shape.by", "facet.by", 
#     "pattern.by", "label.by", "order.by", "stat.by", "limit.by" )
#   
#   
#   
#   #________________________________________________________
#   # Tracks which metadata columns are used.
#   #________________________________________________________
#   md_params <- intersect(by_params, env_names(params))
#   col_names <- c(params[['within']], params[['between']]) # Metadata columns needed in ggdata.
#   revsort   <- c() # For order.by, reverse if '-' prefix.
#   group_by  <- c() # ggdata cols to use for ggplot group.
#   
#   subset_vals <- list()
#   subset_lo   <- list()
#   subset_hi   <- list()
#   
#   
#   #________________________________________________________
#   # Apply subseting defined by 'range' or 'values'.
#   #________________________________________________________
#   for (md_param in md_params) {
#     for (i in seq_along(params[[md_param]])) {
#       
#       col_spec <- params[[md_param]][[i]]
#       col_name <- names(params[[md_param]])[[i]]
#       
#       col_names %<>% c(col_name)
#       revsort   %<>% c(attr(col_spec, 'revsort', exact = TRUE))
#       
#       
#       if (hasName(col_spec, 'range')) {
#         
#         lo <- min(col_spec[['range']])
#         hi <- max(col_spec[['range']])
#         
#         subset_lo[[col_name]] <- max(c(lo, subset_lo[[col_name]]))
#         subset_hi[[col_name]] <- min(c(hi, subset_hi[[col_name]]))
#         if (subset_lo[[col_name]] > subset_hi[[col_name]]) cli_abort("All {col_name} values dropped.")
#       }
#       
#       
#       if (hasName(col_spec, 'values')) {
#         
#         values <- col_spec[['values']]
#         
#         lvls <- names(values) # color.by = list('Sex' = c(Male = "green", Female = "teal"))
#         if (is.null(lvls)) {  # color.by = list('Sex' = c("Male", "Female"))
#           lvls <- as.vector(values)
#           params[[md_param]][[i]][['values']] <- NULL
#         }
#         validate_var_choices( # color.by = list('Sex' = c("m", "f"))
#           var     = 'lvls', 
#           choices = levels(md[[col_name]]), 
#           evar    = paste0(md_param, ": ", col_name), 
#           max     = Inf )
#         
#         if (is.null(subset_vals[[col_name]])) { subset_vals[[col_name]] <- lvls                 }
#         else                                  { subset_vals[[col_name]] %<>% intersect(lvls, .) }
#         if (length(subset_vals[[col_name]]) == 0) cli_abort("All {col_name} values dropped.")
#       }
#       
#     }
#   }
#   
#   
#   #________________________________________________________
#   # Drop rows according to subset specs.
#   #________________________________________________________
#   col_names %<>% unique()
#   md <- md[,intersect(colnames(md), unique(c(".sample", col_names)))]
#   
#   subset_cmds <- c()
#   factor_cmds <- c()
#   omit_cmds   <- c()
#   
#   for (i in col_names) {
#     
#     if (hasName(subset_vals, i) && !identical(levels(md[[i]]), subset_vals[[i]])) {
#       lvls <- subset_vals[[i]]
#       md[[i]] %<>% factor(levels = lvls)
#       lvls_str <- as.args(list(lvls))
#       factor_cmds %<>% c(glue("biom$metadata${coan(i)} %<>% factor(levels = {lvls_str})"))
#     }
#     
#     if (hasName(subset_lo, i) && min(c(Inf, md[[i]]), na.rm = TRUE) < subset_lo[[i]]) {
#       md <- md[md[[i]] >= subset_lo[[i]],]
#       subset_cmds %<>% c(glue("{coan(i)} >= {subset_lo[[i]]}"))
#     }
#     
#     if (hasName(subset_hi, i) && max(c(-Inf, md[[i]]), na.rm = TRUE) > subset_hi[[i]]) {
#       md <- md[md[[i]] <= subset_hi[[i]],]
#       subset_cmds %<>% c(glue("{coan(i)} <= {subset_hi[[i]]}"))
#     }
#     
#   }
#   
#   if (anyNA(md)) {
#     md        <- stats::na.omit(md)
#     omit_cmds <- glue("biom$metadata %<>% na.omit()")
#   }
#   
#   if (length(subset_cmds) > 0)
#     subset_cmds <- glue("biom %<>% subset({paste0(subset_cmds, collapse = ' & ')})")
#   
#   if (length(c(subset_cmds, factor_cmds, omit_cmds)) > 0)
#     params$.subset_code <- paste0(c(subset_cmds, factor_cmds, omit_cmds), collapse = "\n")
#   
#   
#   
#   #________________________________________________________
#   # Discard vestigial factor levels.
#   #________________________________________________________
#   for (i in seq_len(ncol(md)))
#     if (is.factor(md[[i]]))
#       md[[i]] %<>% {factor(., levels = intersect(levels(.), .))}
#   
#   
#   #________________________________________________________
#   # Re-order the samples according to metadata.
#   #________________________________________________________
#   for (i in rev(names(params$order.by)))
#     md <- md[order(md[[i]], decreasing = i %in% revsort),,drop=FALSE]
#   
#   
#   
#   #________________________________________________________
#   # Simplify all *.by params except color/shape/pattern.
#   #________________________________________________________
#   for (param in md_params) {
#     
#     if (length(params[[param]]) == 0) {
#       rlang::env_poke(params, nm = param, value = NULL)
#       
#     } else if (!param %in% c("color.by", "shape.by", "pattern.by", "within", "between")) {
#       params[[param]] %<>% names()
#     }
#   }
#   
#   
#   
#   #________________________________________________________
#   # Explicitly map color/shape/pattern to values.
#   #________________________________________________________
#   for (param in c("color.by", "shape.by", "pattern.by")) {
#     
#     for (i in seq_along(params[[param]])) {
#       
#       col_spec <- params[[param]][[i]]
#       col_name <- attr(col_spec, 'col_name', exact = TRUE)
#       col_type <- attr(col_spec, 'col_type', exact = TRUE)
#       values   <- col_spec[['values']]
#       
#       
#       if (eq(col_type, "cat")) {
#         
#         n <- nlevels(md[[col_name]])
#         
#         
#         # Will need additional colors/shapes/patterns for bdiv box/corr.
#         if (hasName(params, 'within') && n > 1)
#           n <- as.integer(local({
#             if (col_name %in% params$within)  return (n)
#             if (col_name %in% params$between) return ((n * (n - 1) / 2))
#             return ((n * (n - 1) / 2) + n)
#           }))
#         
#         
#         if (is_null(values))
#           values <- switch(
#             EXPR = param,
#             'color.by'   = get_n_colors(n, col_spec[['colors']]),
#             'shape.by'   = get_n_shapes(n),
#             'pattern.by' = get_n_patterns(n) )
#         
#         
#       } else if (eq(col_type, "num")) {
#         if (hasName(col_spec, 'colors')) values <- col_spec[['colors']]
#         if (is_palette(values))          values %<>% get_palette()
#       }
#       
#       col_spec[['colors']] <- NULL
#       col_spec[['values']] <- values
#       params[[param]][[i]] <- col_spec
#     }
#   }
#   
#   
#   
#   return (md)
# }



