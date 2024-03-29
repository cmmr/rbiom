#________________________________________________________
# Converts user's `layers` spec to `layers` environment.
#________________________________________________________
init_layers <- function (params = parent.frame(), choices = NULL, var = "layers", no_init = NULL, do_init = NULL) {
  
  stopifnot(env_has(params, c('.ggdata', '.xcol', '.ycol', '.xmode')))
  
  
  layer_names <- do_init
  
  
  if (!is.null(choices)) {
    
    stopifnot(is_scalar_character(var) && env_has(params, var))
    spec <- get(var, pos = params, inherits = FALSE)
    if (!is_character(spec) || anyNA(spec))
      stop("Invalid '", var, "' value.")
    
    
    for (i in tolower(spec[nchar(spec) > 0]))
      layer_names <- c(layer_names, local({
        
        if (i %in% choices)        return (i)
        if (i %in% names(choices)) return (choices[[i]])
        
        pm <- pmatch(i, choices)
        ii <- strsplit(i, '')[[1]]
        
        if (length(spec) > 1 && !is.na(pm)) return (choices[[pm]])
        if (all(ii %in% names(choices)))    return (choices[ii])
        if (!is.na(pm))                     return (choices[[pm]])
        return (choices[intersect(ii, names(choices))])
      }))
    
    layer_names <- layer_names[!is.na(layer_names)]
    
    
    
    
    #________________________________________________________
    # Ignore shapes/etc without applicable layers.
    #________________________________________________________
    if (!any(c('dot', 'strip', 'pointrange', 'scatter', 'point') %in% layer_names)) 
       params$shape.by <- NULL
    
    if (!any(c('box', 'bar', 'violin')   %in% layer_names)) params$pattern.by <- NULL
    if (!any(c('taxon', 'arrow', 'mean') %in% layer_names)) params$rank <- NULL
    
    
    if (length(layer_names) == 0)
      stop("Invalid 'layers' argument.")
    
    
    #________________________________________________________
    # Add in aesthetic layers as needed.
    #________________________________________________________
    if (!is_null(params$color.by))   layer_names %<>% c('color')
    if (!is_null(params$shape.by))   layer_names %<>% c('shape')
    if (!is_null(params$pattern.by)) layer_names %<>% c('pattern')
    if (!is_null(params$facet.by))   layer_names %<>% c('facet')
    
  }
  
  
  
  init_layer_names <- layer_names %>%
    c('ggplot', ., 'xaxis', 'yaxis', 'labs', 'theme', 'theme_bw') %>%
    setdiff(., no_init) %>% 
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
  hasName(params$layers, layer)
}



#________________________________________________________
# Deletes an entire layer.
#________________________________________________________
del_layer <- function (params, layer) {
  
  if (hasName(params$layers, layer))
    remove(list = layer, pos = layers)
  
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
  # set_layer(params, 'box', fill="white", list(size=2))
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
      stop("Key vector is too long. Key '", key, "' split on '|'.")
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
    stop("Cannot add layer '", layer, "' twice.")
  
  
  xmode     <- params$.xmode
  patterned <- !is_null(params$pattern.by)
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
  # plot(yaxis.expand=1:4, yaxis.trans=sqrt)
  #________________________________________________________
  for (i in seq_along(layer_params))
    if (length(layer_params[[i]]) > 1 || is.function(layer_params[[i]]))
      layer_params[[i]] <- list(layer_params[[i]])
  
  
  
  
  
  #________________________________________________________
  # color.by, pattern.by, and shape.by arguments.
  #________________________________________________________
  
  result <- list()
  
  if (layer %in% c('shape', 'pattern'))
    result <- as.list(params[[paste0(layer, ".by")]][[1]])
  
  if (layer %in% c('color', 'fill')) {
    result <- as.list(params[["color.by"]][[1]])
  }
  
  
  
  attr(result, 'function') <- switch(
    EXPR = layer,
    'ggplot'       = ggplot,
    'ggtree'       = ggtree,
    'arrow'        = geom_segment,
    'bar'          = if (patterned) geom_bar_pattern else geom_bar,
    'box'          = if (patterned) geom_boxplot_pattern else geom_boxplot,
    'brackets'     = geom_segment,
    'cladelab'     = geom_cladelab,
    'color'        = if (patterned) scale_pattern_color_manual else scale_color_manual,
    'crossbar'     = if (patterned) geom_crossbar_pattern else geom_crossbar,
    'density'      = geom_hdr,
    'dot'          = geom_beeswarm,
    'errorbar'     = geom_errorbar,
    'ellipse'      = stat_ellipse,
    'facet'        = if (facetDims == 2) facet_grid else facet_wrap,
    'fill'         = if (patterned) scale_pattern_fill_manual else scale_fill_manual,
    'flip'         = coord_flip,
    # 'free_y'     = facetted_pos_scales,
    'hline'        = geom_hline,
    'label'        = geom_label,
    'labs'         = labs,
    'linerange'    = geom_linerange,
    'mean'         = geom_point,
    'name'         = geom_text,
    'pattern'      = scale_pattern_type_manual,
    'pointrange'   = geom_pointrange,
    'point'        = geom_point,
    'rect'         = geom_rect,
    'scatter'      = geom_point,
    'shape'        = scale_shape_manual,
    'smooth'       = stat_smooth,
    'spider'       = geom_segment,
    'stack'        = if (patterned) geom_col_pattern else geom_col,
    'stats_bg'     = geom_rect,
    'stats_text'   = geom_text,
    'stats_ggtext' = geom_textbox,
    'strip'        = geom_quasirandom,
    'stripe'       = geom_rect,
    'taxon'        = geom_label_repel,
    'theme'        = theme,
    'theme_bw'     = theme_bw,
    'tiplab'       = geom_tiplab,
    'topo'         = geom_hdr_lines,
    'violin'       = if (patterned) geom_violin_pattern else geom_violin,
    'vline'        = geom_vline,
    'xaxis'        = if (xmode == "factor") scale_x_discrete else scale_x_continuous,
    'yaxis'        = scale_y_continuous,
    if (is.null(fn)) get(layer) else get(fn) )
  
  
  attr(result, 'src') <- switch(
    EXPR = layer,
    'arrow'        = "taxa_coords",
    'brackets'     = "stat_brackets",
    'crossbar'     = "vline",
    'errorbar'     = "vline",
    'linerange'    = "vline",
    'mean'         = "taxa_coords",
    'pointrange'   = "vline",
    'spider'       = "spider",
    'stats_text'   = "stat_labels",
    'stats_ggtext' = "stat_labels",
    'taxon'        = "taxa_coords" )
  
  
  result[['data']] <- switch(
    EXPR = layer,
    'arrow'        = ~ attr(., "taxa_coords"),
    'brackets'     = ~ attr(., "stat_brackets"),
    'crossbar'     = ~ attr(., "vline"),
    'errorbar'     = ~ attr(., "vline"),
    'linerange'    = ~ attr(., "vline"),
    'mean'         = ~ attr(., "taxa_coords"),
    'pointrange'   = ~ attr(., "vline"),
    'spider'       = ~ attr(., "spider"),
    'stats_text'   = ~ attr(., "stat_labels"),
    'stats_ggtext' = ~ attr(., "stat_labels"),
    'taxon'        = ~ attr(., "taxa_coords") )
  
  
  regex <- switch(
    EXPR = layer,
    'ggplot'       = "^g(|gplot)\\.",
    'ggtree'       = "^g(|gtree)\\.",
    'arrow'        = "^a(|rrow)\\.",
    'bar'          = "^(r|bar)\\.",
    'box'          = "^b(|ox)\\.",
    'brackets'     = "^h(|line)\\.",
    'cladelab'     = "^c(|lade|ladelab)\\.",
    'color'        = "^color\\.",
    'crossbar'     = "^c(|rossbar)\\.",
    'density'      = "^d(|ensity)\\.",
    'dot'          = "^(pt|d|dot)\\.",
    'errorbar'     = "^e(|rrorbar)\\.",
    'ellipse'      = "^e(|llipse)\\.",
    'facet'        = "^f(|acet)\\.",
    'fill'         = "^fill\\.",
    'flip'         = "^flip\\.",
    'heat'         = "^h(|eat|eatmap)\\.",
    'hexpand'      = "^h(|expand)\\.",
    'hline'        = "^[rh](|line)\\.",
    'label'        = "^l(|abel)\\.",
    'labs'         = "^labs\\.",
    'linerange'    = "^l(|inerange)\\.",
    'mean'         = "^m(|ean)\\.",
    'name'         = "^n(|ame)\\.",
    'pattern'      = "^pattern\\.",
    'pointrange'   = "^(pt|p|pointrange)\\.",
    'point'        = "^(pt|p|point)\\.",
    'shape'        = "^shape\\.",
    'scatter'      = "^(pt|s|scatter)\\.",
    'size'         = "^size\\.",
    'smooth'       = "^(t|c|trend|confidence)\\.",
    'spider'       = "^s(|pider)\\.",
    'stack'        = "^s(|tack)\\.",
    'stats_text'   = "^stats\\.",
    'stats_ggtext' = "^stats\\.",
    'strip'        = "^(pt|s|strip)\\.",
    'taxon'        = "^taxon\\.",
    'theme'        = "^theme\\.",
    'theme_bw'     = "^theme_bw\\.",
    'tiplab'       = "^tip(|lab)\\.",
    'topo'         = "^topo\\.",
    'vline'        = "^[rv](|line)\\.",
    'violin'       = "^v(|iolin)\\.",
    'xaxis'        = "^x(|axis)\\.",
    'yaxis'        = "^y(|axis)\\.",
    sprintf("^%s\\.", layer) )
  
  
  
  
  # Unprefixed dot arguments, e.g. 'scales'="free_x"
  #________________________________________________________
  layer_fun <- attr(result, 'function', exact = TRUE)
  for (i in intersect(names(layer_params), formalArgs(layer_fun)))
    if (!isFALSE(attr(layer_params[[i]], 'display')))
      result[i] <- layer_params[[i]]
  
  # Prefixed dot arguments, e.g. 'facet.scales'="free_x"
  #________________________________________________________
  for (i in grep(regex, names(layer_params), value = TRUE))
    if (!isFALSE(attr(layer_params[[i]], 'display')))
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



#________________________________________________________
#' Subset `biom` by `*.by` args. Assign colors, etc.
#' 
#' @noRd
#' @keywords internal
sync_metadata <- function (params = parent.frame()) {
  
  
  md <- params$biom$metadata
  
  
  by_params <- c(
    "x", "color.by", "shape.by", "facet.by", 
    "pattern.by", "label.by", "order.by", "stat.by", "limit.by" )
  
  
  
  #________________________________________________________
  # Tracks which metadata columns are used.
  #________________________________________________________
  md_params <- intersect(by_params, env_names(params))
  col_names <- c(params[['within']], params[['between']]) # Metadata columns needed in ggdata.
  revsort   <- c() # For order.by, reverse if '-' prefix.
  group_by  <- c() # ggdata cols to use for ggplot group.
  
  subset_vals <- list()
  subset_lo   <- list()
  subset_hi   <- list()
  
  
  #________________________________________________________
  # Apply subseting defined by 'range' or 'values'.
  #________________________________________________________
  for (md_param in md_params) {
    for (i in seq_along(params[[md_param]])) {
      
      col_spec <- params[[md_param]][[i]]
      col_name <- names(params[[md_param]])[[i]]
      
      col_names %<>% c(col_name)
      revsort   %<>% c(attr(col_spec, 'revsort', exact = TRUE))
      
      
      if (hasName(col_spec, 'range')) {
        
        lo <- min(col_spec[['range']])
        hi <- max(col_spec[['range']])
        
        subset_lo[[col_name]] <- max(c(lo, subset_lo[[col_name]]))
        subset_hi[[col_name]] <- min(c(hi, subset_hi[[col_name]]))
        if (subset_lo[[col_name]] > subset_hi[[col_name]]) cli_abort("All {col_name} values dropped.")
      }
      
      
      if (hasName(col_spec, 'values')) {
        
        values <- col_spec[['values']]
        
        lvls <- names(values) # color.by = list('Sex' = c(Male = "green", Female = "teal"))
        if (is.null(lvls)) {  # color.by = list('Sex' = c("Male", "Female"))
          lvls <- as.vector(values)
          params[[md_param]][[i]][['values']] <- NULL
        }
        validate_var_choices( # color.by = list('Sex' = c("m", "f"))
          var     = 'lvls', 
          choices = levels(params$biom$metadata[[col_name]]), 
          evar    = paste0(md_param, ": ", col_name), 
          max     = Inf )
        
        if (is.null(subset_vals[[col_name]])) { subset_vals[[col_name]] <- lvls                 }
        else                                  { subset_vals[[col_name]] %<>% intersect(lvls, .) }
        if (length(subset_vals[[col_name]]) == 0) cli_abort("All {col_name} values dropped.")
      }
      
    }
  }
  
  
  #________________________________________________________
  # Drop rows according to subset specs.
  #________________________________________________________
  col_names %<>% unique()
  md <- md[,unique(c(".sample", col_names))]
  
  subset_cmds <- c()
  factor_cmds <- c()
  omit_cmds   <- c()
  
  for (i in col_names) {
    
    if (hasName(subset_vals, i) && !identical(levels(md[[i]]), subset_vals[[i]])) {
      lvls <- subset_vals[[i]]
      md[[i]] %<>% factor(levels = lvls)
      lvls_str <- as.args(list(lvls))
      factor_cmds %<>% c(glue("biom$metadata${coan(i)} %<>% factor(levels = {lvls_str})"))
    }
    
    if (hasName(subset_lo, i) && min(c(Inf, md[[i]]), na.rm = TRUE) < subset_lo[[i]]) {
      md <- md[md[[i]] >= subset_lo[[i]],]
      subset_cmds %<>% c(glue("{coan(i)} >= {subset_lo[[i]]}"))
    }
    
    if (hasName(subset_hi, i) && max(c(-Inf, md[[i]]), na.rm = TRUE) > subset_hi[[i]]) {
      md <- md[md[[i]] <= subset_hi[[i]],]
      subset_cmds %<>% c(glue("{coan(i)} <= {subset_hi[[i]]}"))
    }
    
  }
  
  if (anyNA(md)) {
    md        <- stats::na.omit(md)
    omit_cmds <- glue("biom$metadata %<>% na.omit()")
  }
  
  if (length(subset_cmds) > 0)
    subset_cmds <- glue("biom %<>% subset({paste0(subset_cmds, collapse = ' & ')})")
  
  if (length(c(subset_cmds, factor_cmds, omit_cmds)) > 0)
    params$.subset_code <- paste0(c(subset_cmds, factor_cmds, omit_cmds), collapse = "\n")
  
  
  
  #________________________________________________________
  # Discard vestigial factor levels.
  #________________________________________________________
  for (i in seq_len(ncol(md)))
    if (is.factor(md[[i]]))
      md[[i]] %<>% {factor(., levels = intersect(levels(.), .))}
  
  
  #________________________________________________________
  # Re-order the samples according to metadata.
  #________________________________________________________
  for (i in rev(names(params$order.by)))
    md <- md[order(md[[i]], decreasing = i %in% revsort),,drop=FALSE]
  
  
  
  #________________________________________________________
  # Simplify all *.by params except color/shape/pattern.
  #________________________________________________________
  for (param in md_params) {
    
    if (length(params[[param]]) == 0) {
      rlang::env_poke(params, nm = param, value = NULL)
      
    } else if (!param %in% c("color.by", "shape.by", "pattern.by", "within", "between")) {
      params[[param]] %<>% names()
    }
  }
  
  
  
  #________________________________________________________
  # Explicitly map color/shape/pattern to values.
  #________________________________________________________
  for (param in c("color.by", "shape.by", "pattern.by")) {
    
    for (i in seq_along(params[[param]])) {
      
      col_spec <- params[[param]][[i]]
      col_name <- attr(col_spec, 'col_name', exact = TRUE)
      col_type <- attr(col_spec, 'col_type', exact = TRUE)
      values   <- col_spec[['values']]
      
      
      if (eq(col_type, "cat")) {
        
        n <- nlevels(md[[col_name]])
        
        
        # Will need additional colors/shapes/patterns for bdiv box/corr.
        if (hasName(params, 'within') && n > 1)
          n <- as.integer(local({
            if (col_name %in% params$within)  return (n)
            if (col_name %in% params$between) return ((n * (n - 1) / 2))
            return ((n * (n - 1) / 2) + n)
          }))
        
        
        if (is_null(values))
          values <- switch(
            EXPR = param,
            'color.by'   = get_n_colors(n, col_spec[['colors']]),
            'shape.by'   = get_n_shapes(n),
            'pattern.by' = get_n_patterns(n) )
        
        
      } else if (eq(col_type, "num")) {
        if (hasName(col_spec, 'colors')) values <- col_spec[['colors']]
        if (is_palette(values))          values %<>% get_palette()
      }
      
      col_spec[['colors']] <- NULL
      col_spec[['values']] <- values
      params[[param]][[i]] <- col_spec
    }
  }
  
  
  
  
  # `params` and `biom` are both environments
  params$biom <- params$biom$clone()
  params$biom$metadata <- md
  
  
  return (invisible(NULL))
}



