

#________________________________________________________
# Finds parent's `layers` variable for a particular layer
#________________________________________________________
hasLayer <- function (layer=get("layer", pos = parent.frame())) {
  hasName(get("layers", pos = parent.frame()), layer)
}


#________________________________________________________
# Finds `layers` variable in parent env, then updates it.
#________________________________________________________
setLayer <- function (layer=get("layer", pos = parent.frame()), ..., .fn = NULL) {
  
  
  layers <- get("layers", pos = parent.frame())
  
  if (!is(layer, 'character')) browser()
  if (!is(layers, 'list')) browser()
  
  if (!hasName(layers, layer))
    initLayer(layer, fn = .fn)
  
  dots <- list(...)
  
  
  # Expand lists. For cases of e.g.
  # setLayer('box', fill="white", list(size=2))
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
    
    # setLayer("bar", fill = "green")
    if (length(key) == 1) {
      if (!hasName(layers[[layer]], key))
        if (is_null(val)) {
          layers[[layer]][key] <- list(NULL)
        } else {
          layers[[layer]][[key]] <- val
        }
      
    # setLayer("bar", "mapping|fill" = "Body Site")
    } else if (length(key) == 2) {
      key1 <- key[[1]]
      key2 <- key[[2]]
      
      if (!hasName(layers[[layer]], key1))
        layers[[layer]][[key1]] <- list()
        
      if (!hasName(layers[[layer]][[key1]], key2))
        if (is_null(val)) {
          layers[[layer]][[key1]][key2] <- list(NULL)
        } else {
          layers[[layer]][[key1]][[key2]] <- val
        }
      
    } else {
      stop("Key vector is too long. Key '", key, "' split on '|'.")
    }
    
  }
  
  
  assign("layers", layers, pos = parent.frame())
  
  return (invisible(NULL))
}


#________________________________________________________
# Set up a layer with its function, function name, etc.
# Uses attr(layers, 'params') to add in user's params.
# Can set the underlying function using 'fn' parameter.
#________________________________________________________
initLayer <- function (layer_names, fn = NULL) {
  
  layers <- get("layers", pos = parent.frame())
  
  params <- attr(layers, 'params', exact = TRUE)
  xmode  <- attr(layers, 'xmode',  exact = TRUE)
  
  patterned <- !is_null(params[['pattern.by']])
  facetDims <- length(params[['facet.by']])
  
  
  
  for (layer in layer_names) {
  
    if (hasName(layers, layer))
      stop("Cannot initialize layer '", layer, "' twice.")
  
    result <- list()
    params <- attr(layers, 'params', exact = TRUE)
    
    
    #________________________________________________________
    # color.by, pattern.by, and shape.by arguments.
    #________________________________________________________
    
    if (layer %in% c('shape', 'pattern'))
      result <- as.list(params[[paste0(layer, ".by")]][[1]])
    
    if (layer %in% c('color', 'fill'))
      result <- as.list(params[["color.by"]][[1]])
    
    
    attr(result, 'function') <- switch(
      EXPR = layer,
      'ggplot'       = ggplot,
      'ggtree'       = ggtree,
      'arrow'        = geom_segment,
      'bar'          = if (patterned) geom_bar_patternelse else geom_bar,
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
      'residual'     = geom_segment,
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
      'trend'        = stat_smooth,
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
      'residual'     = "residual",
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
      'residual'     = ~ attr(., "residual"),
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
      'residual'     = "^r(|esidual)\\.",
      'shape'        = "^shape\\.",
      'scatter'      = "^(pt|s|scatter)\\.",
      'size'         = "^size\\.",
      'smooth'       = "^s(|mooth)\\.",
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
      'trend'        = "^(t|trend|c|confidence)\\.",
      'vline'        = "^[rv](|line)\\.",
      'violin'       = "^v(|iolin)\\.",
      'xaxis'        = "^x(|axis)\\.",
      'yaxis'        = "^y(|axis)\\.",
      sprintf("^%s\\.", layer) )
    
    
    
    layer_func <- attr(result, 'function', exact = TRUE)
    params     <- attr(layers, 'params',   exact = TRUE)
    
    
    
    # Ignore formal arguments in the called plotting function
    #________________________________________________________
    plot_func <- attr(layers, 'function', exact = TRUE)
    params    <- params[setdiff(names(params), formalArgs(plot_func))]
    
    
    # Ensure NULLs can be inserted into a list
    #________________________________________________________
    for (i in seq_along(params))
      if (is_null(params[[i]]))
        params[[i]] <- list(NULL)
    
    
    # Arrays and functions need list wrappers too.
    # plot(yaxis.expand=1:4, yaxis.trans=sqrt)
    #________________________________________________________
    for (i in seq_along(params))
      if (length(params[[i]]) > 1 || is.function(params[[i]]))
        params[[i]] <- list(params[[i]])
    
    
    # Unprefixed dot arguments, e.g. 'scales'="free_x"
    #________________________________________________________
    for (i in intersect(names(params), formalArgs(layer_func)))
      if (!isFALSE(attr(params[[i]], 'display')))
        result[i] <- params[[i]]
    
    # Prefixed dot arguments, e.g. 'facet.scales'="free_x"
    #________________________________________________________
    for (i in grep(regex, names(params), value = TRUE))
      if (!isFALSE(attr(params[[i]], 'display')))
        result[sub(regex, "", i, perl = TRUE)] <- params[[i]]
    
    
    layers[[layer]] <- result
  }
  
  assign("layers", layers, pos = parent.frame())
  return (invisible(NULL))
}

