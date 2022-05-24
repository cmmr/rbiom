

#________________________________________________________
# Finds parent's `layers` variable for a particular layer
#________________________________________________________
hasLayer <- function (layer=get("layer", pos = parent.frame())) {
  hasName(get("layers", pos = parent.frame()), layer)
}


#________________________________________________________
# Finds `layers` variable in parent env, then updates it.
#________________________________________________________
setLayer <- function (layer=get("layer", pos = parent.frame()), ...) {
  
  
  layers <- get("layers", pos = parent.frame())
  
  if (!is(layer, 'character')) browser()
  if (!is(layers, 'list')) browser()
  
  if (!hasName(layers, layer))
    initLayer(layer)
  
  dots <- list(...)
  
  
  # Expand lists. For cases of e.g.
  # setLayer('box', fill="white", list(size=2))
  keyvals <- list()
  for (i in seq_along(dots)) {
    if (isTRUE(nzchar(key <- names(dots)[[i]]))) {
      if (is.null(dots[[i]])) {
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
        if (is.null(val)) {
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
        if (is.null(val)) {
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
#________________________________________________________
initLayer <- function (layer_names) {
  
  layers <- get("layers", pos = parent.frame())
  
  params <- attr(layers, 'params', exact = TRUE)
  xmode  <- attr(layers, 'xmode',  exact = TRUE)
  
  patterned <- !is.null(params[['pattern.by']])
  facetDims <- length(params[['facet.by']])
  
  
  for (layer in layer_names) {
  
    if (hasName(layers, layer))
      stop("Cannot initialize layer '", layer, "' twice.")
  
    result <- list()
    
    attr(result, 'fn') <- switch(
      EXPR = layer,
      'ggplot'     = "ggplot",
      'arrow'      = "geom_segment",
      'bar'        = ifelse(patterned, "ggpattern::geom_bar_pattern", "geom_bar"),
      'box'        = ifelse(patterned, "ggpattern::geom_boxplot_pattern", "geom_boxplot"),
      'brackets'   = "geom_segment",
      'color'      = ifelse(patterned, "ggpattern::scale_pattern_color_manual", "scale_color_manual"),
      'crossbar'   = ifelse(patterned, "ggpattern::geom_crossbar_pattern", "geom_crossbar"),
      'density'    = "ggdensity::geom_hdr",
      'dot'        = "ggbeeswarm::geom_beeswarm",
      'errorbar'   = "geom_errorbar",
      'ellipse'    = "stat_ellipse",
      'facet'      = ifelse(facetDims == 2, "facet_grid", "facet_wrap"),
      'fill'       = ifelse(patterned, "ggpattern::scale_pattern_fill_manual", "scale_fill_manual"),
      'labs'       = "labs",
      'linerange'  = "geom_linerange",
      'mean'       = "geom_point",
      'name'       = "geom_text",
      'pattern'    = "ggpattern::scale_pattern_type_manual",
      'pointrange' = "geom_pointrange",
      'point'      = "geom_point",
      'regression' = "stat_smooth",
      'shape'      = "scale_shape_manual",
      'spider'     = "geom_segment",
      'stack'      = ifelse(patterned, "ggpattern::geom_col_pattern", "geom_col"),
      'stats_text' = "geom_text",
      'strip'      = "ggbeeswarm::geom_quasirandom",
      'taxon'      = "ggrepel::geom_label_repel",
      'theme'      = "theme",
      'theme_bw'   = "theme_bw",
      'topo'       = "ggdensity::geom_hdr_lines",
      'violin'     = ifelse(patterned, "ggpattern::geom_violin_pattern", "geom_violin"),
      'xaxis'      = ifelse(xmode == "factor", "scale_x_discrete", "scale_x_continuous"),
      'yaxis'      = "scale_y_continuous",
      layer )
    
    attr(result, 'function') <- eval(parse(text = attr(result, 'fn', exact = TRUE)))
    
    
    attr(result, 'src') <- switch(
      EXPR = layer,
      'arrow'      = "biplot",
      'brackets'   = "stat_brackets",
      'crossbar'   = "vline",
      'errorbar'   = "vline",
      'linerange'  = "vline",
      'mean'       = "biplot",
      'pointrange' = "vline",
      'spider'     = "spider",
      'stats_text' = "stat_labels",
      'taxon'      = "biplot" )
    
    
    result[['data']] <- switch(
      EXPR = layer,
      'arrow'      = ~ attr(., "biplot"),
      'brackets'   = ~ attr(., "stat_brackets"),
      'crossbar'   = ~ attr(., "vline"),
      'errorbar'   = ~ attr(., "vline"),
      'linerange'  = ~ attr(., "vline"),
      'mean'       = ~ attr(., "biplot"),
      'pointrange' = ~ attr(., "vline"),
      'spider'     = ~ attr(., "spider"),
      'stats_text' = ~ attr(., "stat_labels"),
      'taxon'      = ~ attr(., "biplot") )
    
    
    regex <- switch(
      EXPR = layer,
      'ggplot'     = "^g(|gplot)\\.",
      'arrow'      = "^a(|rrow)\\.",
      'bar'        = "^(|ba)r\\.",
      'box'        = "^b(|ox)\\.",
      'brackets'   = "^h(|line)\\.",
      'color'      = "^color\\.",
      'crossbar'   = "^c(|rossbar)\\.",
      'density'    = "^d(|ensity)\\.",
      'dot'        = "^d(|ot)\\.",
      'errorbar'   = "^e(|rrorbar)\\.",
      'ellipse'    = "^e(|llipse)\\.",
      'facet'      = "^f(|acet)\\.",
      'fill'       = "^fill\\.",
      'labs'       = "^labs\\.",
      'linerange'  = "^l(|inerange)\\.",
      'mean'       = "^m(|ean)\\.",
      'name'       = "^n(|ame)\\.",
      'pattern'    = "^pattern\\.",
      'pointrange' = "^p(|ointrange)\\.",
      'point'      = "^p(|oint)\\.",
      'regression' = "^r(|egression)\\.",
      'shape'      = "^shape\\.",
      'size'       = "^size\\.",
      'spider'     = "^s(|pider)\\.",
      'stack'      = "^s(|tack)\\.",
      'stats_text' = "^a(|nno)\\.",
      'strip'      = "^s(|trip)\\.",
      'taxon'      = "^taxon\\.",
      'topo'       = "^topo\\.",
      'theme'      = "^t(|heme)\\.",
      'theme_bw'   = "^theme_bw\\.",
      'violin'     = "^v(|iolin)\\.",
      'xaxis'      = "^x(|axis)\\.",
      'yaxis'      = "^y(|axis)\\.",
      sprintf("^%s\\.", layer) )
    
    
    
    params <- attr(layers, 'params', exact = TRUE)
    func   <- attr(result, 'function', exact = TRUE)
    
    
    
    # Blacklist named parameters for rbiom's plot() function
    #________________________________________________________
    params <- params[setdiff(names(params), formalArgs(plot.BIOM))]
    
    
    # Ensure NULLs can be inserted into a list
    #________________________________________________________
    for (i in seq_along(params))
      if (is.null(params[[i]]))
        params[[i]] <- list(NULL)
    
    
    # Arrays need a list wrapper too. plot(yaxis.expand=1:4)
    #________________________________________________________
    for (i in seq_along(params))
      if (length(params[[i]]) > 1)
        params[[i]] <- list(params[[i]])
    
    
    # Unprefixed dot arguments, e.g. 'scales'="free_x"
    #________________________________________________________
    for (i in intersect(names(params), formalArgs(func)))
      result[i] <- params[[i]]
    
    
    # Prefixed dot arguments, e.g. 'facet.scales'="free_x"
    #________________________________________________________
    for (i in grep(regex, names(params), value = TRUE))
      result[sub(regex, "", i, perl = TRUE)] <- params[[i]]
    
    
    layers[[layer]] <- result
  }
  
  assign("layers", layers, pos = parent.frame())
  return (invisible(NULL))
}

