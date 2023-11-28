
#________________________________________________________
# Computes p-values for categorical differences
#________________________________________________________
boxplot_stats <- function (params) {
  
  stopifnot(is_bare_environment(params))
  
  
  layers <- params$layers
  ggdata <- params$.ggdata
  xcol   <- params$.xcol
  ycol   <- params$.ycol
  
  
  color.by    <- names(params$color.by)
  shape.by    <- names(params$shape.by)
  pattern.by  <- names(params$pattern.by)
  facet.by    <- params$facet.by
  p.label     <- params$p.label
  
  subgroups <- c(
    names(params$color.by),
    names(params$shape.by),
    names(params$pattern.by) ) %>%
    unique() %>%
    setdiff(c(xcol, facet.by))
  
  
  
  #________________________________________________________
  # No work to do.
  #________________________________________________________
  if (plyr::empty(ggdata) || all(p.label == 0)) {
    set_layer(params, 'yaxis', expand = c(0.02, 0, 0.02, 0))
    return (invisible(params))
  }
  
  
  
  #________________________________________________________
  # Pairwise requires a valid `x`, and no sub-groupings.
  #________________________________________________________
  if (!xcol %in% c(".all", ".taxa") && length(subgroups) == 0) {
    test  <- "pw_means"
    stats <- stats_table(
      df       = ggdata, 
      stat.by  = xcol, 
      resp     = ycol,
      test     = test, 
      level    = params$level, 
      split.by = facet.by, 
      p.adj    = params$p.adj )
  } 
  
  #________________________________________________________
  # Groupwise requires a sub-grouping.
  #________________________________________________________
  else if (length(subgroups) == 1) {
    test  <- "means"
    stats <- stats_table(
      df       = ggdata, 
      stat.by  = subgroups, 
      resp     = ycol,
      test     = test, 
      level    = params$level, 
      split.by = c(xcol, facet.by) %>% setdiff(".all"), 
      p.adj    = params$p.adj )
  }
  
  #________________________________________________________
  # Plot definition isn't compatible with stats.
  #________________________________________________________
  else {
    set_layer(params, 'yaxis', expand = c(0.02, 0, 0.02, 0))
    return (invisible(params))
  }
  
  params$.plot_attrs[['stats']] <- stats
  
  
  
  
  #________________________________________________________
  # Implement `p.top` for taxa_boxplot().
  #________________________________________________________
  if (isTRUE(is.finite(params$p.top))) {
    
    apply_p.top(params)
    
    if (plyr::empty(params$.ggdata))
      return (invisible(params))
  }
  
  
  #________________________________________________________
  # Drop p-values that are NA or below the display threshold.
  #________________________________________________________
  if (hasName(stats, '.p.val')) stats <- stats[!is.na(stats[['.p.val']]),]
  if (hasName(stats, '.adj.p')) stats <- stats[stats[['.adj.p']] <= max(p.label),]
  
  
  #________________________________________________________
  # No p-values to display.
  #________________________________________________________
  if (plyr::empty(stats)) {
    set_layer(params, 'yaxis', expand = c(0.02, 0, 0.02, 0))
    return (invisible(params))
  }
  
  
  
  #________________________________________________________
  # The y-position needed to not overlap layer content.
  #________________________________________________________
  
  y_pos_fun <- local({
    
    y_pos_funs <- c(
        if (has_layer(params, 'violin'))     "violin" else NULL,
        if (has_layer(params, 'dot'))        "max"    else NULL,
        if (has_layer(params, 'strip'))      "max"    else NULL,
        if (has_layer(params, 'box'))        "box"    else NULL,
        if (has_layer(params, 'errorbar'))   "vline"  else NULL,
        if (has_layer(params, 'pointrange')) "vline"  else NULL,
        if (has_layer(params, 'crossbar'))   "vline"  else NULL,
        if (has_layer(params, 'linerange'))  "vline"  else NULL,
        if (has_layer(params, 'bar'))        "mean"   else NULL ) %>% 
      unique() %>%
      lapply(switch,
        max    = {function (v) { max(v)                    }},
        mean   = {function (v) { mean(v)                   }},
        box    = {function (v) { boxplot.stats(v)$stats[5] }},
        violin = {function (v) { max(density(v)[['x']])    }},
        vline  = switch(params$ci, 
          range = {function (v) { max(v)                                                }},
          mad   = {function (v) { median(v) + mad(v, median(v))                         }},
          sd    = {function (v) { mean(v) + sd(v)                                       }},
          se    = {function (v) { mean(v) + sqrt(var(v)/length(v))                      }},
          ci    = {function (v) { t.test(v, conf.level = params$level)$conf.int[2] }} ))
    
    function (df) {
      
      vals <- df[[ycol]]
      
      y_pos_funs %>%
        lapply(function (f) tryCatch(f(vals), error = function (e) NA)) %>%
        unlist() %>%
        max(0, na.rm = TRUE) %>%
        labeling::extended(dmin = 0, dmax = ., m = 5, only.loose = TRUE) %>%
        max(0, na.rm = TRUE)
    }
    
  })
  
  
  
  #________________________________________________________
  # Minimum positions for brackets/labels, per facet.
  #________________________________________________________
  if (is.null(facet.by) || isFALSE(params$.free_y)) {
    
    stats[['.ypos']] <- ggdata %>%
      plyr::ddply(
        .variables = ply_cols(c(xcol, subgroups, facet.by)), 
        .fun       = function (df) data.frame(.ypos = y_pos_fun(df)) ) %>%
      dplyr::pull('.ypos') %>%
      max(0, na.rm = TRUE)
    
  } else {
    
    stats %<>% dplyr::left_join(
      by = facet.by, 
      y  = plyr::ddply(ggdata, ply_cols(facet.by), function (df) {
        plyr::daply(df, ply_cols(c(xcol, subgroups)), y_pos_fun) %>%
        max(0, na.rm = TRUE) %>%
        data.frame(.ypos = .) }) )
  }
  
  
  
  
  
  #________________________________________________________
  # Format how p-values are displayed.
  #________________________________________________________
  if (length(p.label) == 1) {
    
    stats[['.label']] <- stats[['.adj.p']] %>%
      formatC(format="g", digits = 1) %>%
      paste("italic(p)==", .)
    
  } else {
    for (i in p.label)
      stats[['.label']] %<>% paste0(ifelse(stats[['.adj.p']] <= i, "*", ""))
  }
  
  
  
  #________________________________________________________
  # Pairwise Brackets ====
  #________________________________________________________
  if (test == "pw_means") {
    
    
    # Convert e.g. 'Saliva - Stool' to 'Saliva' and 'Stool'
    #________________________________________________________
    pairs <- attr(params$.plot_attrs[['stats']], 'pairs', exact = TRUE)
    stats[['.cat1']] <- sapply(stats[[xcol]], function (p) pairs[[p]][[1]])
    stats[['.cat2']] <- sapply(stats[[xcol]], function (p) pairs[[p]][[2]])
    remove("pairs")
    
    
    # Convert e.g. c('Saliva', 'Stool') to c(4, 5)
    #________________________________________________________
    xpos  <- as.factor(ggdata[[xcol]])
    xpos  <- setNames(seq_along(levels(xpos)), levels(xpos))
    stats[['.xmin']]  <- unname(xpos[stats[['.cat1']]])
    stats[['.xmax']]  <- unname(xpos[stats[['.cat2']]])
    
    
    if (is.null(facet.by)) {
      stats[['.step']] <- stats[['.ypos']] * .13
      stats[['.ypos']] <- stats[['.ypos']] + (stats[['.step']] * seq_len(nrow(stats)))
      
    } else {
      
      # When some x-values are absent from some facets, we'll need
      # to re-number the remaining x-positions.
      
      stats <- plyr::ddply(stats, ply_cols(rev(facet.by)), function (z) {
        
        # Position bracket height on a per-facet basis
        z[['.step']] <- z[['.ypos']] * .13
        z[['.ypos']] <- z[['.ypos']] + (z[['.step']] * seq_len(nrow(z)))
        
        # Drop x categories that are absent from this facet
        if (params$.free_x) {
          ggdata2      <- plyr::match_df(ggdata, z, on=facet.by)
          xpos2        <- intersect(names(xpos), ggdata2[[xcol]])
          xpos2        <- setNames(seq_along(xpos2), xpos2)
          z[['.xmin']] <- unname(xpos2[as.character(z[['.cat1']])])
          z[['.xmax']] <- unname(xpos2[as.character(z[['.cat2']])])
        }
        
        return (z)
      })
    }
    stats[['.tick']] <- stats[['.ypos']] - stats[['.step']] / 4
    remove("xpos")
    
    
    stats_df <- stats[,intersect(names(stats), names(ggdata)),drop=FALSE]
    
    attr(ggdata, 'stat_labels') <- data.frame(
      check.names = FALSE, 
      stats_df,
      .x     = (stats[['.xmin']] + stats[['.xmax']]) / 2,
      .y     = stats[['.ypos']] + (stats[['.step']] / 15),
      .label = stats[['.label']] )
    
    attr(ggdata, 'stat_brackets') <- data.frame(
      check.names = FALSE, 
      dplyr::bind_rows(stats_df, stats_df, stats_df),
      .x    = c(stats[['.xmin']], stats[['.xmin']], stats[['.xmax']]),
      .xend = c(stats[['.xmax']], stats[['.xmin']], stats[['.xmax']]),
      .y    = c(stats[['.ypos']], stats[['.ypos']], stats[['.ypos']]),
      .yend = c(stats[['.ypos']], stats[['.tick']], stats[['.tick']]) )
    
    
    set_layer(
      params = params, 
      layer  = 'brackets', 
      'mapping|x'    = ".x",
      'mapping|xend' = ".xend",
      'mapping|y'    = ".y",
      'mapping|yend' = ".yend" )
    
    set_layer(
      params = params, 
      layer  = 'stats_text', 
      'size'          = 3,
      'vjust'         = ifelse(isTRUE(params$flip), 0.5, 0),
      'parse'         = length(p.label) == 1,
      'mapping|x'     = ".x",
      'mapping|y'     = ".y",
      'mapping|label' = ".label" )
    
  }
  
  
  #________________________________________________________
  # Group-wise Brackets ====
  #________________________________________________________
  else if (test == "means") {
    
    if (xcol == ".all") stats[['.all']] <- "all"
    stats[[xcol]] %<>% factor(levels = levels(ggdata[[xcol]]))
    stats[['.xpos']] <- as.numeric(stats[[xcol]])
    
    
    # When some x-values are absent from some facets, we'll need
    # to re-number the remaining x-positions.
    
    if (!is.null(facet.by) && isTRUE(params$.free_x)) {
      
      all_x <- levels(stats[[xcol]])
      
      stats %<>% plyr::ddply(ply_cols(facet.by), function (z) {
        
        facet_x <- plyr::match_df(ggdata, z, on=facet.by)[[xcol]] %>%
          as.character() %>%
          unique()
        
        map_x <- sapply(all_x, function (x) which(facet_x == x) %>% ifelse(is_null(.), NA, .))
        z[['.xpos']] <- map_x[as.numeric(z[[xcol]])]
        
        return (z)
      })
    }
    
    
    stats_df <- stats[,intersect(names(stats), names(ggdata)),drop=FALSE]
    
    attr(ggdata, 'stat_labels') <- data.frame(
      check.names = FALSE, 
      stats_df,
      .x     = stats[['.xpos']],
      .y     = stats[['.ypos']] * 1.10,
      .label = stats[['.label']] )
    
    
    if (isFALSE(params$flip)) {
      
      attr(ggdata, 'stat_brackets') <- data.frame(
        check.names = FALSE, 
        stats_df,
        .x    = stats[['.xpos']] - .4,
        .xend = stats[['.xpos']] + .4,
        .y    = stats[['.ypos']] * 1.08,
        .yend = stats[['.ypos']] * 1.08 )
      
      set_layer(
        params = params, 
        layer  = 'brackets', 
        'mapping|x'    = ".x",
        'mapping|xend' = ".xend",
        'mapping|y'    = ".y",
        'mapping|yend' = ".yend" )
    }
    
    set_layer(
      params = params, 
      layer  = 'stats_text', 
      'size'          = 3,
      'hjust'         = ifelse(isTRUE(params$flip), 0.1, 0.5),
      'vjust'         = ifelse(isTRUE(params$flip), 0.5, 0),
      'parse'         = length(p.label) == 1,
      'mapping|x'     = ".x",
      'mapping|y'     = ".y",
      'mapping|label' = ".label" )
    
  }
  
  
  
  if (isTRUE(params$caption)) {
    
    set_layer(params, 'theme', plot.caption = element_text(size = 9, face = "italic"))
    set_layer(params, 'labs',  caption      = local({
      
      test <- switch(
        EXPR = test,
        pw_means = "Mann-Whitney",
        means    = "Kruskal-Wallis" )
      
      meth <- switch(
        EXPR = params$p.adj,
        holm       = "Holm",                  # (1979)
        hochberg   = "Hochberg",              # (1988)
        hommel     = "Hommel",                # (1988)
        BH         = "Benjamini & Hochberg",  # (1995)
        fdr        = "Benjamini & Hochberg",  # (1995)
        BY         = "Benjamini & Yekutieli", # (2001)
        bonferroni = "Bonferroni",
        none       = "no" )
      
      return(glue("{test} p-values, with {meth} FDR correction."))
      
    }))
  }
  
  
  
  
  #________________________________________________________
  # Add extra space to the right for the p-values.
  #________________________________________________________
  if (params$flip) { set_layer(params, 'yaxis', expand = c(0.02, 0, 0.15, 0) )
  } else           { set_layer(params, 'yaxis', expand = c(0.02, 0, 0.08, 0) ) }
  
  
  #________________________________________________________
  # Don't label the y-axis beyond the data range.
  #________________________________________________________
  if (has_layer(params, 'facet')) {
    
    set_layer(params, 'yaxis', limits = c(0, NA))
    
  } else {
    
    set_layer(
      params = params, 
      layer  = 'yaxis',
      breaks = local({
        ymax   <- min(stats[[".ypos"]])
        breaks <- labeling::extended(0, ymax, 5, only.loose = TRUE)
        return (breaks[breaks <= ymax])
      }))
    
    set_layer(
      params = params, 
      layer  = 'stats_bg',
      'color'   = NA,
      'fill'    = "white",
      'mapping' = aes(
        xmin = !!-Inf, 
        xmax = !!Inf, 
        ymin = !!signif(min(stats[['.ypos']]), 3), 
        ymax = !!Inf ))
  }
  
  
  
  params$.ggdata <- ggdata
  
  return (invisible(params))
}



