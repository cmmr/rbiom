
#________________________________________________________
# Computes p-values for categorical differences
#________________________________________________________
boxplot_stats <- function (params) {
  
  stopifnot(is_bare_environment(params))
  
  
  layers   <- params$layers
  ggdata   <- params$.ggdata
  xcol     <- params$.xcol
  ycol     <- params$.ycol
  stat.by  <- params$stat.by
  facet.by <- params$facet.by
  p.label  <- params$p.label
  test     <- params$test
  
  
  
  #________________________________________________________
  # No work to do.
  #________________________________________________________
  if (plyr::empty(ggdata) || is.null(stat.by) || all(p.label == 0) || eq(test, 'none')) {
    set_layer(params, 'yaxis', expand = c(0.02, 0, 0.02, 0))
    return (invisible(params))
  }
  
  
  #________________________________________________________
  # Groupwise when xcol != stat.by.
  #________________________________________________________
  else if (!eq(xcol, stat.by)) {
    test  <- "kruskal"
    stats <- stats_table(
      df       = ggdata, 
      stat.by  = stat.by, 
      resp     = ycol,
      test     = test, 
      split.by = c(xcol, facet.by), 
      p.adj    = params$p.adj )
  }
  
  
  #________________________________________________________
  # Pairwise .
  #________________________________________________________
  else {
    test  <- "wilcox"
    stats <- stats_table(
      df       = ggdata, 
      stat.by  = xcol, 
      resp     = ycol,
      test     = test, 
      split.by = facet.by, 
      p.adj    = params$p.adj )
  }
  
  
  params$.plot_attrs[['stats']] <- stats
  
  
  #________________________________________________________
  # Implement `p.top` for taxa_boxplot().
  #________________________________________________________
  params <- apply_p.top(params)
  ggdata <- params$.ggdata
  stats  <- params$.plot_attrs[['stats']]
  if (plyr::empty(ggdata)) return (invisible(params))
  
  
  
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
        max    = {function (v) { max(v)                        }},
        mean   = {function (v) { mean(v)                       }},
        violin = {function (v) { max(stats::density(v)[['x']]) }},
        box    = switch(as.character(isFALSE(params$outliers)), 
          'TRUE'  = {function (v) { grDevices::boxplot.stats(v)$stats[5]             }},
          'FALSE' = {function (v) { max(v)                                           }} ),
        vline  = switch(params$ci, 
          'range' = {function (v) { max(v)                                           }},
          'mad'   = {function (v) { median(v) + stats::mad(v, median(v))             }},
          'sd'    = {function (v) { mean(v) + sd(v)                                  }},
          'se'    = {function (v) { mean(v) + sqrt(var(v)/length(v))                 }},
          'ci'    = {function (v) { t.test(v, conf.level = params$level)$conf.int[2] }} ))
    
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
  if (!isTRUE(params$.free_y)) {
    
    stats[['.ypos']] <- ggdata %>%
      plyr::ddply(
        .variables = ply_cols(unique(c(xcol, stat.by, facet.by))), 
        .fun       = function (df) data.frame(.ypos = y_pos_fun(df)) ) %>%
      dplyr::pull('.ypos') %>%
      max(0, na.rm = TRUE)
    
  } else {
    
    stats %<>% dplyr::left_join(
      by = facet.by, 
      y  = ggdata %>% 
        plyr::ddply(ply_cols(facet.by), function (df) {
          plyr::daply(df, ply_cols(unique(c(xcol, stat.by))), y_pos_fun) %>%
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
  if (test == "wilcox") {
    
    
    # Convert e.g. 'Saliva - Stool' to 'Saliva' and 'Stool'.
    #________________________________________________________
    pairs <- combn(levels(ggdata[[xcol]]), m = 2, simplify = FALSE)
    names(pairs) <- sapply(pairs, paste, collapse = ' - ')
    
    stats[['.cat1']]  <-  unname(sapply(as.character(stats[[xcol]]), function (p) pairs[[p]][[1]]))
    stats[['.cat2']]  <-  unname(sapply(as.character(stats[[xcol]]), function (p) pairs[[p]][[2]]))
    stats[['.cat1']] %<>% factor(levels(ggdata[[xcol]]))
    stats[['.cat2']] %<>% factor(levels(ggdata[[xcol]]))
    remove("pairs")
    
    
    # Line segments connecting category pairs.
    #________________________________________________________
    stats <- stats %>% 
      plyr::ddply(ply_cols(rev(facet.by)), function (z) {
      
        # Drop x categories that are absent from this facet
        if (isTRUE(params$.free_x)) {
          
          # xcats <- dplyr::left_join(
          #     x  = z[,facet.by,FALSE],
          #     y  = ggdata[,c(xcol, facet.by),FALSE],
          #     by = facet.by ) %>%
          #   dplyr::pull(xcol) %>%
          #   unique() %>%
          #   intersect(x = levels(ggdata[[xcol]]))
          
          xcats <- c(as.character(z[['.cat1']]), as.character(z[['.cat2']]))
          xcats <- intersect(levels(ggdata[[xcol]]), xcats)
          
          z[['.cat1']] %<>% factor(levels = xcats)
          z[['.cat2']] %<>% factor(levels = xcats)
        }
        
        # Convert e.g. c('Saliva', 'Stool') to c(4, 5)
        z[['.xmin']] <- as.numeric(z[['.cat1']])
        z[['.xmax']] <- as.numeric(z[['.cat2']])
        
        # Keep only the fields we need.
        z <- z[,c('.xmin', '.xmax', '.ypos', '.label'), drop = FALSE]
        
        # Position bracket height on a per-facet basis
        z[['.step']] <- z[['.ypos']] * .13
        z[['.ypos']] <- z[['.ypos']] + (z[['.step']] * (seq_len(nrow(z)) - 1))
        z[['.tick']] <- z[['.ypos']] - z[['.step']] / 4
        
        return (z)
      }) %>%
      dplyr::select(-any_of('.id'))
    
    
    .xmin <- .xmax <- .step <- .ypos <- .tick <- .label <- NULL # for CRAN check only
    
    attr(ggdata, 'stat_labels') <- tibble(
        stats,
        .x = (stats[['.xmin']] + stats[['.xmax']]) / 2,
        .y = stats[['.ypos']] + (stats[['.step']] / 15) ) %>%
      dplyr::select(-c(.xmin, .xmax, .step, .ypos, .tick))
    
    attr(ggdata, 'stat_brackets') <- tibble(
        dplyr::bind_rows(stats, stats, stats),
        .x    = c(stats[['.xmin']], stats[['.xmin']], stats[['.xmax']]),
        .xend = c(stats[['.xmax']], stats[['.xmin']], stats[['.xmax']]),
        .y    = c(stats[['.ypos']], stats[['.ypos']], stats[['.ypos']]),
        .yend = c(stats[['.ypos']], stats[['.tick']], stats[['.tick']]) ) %>%
      dplyr::select(-c(.xmin, .xmax, .step, .ypos, .tick, .label))
    
    remove('.xmin', '.xmax', '.step', '.ypos', '.tick', '.label')
    
    
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
  else if (test == "kruskal") {
    
    # Vertically align all group-wise brackets.
    if (isFALSE(params$.free_y)) stats$.ypos %<>% max()
    
    
    attr(ggdata, 'stat_labels') <- stats %>% 
      plyr::ddply(ply_cols(c(xcol, facet.by)), function (z) {
        tibble(!!ycol := z$.ypos[[1]] * 1.10, .label = z$.label[[1]]) })
    
    set_layer(
      params = params, 
      layer  = 'stats_text', 
      'size'          = 3,
      'hjust'         = ifelse(isTRUE(params$flip), 0.1, 0.5),
      'vjust'         = ifelse(isTRUE(params$flip), 0.5, 0),
      'parse'         = length(p.label) == 1,
      'mapping|label' = ".label" )
    
    
    
    if (isFALSE(params$flip)) {
      
      attr(ggdata, 'stat_brackets') <- stats %>% 
        plyr::ddply(ply_cols(facet.by), function (z) {
          
          if (is.null(xcol)) {
            z <- head(z, 1)
            z[['.x']]    <- -0.4
            z[['.xend']] <-  0.4
          } else {
            
            # Drop x categories that are absent from this facet
            if (isTRUE(params$.free_x)) {
              
              # xcats <- dplyr::left_join(
              #   x  = z[,facet.by,FALSE],
              #   y  = ggdata[,c(xcol, facet.by),FALSE],
              #   by = facet.by ) %>%
              #   dplyr::pull(xcol) %>%
              #   unique() %>%
              #   intersect(x = levels(ggdata[[xcol]]))
              
              xcats <- as.character(z[[xcol]])
              xcats <- intersect(levels(ggdata[[xcol]]), xcats)
              
              z[[xcol]] %<>% factor(levels = xcats)
            }
            
            z <- z[which(!duplicated(z[[xcol]])),,drop=FALSE]
            z[['.x']]    <- as.numeric(z[[xcol]]) - .4
            z[['.xend']] <- as.numeric(z[[xcol]]) + .4
          }
          
          z[['.y']] <- z[['.ypos']] * 1.08
          
          .x <- .xend <- .y <- NULL # for CRAN check only
          
          return (dplyr::select(z, c(.x, .xend, .y)))
        }) %>%
        dplyr::select(-any_of('.id'))
      
      set_layer(
        params = params, 
        layer  = 'brackets', 
        'mapping|x'    = ".x",
        'mapping|xend' = ".xend",
        'mapping|y'    = ".y",
        'mapping|yend' = ".y" )
        
    }
    
    
  }
  
  
  
  if (isTRUE(params$caption)) {
    
    set_layer(params, 'theme', plot.caption = element_text(size = 9, face = "italic"))
    set_layer(params, 'labs',  caption      = local({
      
      test <- switch(
        EXPR = test,
        wilcox  = "Mann-Whitney",
        kruskal = "Kruskal-Wallis" )
      
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
  if (isTRUE(params$.free_y)) {
    
    set_layer(params, 'yaxis', limits = c(0, NA))
    
  } else {
    
    ymax   <- min(stats[[".ypos"]])
    breaks <- labeling::extended(0, ymax, 10, only.loose = TRUE)
    breaks <- breaks[breaks < ymax]
    
    set_layer(
      params = params, 
      layer  = 'yaxis',
      breaks       = breaks[seq_along(breaks) %% 2 == 1],
      minor_breaks = breaks[seq_along(breaks) %% 2 == 0] )
    
    set_layer(
      params = params, 
      layer  = 'stats_bg',
      'color'   = NA,
      'fill'    = "white",
      'mapping' = aes(
        xmin = !!-Inf, 
        xmax = !!Inf, 
        ymin = !!signif(ymax, 3), 
        ymax = !!Inf ))
  }
  
  
  
  params$.ggdata <- ggdata
  
  return (invisible(params))
}



