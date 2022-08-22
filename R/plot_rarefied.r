plot_rarefied <- function(biom, x, color.by = NULL, shape.by = NULL, facet.by = NULL, colors = NULL, shapes = NULL, rline = NULL, ...) {
  
  dots <- list(...)
  
  
  #===============================================
  # Validate metadata fields
  #===============================================
  if (!is.null(color.by)) color.by <- validate_metrics(biom, color.by, mode="meta")
  if (!is.null(shape.by)) shape.by <- validate_metrics(biom, shape.by, mode="meta")
  if (!is.null(facet.by)) facet.by <- validate_metrics(biom, facet.by, mode="meta")
  
  
  #===============================================
  # Collect ggplot aesthetics and elements here
  #===============================================
  aes_args <- list()
  elements <- list('theme_bw' = list())
  
  
  
  if (isTRUE(tolower(x) == "reads")) {
  
    #===============================================
    # Reads Retained: Sample (x) vs Depth (y)
    #===============================================
    
    ss <- sort(sample_sums(biom))
    df <- data.frame(
      .sample = factor(names(ss), levels=names(ss)),
      .all    = unname(ss)
    )
    
    aes_args                    <- list(x=".sample", y=".all")
    elements[['labs']]          <- list(x="Sample", y="Reads")
    elements[['geom_col']]      <- list()
    elements[['scale_y_log10']] <- loglabels(df[['.all']])
    elements[['theme']]         <- list(
      legend.position    = c(0.2, 0.88),
      legend.direction   = "vertical",
      legend.background  = element_rect(fill=NA),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      panel.grid.major.x = element_blank()  )
    
    if (!is.null(rline)) {
      
      elements[['geom_hline']] <- list('yintercept' = rline, linetype="dotted")
      
      df[['Retained']] <- ifelse(df[['.all']] >= rline, rline, 0)
      df[['Excluded']] <- ifelse(
        test = df[['.all']] > rline,
        yes  = (df[['.all']] - df[['Retained']]) / rline + 1,
        no   = (df[['.all']] - df[['Retained']])
      )
      
      nTotal <- sum(df[['.all']])
      nKept  <- sum(df[['Retained']])
      elements[['ggtitle']] <- list(sprintf(
        fmt = "Reads Retained: %s / %s = %s%%", 
        scales::comma(nKept), 
        scales::comma(nTotal), 
        as.character(signif((nKept / nTotal) * 100, digits=3)) ))
      remove("nTotal", "nKept")
      
      df <- tidyr::pivot_longer(
        data      = df, 
        cols      = c('Retained', 'Excluded'), 
        names_to  = ".color", 
        values_to = ".value") %>%
        dplyr::filter(.value > 0) %>%
        dplyr::arrange(.color)
      
      aes_args[['y']] <- ".value"
      
      if (is.null(colors)) colors <- c('Excluded' = "#ED5F16", 'Retained' = "#00B9EB")
      aes_args[['fill']]  <- ".color"
      aes_args[['color']] <- ".color"
      elements[['scale_fill_manual']]  <- list(values=colors)
      elements[['scale_color_manual']] <- list(values=colors)
      elements[['labs']][['fill']]     <- ""
      elements[['labs']][['color']]    <- ""
      
    }
    
    
  } else if (isTRUE(tolower(x) == "samples")) {
    
    #===============================================
    # Samples Retained: Depth (x) vs Rank (y)
    #===============================================
    
    ss <- sort(sample_sums(biom))
    df <- data.frame(
      .reads  = unname(ss), 
      .rank   = seq_along(ss),
      .sample = names(ss)
    )
    
    aes_args                    <- list(x=".reads", y=".rank")
    elements[['labs']]          <- list(x="Reads", y="Sample Rank")
    elements[['geom_point']]    <- list()
    elements[['scale_x_log10']] <- loglabels(df[['.reads']])
    
    
    #-----------------------------------------------
    # How to color the samples
    #-----------------------------------------------
    if (!is.null(color.by)) {
      
      df[[color.by]]      <- metadata(biom, color.by)[df[['.sample']]]
      aes_args[['color']] <- backtick(color.by)
      
      if (!is.null(colors))
        elements[['scale_color_manual']][['values']] <- colors
      
    } else if (!is.null(rline)) {
      
      df[['.rline']] <- factor(ifelse(ss < max(c(0, rline)), "Excluded", "Retained"))
      if (is.null(colors)) colors <- c("Excluded"="#ED5F16", "Retained"="#00B9EB")
      
      aes_args[['color']]           <- ".rline"
      elements[['labs']][['color']] <- ""
      
      elements <- c(elements, list(
        'scale_color_manual' = list(values=colors),
        'guides'             = list(color=guide_legend(override.aes=list(size=4))),
        'theme'              = list(
          legend.position   = c(0.2, 0.88),
          legend.direction  = "vertical",
          legend.background = element_rect(fill=NA),
          legend.key        = element_rect(fill=NA) )
      ))
      
    } 
    
    
    #-----------------------------------------------
    # How to shape the samples
    #-----------------------------------------------
    if (!is.null(shape.by)) {
      
      df[[shape.by]]      <- metadata(biom, shape.by)[df[['.sample']]]
      aes_args[['shape']] <- backtick(shape.by)
      
      if (!is.null(shapes))
        elements[['scale_shape_manual']][['values']] <- shapes
      
    }
    
    
    #-----------------------------------------------
    # Draw a vertical line at the rarefaction depth
    #-----------------------------------------------
    if (!is.null(rline)) {
      nTotal <- scales::comma(length(ss))
      nKept  <- scales::comma(sum(ss >= rline))
      nPct   <- as.character(signif((sum(ss >= rline) / length(ss)) * 100, digits=3))
      
      elements <- c(elements, list(
       'geom_vline' = list(xintercept=rline, linetype="dotted"),
       'ggtitle'    = list(sprintf("Samples Retained: %s / %s = %s%%", nKept, nTotal, nPct))
      ))
      
      remove("nTotal", "nKept", "nPct")
    }
    
  } else {
    
    adiv_metric <- try(validate_metrics(biom, x, mode="adiv"), silent = TRUE)
    
    if (is.character(adiv_metric)) {
  
      #===============================================
      # Rarefaction Curve
      #===============================================
      
      df <- adiv_table(biom, rarefy = "multi", metrics = adiv_metric, long = TRUE, safe = TRUE)
      
      aes_args                         <- list('x' = ".depth", 'y' = ".value")
      elements[['labs']]               <- list('x' = "Rarefaction Depth", 'y' = adiv_metric)
      elements[['stat_smooth']]        <- list('method' = "lm", 'formula' = "y ~ log(x)")
      elements[['scale_x_continuous']] <- list('labels' = si_units)
      
      if (!is.null(color.by)) {
        df[[color.by]] <- metadata(biom, color.by)[df[['.sample']]]
        aes_args[['fill']]  <- backtick(color.by)
        aes_args[['color']] <- backtick(color.by)
      }
      
      if (!is.null(rline))
        elements[['geom_vline']] = list('xintercept' = rline, 'linetype' = "dotted")
      
    } else {
      stop("Don't know how to plot Rarefied ~ ", backtick(x), ".")
    }
  }
  
  
  
  #===============================================
  # Faceting
  #===============================================
  if (!is.null(facet.by)) {
    
    for (i in facet.by)
      df[[i]] <- metadata(biom, i)[df[['.sample']]]
    
    if (length(facet.by) == 2) {
      elements[['facet_grid']] <- list(
        'rows' = as.formula(paste(
          backtick(facet.by[[2]]), "~", backtick(facet.by[[1]])
        )))
      for (i in intersect(formalArgs(facet_wrap), names(dots)))
        elements[['facet_grid']][[i]] <- dots[[i]]
    }
    
    if (length(facet.by) == 1) {
      elements[['facet_wrap']] <- list(
        'facets' = backtick(facet.by[[1]]) )
      for (i in intersect(formalArgs(facet_wrap), names(dots)))
        elements[['facet_wrap']][[i]] <- dots[[i]]
    }
    
    
    #-----------------------------------------------
    # Undo automatic legend positioning
    #-----------------------------------------------
    del <- grepl("^legend\\.", names(elements[['theme']]))
    elements[['theme']] <- elements[['theme']][!del]
    remove("del")
  }
  
  
  #===============================================
  # Add in additional dots arguments
  #===============================================
  for (fn in names(elements))
    for (i in intersect(names(dots), formalArgs(fn)))
      if (!i %in% c('y', 'formula'))
        elements[[fn]][[i]] <- dots[[i]]
  
  
  #===============================================
  # Assemble the ggplot object
  #===============================================
  gg <- ggplot(df, do.call(aes_string, aes_args, quote = TRUE))
  for (fn in names(elements))
    gg <- gg + do.call(fn, elements[[fn]])
  
  
  return (gg)
 
}
