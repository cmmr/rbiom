
plot_stacked <- function (biom, x, y, facet.by = NULL, label.by = NULL, sort.by = NULL, colors = NULL, taxa = 5, other = FALSE, xlab.angle = "auto", ...) {
  
  dots <- list(...)
  rank <- as.vector(y)
  
  
  #--------------------------------------------------------------
  # Settings for the ggplot object
  #--------------------------------------------------------------
  elements <- list('geom_col' = list(), 'labs' = list(x=NULL, fill=rank))
  
  
  #--------------------------------------------------------------
  # Limit to top n taxa, n relative abundance, or custom list
  #--------------------------------------------------------------
  if (!is.null(taxa)) {
    
    mat <- as.matrix(taxa.rollup(biom, rank = rank, sparse = T))
    
    if (is.numeric(taxa)) {
      rel <- sort(rowMeans(t(t(mat) / colSums(mat))), decreasing = TRUE)
      if (taxa >= 1) { taxa <- head(names(rel), taxa) 
      } else         { taxa <- names(rel)[rel >= taxa] }
      
    } else if (isTRUE(abbr)) {
      
      opts <- tolower(rownames(mat))
      taxa <- tolower(as.character(taxa))
      taxa <- sapply(taxa, startsWith, x = opts)
      taxa <- rownames(mat)[apply(taxa, 1L, any)]
      
    } else {
      taxa <- intersect(as.character(taxa), rownames(mat))
    }
    
    if (length(taxa) == 0)
      stop("No taxa match the criteria: ", capture.output(str(params[['taxa']])))
    
    #--------------------------------------------------------------
    # Sum the remaining taxa into an "Other" category
    #--------------------------------------------------------------
    if (length(taxa) == nrow(mat) || isFALSE(other)) {
      mat <- mat[taxa,,drop=F]
      
    } else {
      other <- matrix(
        data = colSums(mat[!rownames(mat) %in% taxa,,drop=F]), 
        nrow = 1, dimnames = list("Other", colnames(mat)))
      mat <- rbind(mat[taxa,,drop=F], other)
    }
  }
  
  
  #--------------------------------------------------------------
  # Convert counts matrix to a data.frame
  #--------------------------------------------------------------
  df <- data.frame(
    check.names = FALSE,
    '.sample' = colnames(mat)[col(mat)],
    '.taxa'   = rownames(mat)[row(mat)] %>% as.factor(),
    '.value'  = as.numeric(mat)
  )
  
  
  
  #--------------------------------------------------------------
  # Control the values displayed on the y axis
  #--------------------------------------------------------------
  if (is.rarefied(biom)) {
    df[['.value']] <- df[['.value']] / attr(biom, 'rarefaction')
    breaks         <- base::pretty(df[['.value']])
    labels         <- paste0(breaks * 100, "%")
    elements[['scale_y_continuous']] <- list(
      'expand' = c(0,0),
      'limits' = c(0,1),
      'breaks' = breaks,
      'labels' = labels )
    remove("breaks")
    
    elements[['labs']][['y']] <- "Relative Abundance"
  } else {
    elements[['labs']][['y']] <- "Raw Abundance"
  }
  
  
  
  #--------------------------------------------------------------
  # Adjust the ordering and labeling of x-axis values
  #--------------------------------------------------------------
  if (!is.null(label.by) || !is.null(sort.by)) {
    
    if (is.null(sort.by))  { breaks <- sample.names(biom)
    } else                 { breaks <- names(sort(metadata(biom, sort.by))) }
    
    if (is.null(label.by)) { labels <- breaks
    } else                 { labels <- metadata(biom, label.by)[breaks] }
    
    df[['.sample']]                <- factor(df[['.sample']], 'levels' = breaks)
    elements[['scale_x_discrete']] <- list('breaks' = breaks, 'labels' = labels)
    
    elements[['labs']][['x']] <- paste0(
      "Samples (",
      if (!is.null(label.by)) paste0("labeled by '", label.by, "'") else NULL,
      if (!is.null(sort.by) && !is.null(label.by)) " and "          else NULL,
      if (!is.null(sort.by)) paste0("sorted by '",  sort.by, "'")   else NULL,
      ")" )
  }
  
  
  
  #--------------------------------------------------------------
  # Set up the global aesthetics
  #--------------------------------------------------------------
  aes_args <- list(x = ".sample", y = ".value", fill = ".taxa")
  
  
  
  #--------------------------------------------------------------
  # Custom colors for each taxa
  #--------------------------------------------------------------
  colors <- assign_colors(colors, df[['.taxa']])
  elements[['scale_fill_manual']] <- list('values' = colors)
  
  
  
  #--------------------------------------------------------------
  # Facet arguments
  #--------------------------------------------------------------
  if (!is.null(facet.by)) {
    
    for (i in facet.by)
      df[[i]] <- metadata(biom, i)[as.character(df[['.sample']])]
    
    elements[['facet_wrap']] <- list(
      'facets' = backtick(facet.by),
      'scales' = "free_x" )
  }
  
  #--------------------------------------------------------------
  # Theme arguments
  #--------------------------------------------------------------
  elements[['theme']] <- list(
      'plot.margin' = unit(c(1,1,1,1), "lines"), 
      'panel.grid'  = element_blank(),
      'axis.text.x' = element_text(angle=-90, vjust=0.3, hjust=0)
  )
  if (isTRUE(xlab.angle == 30 || tolower(xlab.angle) == "angled"))
    elements[['theme']][['axis.text.x']] <- element_text(angle=-30, vjust=1, hjust=0)
  
  
  
  #--------------------------------------------------------------
  # Create the plot and add each layer with its arguments
  #--------------------------------------------------------------
  p <- ggplot(df, do.call(aes_string, aes_args, quote = TRUE)) +
    theme_bw()
  
  
  for (layer in names(elements)) {
    func <- do.call(`::`, list("ggplot2", layer))
    
    # Unprefixed dot arguments, e.g. 'scales'="free_x"
    #--------------------------------------------------------------
    for (i in intersect(names(dots), formalArgs(func)))
      elements[[layer]][[i]] <- dots[[i]]
    
    p <- p + do.call(func, elements[[layer]])
  }
  
  return (p)
}
