
# anno components:
#  - t Title
#  - m Method (permanova info as caption)
#  - p P-Value
#  - r R-Squared
#  - f F-Statistic

# layers components:
#  - e Ellipses
#  - p Sample Points
#  - i Sample IDs
#  - c Centroids
#  - b BiPlot Circle
#  - a BiPlot Arrow
#  - t BiPlot Text


plot_ordination <- function (biom, x, y, layers = NULL, color.by = NULL, shape.by = NULL, colors = NULL, shapes = NULL, weighted = TRUE, rank = "auto", taxa = 5, p.adj = "fdr", anno = "tmprf", ...) {
  
  dots <- list(...)
  
  
  #-----------------------------------------------
  # Sanity Check
  #-----------------------------------------------
  if (nsamples(biom) < 4)
    stop("At least four samples are needed for an ordination.")
  
  if (is.character(taxa))
    taxa <- validate_metrics(biom, taxa, mode='taxon', multi=TRUE)
  
  
  
  #-----------------------------------------------
  # Computed defaults
  #-----------------------------------------------
  if (is.null(layers) || layers %in% c("", "box"))
    layers <- ifelse(is.null(color.by), "p", "pce")
  
  
  
  #-----------------------------------------------
  # User-editable layers
  #-----------------------------------------------
  layers <- strsplit(layers, '')[[1]] %>% 
    tolower() %>% 
    sapply(simplify = F, function (i) {
      list(
        'fun' = list(
          'p' = ggplot2::geom_point,   'b' = ggplot2::geom_point,
          'i' = ggplot2::geom_text,    't' = ggrepel::geom_label_repel,
          'c' = ggplot2::geom_segment, 'a' = ggplot2::geom_segment,
          'e' = ggplot2::stat_ellipse )[[i]],
        'args' = list()
      )
    })
  
  
  #-----------------------------------------------
  # Non-editable layers
  #-----------------------------------------------
  anno     <- strsplit(anno, '')[[1]] %>% tolower()
  elements <- list(
    'labs'               = list(x = NULL, y = NULL),
    'scale_x_continuous' = list(),
    'scale_y_continuous' = list() )
  
  
  
  #-----------------------------------------------
  # Compute distance matrix and ordination.
  #-----------------------------------------------
  dm  <- beta.div(biom, method = y, weighted = weighted)
  ord <- ordinate(dm, ord = x, k = 2, safe = TRUE)
  
  colnames(ord)[1:2] <- c(".x", ".y")
  ord[['.label']]    <- rownames(ord)
  for (i in unique(c(color.by, shape.by)))
    ord[[i]] <- metadata(biom, i)[rownames(ord)] %>% as.factor()
  
  
  #-----------------------------------------------
  # Common aes arguments to subset as needed
  #-----------------------------------------------
  sample_aes <- function (keys) {
    args <- list(x = ".x", y = ".y", label = ".label")
    if (!is.null(color.by)) args[['color']] <- backtick(color.by)
    if (!is.null(shape.by)) args[['shape']] <- backtick(shape.by)
    do.call(aes_string, args[intersect(keys, names(args))])
  }
  
  
  
  #-----------------------------------------------
  # Plot Title
  #-----------------------------------------------
  if ("t" %in% anno)
    elements[['labs']][['title']] <- paste(
      ifelse(weighted, "Weighted", "Unweighted"), y, x)
  
  
  #-----------------------------------------------
  # Sample points
  #-----------------------------------------------
  if ("p" %in% names(layers))
    layers[['p']][['args']][['mapping']] <- sample_aes(
      keys = c("x", "y", "color", "shape") )
  
  
  #-----------------------------------------------
  # Sample IDs
  #-----------------------------------------------
  if ("i" %in% names(layers))
    layers[['i']][['args']][['mapping']] <- sample_aes(
      keys = c("x", "y", "color", "label") )
  
  
  #-----------------------------------------------
  # Ellipses
  #-----------------------------------------------
  if ("e" %in% names(layers))
    layers[['e']][['args']][['mapping']] <- sample_aes(
      keys = c("x", "y", "color") )
  
  
  #-----------------------------------------------
  # Centroids
  #-----------------------------------------------
  if ("c" %in% names(layers)) {
    
    centroids <- ord[,c('.x', '.y', color.by),drop=F]
    
    if (is.null(color.by)) {
      centroids[['.xend']] <- mean(centroids[['.x']])
      centroids[['.yend']] <- mean(centroids[['.y']])
      
    } else {
      aggData <- aggregate(as.formula(paste("cbind(.x, .y) ~", backtick(color.by))), data=centroids, mean)
      names(aggData) <- c(color.by, ".xend", ".yend")
      centroids <- merge(centroids, aggData, by=color.by, all.x=TRUE)
      remove("aggData")
    }
    
    aes_args <- list(x = ".x", y = ".y", xend='.xend', yend='.yend')
    if (!is.null(color.by)) aes_args[['color']] <- backtick(color.by)
    
    layers[['c']][['args']][['data']]    <- centroids
    layers[['c']][['args']][['mapping']] <- do.call(aes_string, aes_args)
    layers[['c']][['args']][['size']]    <- 0.75
    layers[['c']][['args']][['alpha']]   <- 0.4
    
    remove("centroids", "aes_args")
  }
  
  
  #-----------------------------------------------
  # BiPlot
  #-----------------------------------------------
  if (any(c("b", "a", "t") %in% names(layers))) {
    
    # Determine which taxonomic rank to overlay
    if (is.null(rank) || rank %in% c("", "auto")) {
      if (is.character(taxa)) {
        rank <- names(which.max(apply(taxonomy(biom), 2L, function (x) sum(x %in% taxa))))
      } else {
        rank <- tail(c('OTU', taxa.ranks(biom)), 1)
      }
    } else {
      rank <- validate_metrics(biom, rank, mode="rank")
    }
    
    # Calculate Weighted average for each taxon
    biplot <- distill(biom, rank, safe = TRUE, md = FALSE)
    biplot[['.value']] <- biplot[['.value']] / sample.sums(biom)[biplot[['.sample']]]
    biplot[['.x']] <- ord[biplot[['.sample']], '.x']
    biplot[['.y']] <- ord[biplot[['.sample']], '.y']
    center.x <- mean(biplot[['.x']])
    center.y <- mean(biplot[['.y']])
    
    # Limit to specific list of taxa
    if (is.character(taxa))
      biplot <- biplot[biplot[['.taxa']] %in% taxa,,drop=F]
    
    # Get x,y,size and p-value for each taxon
    biplot <- plyr::ddply(biplot, ".taxa", function(x) {
      
      # Calculate Weighted average for each taxon
      Abundance <- sum(x[['.value']])
      if (Abundance == 0) return(NULL)
      Axis.1 <- sum(x[['.x']] * x[['.value']]) / Abundance
      Axis.2 <- sum(x[['.y']] * x[['.value']]) / Abundance
      
      result <- data.frame(
        '.x'     = Axis.1, 
        '.y'     = Axis.2,
        '.value' = Abundance/nrow(x), 
        '.pVal'  = NA )
      
      # Only compute p-value if we're filtering by p-value.
      # Use random reassignment of abundances to calculate 
      # p-value. Known as approximate/random/Monte Carlo
      # permutation test.
      
      if (isTRUE(taxa < 1)) {
        
        ptest <- boot::boot(data=x, R=1000, statistic=function(data, ind) {
          r.x <- sum(data[['.x']] * data[ind,'.value']) / Abundance
          r.y <- sum(data[['.y']] * data[ind,'.value']) / Abundance
          d   <- sqrt((r.x - center.x)^2 + (r.y - center.y)^2)
          return (d)
        })
  
        d <- sqrt((Axis.1 - center.x)^2 + (Axis.2 - center.y)^2)
        result[['.pVal']] <- length(which(ptest[['t']] >= d)) / ptest[['R']]
      }
  
      return (result)
    })
    
    # Limit to taxa with adj. p-value <= n
    if (isTRUE(taxa < 1)) {
      biplot[['.pVal']] <- p.adjust(biplot[['.pVal']], method=p.adj)
      biplot            <- biplot[biplot[['.pVal']] <= taxa,,drop=F]
    }
    
    # Limit to top n most abundant taxa
    if (isTRUE(taxa >= 1)) {
      biplot <- biplot[tail(order(biplot[['.value']]), taxa),,drop=F]
    }
  
  
    # BiPlot points and text share the same size scale
    biplot[['point.size']] <- biplot[['.value']]
    biplot[['label.size']] <- biplot[['.value']] / 20
    if (!"b" %in% names(layers)) biplot[['point.size']] <- 0
    
    attr(biplot, 'center.x') <- center.x
    attr(biplot, 'center.y') <- center.y
    
    remove("center.x", "center.y")
  }
  
  # browser()
  
  
  #-----------------------------------------------
  # BiPlot Arrows
  #-----------------------------------------------
  if ("a" %in% names(layers))
    layers[['a']][['args']] <- list(
      'data'    = biplot[,c(".x", ".y"),drop=F],
      'mapping' = aes(x=.x, y=.y),
      'xend'    = attr(biplot, 'center.x'),
      'yend'    = attr(biplot, 'center.y'),
      'color'   = "darkgray", 
      'size'    = 0.75, 
      'alpha'   = 0.4,
      'arrow'   = arrow(ends="first", length=unit(.5,"cm"))
    )

  
  #-----------------------------------------------
  # BiPlot Points
  #-----------------------------------------------
  if ("b" %in% names(layers))
    layers[['b']][['args']] <- list(
      'data'    = biplot[,c(".x", ".y", "point.size"),drop=F],
      'mapping' = aes(x=.x, y=.y, size=point.size),
      'color'   = "darkgray", 
      'alpha'   = 0.5
    )
  

  #-----------------------------------------------
  # BiPlot Taxa Names
  #-----------------------------------------------
  if ("t" %in% names(layers))
    layers[['t']][['args']] <- list(
      'data'    = biplot[,c(".x", ".y", ".taxa", "label.size", "point.size"),drop=F],
      'mapping' = aes(
        'x'          = .x, 
        'y'          = .y, 
        'label'      = .taxa, 
        'size'       = label.size,
        'point.size' = point.size * 2 ),
      'show.legend'        = FALSE,
      'fill'               = alpha(c("white"), 0.8),
      'box.padding'        = 1,
      'segment.curvature'  = -0.1, 
      'segment.linetype'   = 8, 
      'seed'               = 0
    )
  
  
  #-----------------------------------------------
  # Colors, shapes, and sizes.
  #-----------------------------------------------
  
  if (!is.null(color.by)) {
    colors <- assign_colors(colors, ord[[color.by]])
    elements[['scale_color_manual']] <- list('values' = colors)
    elements[['scale_fill_manual']]  <- list('values' = colors)
  }
  
  if (!is.null(shape.by)) {
    shapes <- assign_shapes(shapes, ord[[shape.by]])
    elements[['scale_shape_manual']] <- list('values' = shapes)
  }
  
  if (any(c("b", "t") %in% names(layers)))
    elements[['scale_size_continuous']] <- list(
      'range'  = c(3, ifelse("b" %in% names(layers), 15, 5)),
      'name'   = "Taxa Abundance",
      'labels' = function (x) paste0(x * 100, "%") )
  
  
  
  #-----------------------------------------------
  # Adonis / PERMANOVA
  #-----------------------------------------------
  if (!is.null(color.by) && any(c("p", "r", "f") %in% anno)) {

    mtx   <- as.matrix(dm)[rownames(ord), rownames(ord)]
    perms <- t(replicate(999, sample(1:nrow(mtx), nrow(mtx))))
    obj   <- try(vegan::adonis(mtx ~ ord[[color.by]], permutations=perms), silent=TRUE)
    
    if (is(obj, "adonis")) {
      
      elements[['labs']][['subtitle']] <- paste(collapse="; ", c(
        if ("p" %in% anno) paste("P-Value:",     signif(obj$aov.tab[['Pr(>F)']][[1]],  3))  else NULL,
        if ("r" %in% anno) paste("R-Squared:",   signif(obj$aov.tab[['R2']][[1]],      3))  else NULL,
        if ("f" %in% anno) paste("F-Statistic:", signif(obj$aov.tab[['F.Model']][[1]], 3)) else NULL
      ))
      
      if ("m" %in% anno)
        elements[['labs']][['caption']] <- "Statistics computed with Adonis (1000 permutations)."
    }
    remove("mtx", "perms", "obj")
  }
  
  
  
  #-----------------------------------------------
  # Load sample ordination into ggplot
  #-----------------------------------------------
  p <- ggplot2::ggplot(ord) +
    theme_bw() +
    theme(
      panel.border     = element_rect(colour="black", fill=F, size=1),
      panel.grid.major = element_blank(), # remove major grid
      panel.grid.minor = element_blank(), # remove minor grid
      panel.background = element_rect(fill="white") )
  
  
  #-----------------------------------------------
  # Add customizable layers - points, lines, ...
  #-----------------------------------------------
  for (layer in names(layers)) {
    
    qry <- paste0(layer, ".")
    for (i in names(dots)[startsWith(names(dots), qry)])
      layers[[layer]][['args']][[sub(qry, "", i)]] <- dots[[i]]
    
    fun <- layers[[layer]][['fun']]
    p   <- p + do.call(fun, layers[[layer]][['args']])
  }
  
  
  #-----------------------------------------------
  # Add other plot settings - scales, labels, ...
  #-----------------------------------------------
  for (layer in names(elements)) {
    fun <- do.call(`::`, list("ggplot2", layer))
    p   <- p + do.call(fun, elements[[layer]])
  }
  
  
  return(p)
  
}

