
# anno components:
#  - t Title
#  - m Method (permanova info as caption)
#  - p P-Value
#  - r R-Squared
#  - f F-Statistic
#  - a Axes

# layers components:
#  - e Ellipses
#  - p Sample Points
#  - i Sample IDs
#  - c Centroids
#  - b BiPlot Circle
#  - a BiPlot Arrow
#  - t BiPlot Text


plot_ordination <- function (
  biom, x, y, layers = NULL, 
  color.by = NULL, shape.by = NULL, facet.by = NULL, 
  colors = NULL, shapes = NULL, 
  weighted = TRUE, rank = "auto", taxa = 5, p.adj = "fdr", anno = "tmprf", ...) {
  
  dots <- list(...)
  
  
  #-----------------------------------------------
  # Sanity Check
  #-----------------------------------------------
  if (nsamples(biom) < 4)
    stop("At least four samples are needed for an ordination.")
  
  if (is.character(taxa))
    taxa <- validate_metrics(biom, taxa, mode='taxon', multi=TRUE)
  
  
  
  #-----------------------------------------------
  # User-editable layers
  #-----------------------------------------------
  layers <- tolower(layers)
  layerlist <- c(
    "point" = "p", "centroid" = "c", "ellipse" = "e", "name" = "n", 
    "mean"  = "m", "taxon" = "t",    "arrow" = "a" )
  if (all(is.na(pmatch(layers, names(layerlist)))))
    layers <- strsplit(layers, '')[[1]]
  layers <- unname(layerlist[pmatch(tolower(layers), names(layerlist))])
  layers <- layers[!is.na(layers)]
  
  if (is.null(layers))
    layers <- if(is.null(color.by)) "p" else c("p", "c", "e")
  
  layers <- sapply(layers, simplify = F, function (i) {
      list(
        'fun' = list(
          'p' = ggplot2::geom_point,   'm' = ggplot2::geom_point,
          'n' = ggplot2::geom_text,    't' = ggrepel::geom_label_repel,
          'c' = ggplot2::geom_segment, 'a' = ggplot2::geom_segment,
          'e' = ggplot2::stat_ellipse )[[i]],
        'args' = list()
      )
    })
  
  
  #-----------------------------------------------
  # Annotations
  #-----------------------------------------------
  anno <- tolower(anno)
  annolist <- c(
    "title" = "t", "method" = "m", "p-value" = "p", 
    "r-squared"  = "r", "f-statistic" = "f", "axes" = "a" )
  if (all(is.na(pmatch(anno, names(annolist)))))
    anno <- strsplit(anno, '')[[1]]
  anno <- unname(annolist[pmatch(tolower(anno), names(annolist))])
  anno <- anno[!is.na(anno)]
  
  if (is.null(anno))
    anno <- c("t", "m", "p", "r", "f")
  
  
  elements <- list(
    'scale_x_continuous' = list(),
    'scale_y_continuous' = list(),
    'labs'               = list('x' = NULL, 'y' = NULL),
    'theme'              = list(
      'panel.border'     = element_rect(color = "black", fill = FALSE, size = 1),
      'panel.grid.major' = element_blank(), # remove major grid
      'panel.grid.minor' = element_blank(), # remove minor grid
      'panel.background' = element_rect(fill = "white") )
  )
  if (isFALSE("a" %in% anno)) {
    elements[['theme']][['axis.text']]  <- element_blank()
    elements[['theme']][['axis.ticks']] <- element_blank()
  }
  
  
  #-----------------------------------------------
  # Determine which taxonomic rank to overlay
  #-----------------------------------------------
  if (any(c("m", "a", "t") %in% names(layers))) {
    
    if (is.null(rank) || rank %in% c("", "auto")) {
      if (is.character(taxa)) {
        rank <- names(which.max(apply(taxonomy(biom), 2L, function (x) sum(x %in% taxa))))
      } else {
        rank <- tail(c('OTU', taxa.ranks(biom)), 1)
      }
    } else {
      rank <- validate_metrics(biom, rank, mode="rank")
    }
  }
  
  
  
  #-----------------------------------------------
  # All the plot layers share a single data.frame
  #-----------------------------------------------
  ggdata <- bdply(biom, facet.by, function (b) {
    
    #-----------------------------------------------
    # Compute distance matrix and ordination.
    #-----------------------------------------------
    dm  <- beta.div(b, method = y, weighted = weighted)
    ord <- ordinate(dm, ord = x, k = 2, safe = TRUE)
    dat <- data.frame(
      '.src'    = "ord",
      '.sample' = rownames(ord),
      '.label'  = rownames(ord),
      '.x'      = ord[['.axis.1']], 
      '.y'      = ord[['.axis.2']],
      '.xend' = NA, '.yend' = NA, '.size' = NA
    )
      
    
    #-----------------------------------------------
    # Adonis / PERMANOVA
    #-----------------------------------------------
    if (!is.null(color.by) && any(c("p", "r", "f") %in% anno)) {
      
      ids   <- sample.names(b)
      mtx   <- as.matrix(dm)[ids, ids]
      perms <- t(replicate(999, sample(1:nrow(mtx), nrow(mtx))))
      color <- unname(metadata(b, color.by))
      obj   <- try(vegan::adonis(mtx ~ color, permutations=perms), silent=TRUE)
      
      if (is(obj, "adonis"))
        dat <- rbind(dat,
          data.frame(
            '.src'    = "stats",
            '.sample' = NA,
            '.label'  = paste(collapse="; ", c(
                if ("p" %in% anno) paste("P-Value:",     signif(obj$aov.tab[['Pr(>F)']][[1]],  3)) else NULL,
                if ("r" %in% anno) paste("R-Squared:",   signif(obj$aov.tab[['R2']][[1]],      3)) else NULL,
                if ("f" %in% anno) paste("F-Statistic:", signif(obj$aov.tab[['F.Model']][[1]], 3)) else NULL
              )),
            '.x' = NA, '.y' = NA, '.xend' = NA, '.yend' = NA, '.size' = NA
          ))
      
      remove("ids", "mtx", "perms", "color", "obj")
    }
    
    
    #-----------------------------------------------
    # Centroids
    #-----------------------------------------------
    if ("c" %in% names(layers))
      dat <- rbind(dat,
        plyr::ddply(
          .data      = data.frame(
            ord, 
            '.sample' = rownames(ord), 
            '.color'  = metadata(b, color.by)[rownames(ord)] ), 
          .variables = '.color', 
          .fun       = function (o) {
            data.frame(
              '.src'    = 'centroids',
              '.sample' = o[['.sample']],
              '.label'  = o[['.sample']],
              '.x'      = o[['.axis.1']],
              '.y'      = o[['.axis.2']],
              '.xend'   = o[['.axis.1']] %>% mean(),
              '.yend'   = o[['.axis.2']] %>% mean(),
              '.size'   = NA
            )
        })[colnames(dat)]
      )
      
    
    #-----------------------------------------------
    # BiPlot
    #-----------------------------------------------
    if (any(c("m", "a", "t") %in% names(layers))) {
      
      # Calculate Weighted average for each taxon
      biplot             <- distill(b, rank, safe = TRUE, md = FALSE)
      biplot[['.value']] <- biplot[['.value']] / sample.sums(biom)[biplot[['.sample']]]
      biplot[['.x']]     <- ord[biplot[['.sample']], '.axis.1']
      biplot[['.y']]     <- ord[biplot[['.sample']], '.axis.2']
      center.x           <- mean(biplot[['.x']])
      center.y           <- mean(biplot[['.y']])
      
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
      
      
      dat <- rbind(dat,
        data.frame(
          '.src'    = "biplot",
          '.sample' = NA,
          '.label'  = biplot[['.taxa']],
          '.x'      = biplot[['.x']],
          '.y'      = biplot[['.y']],
          '.xend'   = center.x,
          '.yend'   = center.y,
          '.size'   = biplot[['.value']]
        ))
      
      remove("biplot", "center.x", "center.y")
    }
    
    
    
    return (dat)
  })
  
  
  
  #-----------------------------------------------
  # Add in any missing metadata; preserve levels
  #-----------------------------------------------
  for (i in unique(c(color.by, shape.by, facet.by))) {
    
    if (i %in% names(ggdata)) {
      if (is.factor(biom[['metadata']][[i]])) {
        ggdata[[i]] <- factor(ggdata[[i]], levels = levels(biom[['metadata']][[i]]))
      } else {
        ggdata[[i]] <- as.factor(ggdata[[i]])
      }
      
    } else {
      ggdata[[i]] <- metadata(biom, i)[ggdata[['.sample']]] %>% as.factor()
    }
  }
  
  
  
  #-----------------------------------------------
  # How to fetch specific sub datasets from ggdata
  #-----------------------------------------------
  ord_data       <- function (x) subset(x, .src == "ord")
  stats_data     <- function (x) subset(x, .src == "stats")
  centroids_data <- function (x) subset(x, .src == "centroids")
  biplot_data    <- function (x) {
    x <- subset(x, .src == "biplot")
    x[['point.size']] <- x[['.size']]
    x[['label.size']] <- x[['.size']] / 20
    if (!"m" %in% names(layers)) x[['point.size']] <- 0
    return (x)
  }
  
  
  
  #-----------------------------------------------
  # Common aes arguments to subset as needed
  #-----------------------------------------------
  curr_aes <- function (...) {
    keys <- lapply(substitute(list(...))[-1], deparse)
    args <- list(
      x = ".x", xend = ".xend", size = ".size", 
      y = ".y", yend = ".yend", label = ".label" )
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
    layers[['p']][['args']] <- list(
      data    = ord_data,
      mapping = curr_aes(x, y, color, shape)
    )
  
  
  #-----------------------------------------------
  # Sample Names (IDs)
  #-----------------------------------------------
  if ("n" %in% names(layers))
    layers[['n']][['args']] <- list(
      data    = ord_data,
      mapping = curr_aes(x, y, color, label)
    )
  
  
  #-----------------------------------------------
  # Ellipses
  #-----------------------------------------------
  if ("e" %in% names(layers))
    layers[['e']][['args']] <- list(
      data    = ord_data,
      mapping = curr_aes(x, y, color)
    )
  
  
  #-----------------------------------------------
  # BiPlot Arrows
  #-----------------------------------------------
  if ("a" %in% names(layers))
    layers[['a']][['args']] <- list(
      'data'    = biplot_data,
      'mapping' = aes(x=.x, xend=.xend, y=.y, yend=.yend),
      'color'   = "darkgray", 
      'size'    = 0.75, 
      'alpha'   = 0.4,
      'arrow'   = arrow(ends="first", length=unit(.5,"cm"))
    )

  
  #-----------------------------------------------
  # BiPlot Means (Points)
  #-----------------------------------------------
  if ("m" %in% names(layers))
    layers[['m']][['args']] <- list(
      'data'    = biplot_data,
      'mapping' = aes(x=.x, y=.y, size=point.size),
      'color'   = "darkgray", 
      'alpha'   = 0.5
    )
  

  #-----------------------------------------------
  # BiPlot Taxa Names
  #-----------------------------------------------
  if ("t" %in% names(layers))
    layers[['t']][['args']] <- list(
      'data'    = biplot_data(ggdata),
      'mapping' = aes(
        'x'          = .x, 
        'y'          = .y, 
        'label'      = .label, 
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
  # Centroids Layer
  #-----------------------------------------------
  if (any(ggdata[['.src']] == "centroids"))
    layers[['c']][['args']] <- list (
      'data'    = centroids_data,
      'mapping' = curr_aes(x, y, xend, yend, color),
      'size'    = 0.75,
      'alpha'   = 0.4
    )
  
  
  #-----------------------------------------------
  # Stats Layer
  #-----------------------------------------------
  if (any(ggdata[['.src']] == "stats")) {
    
    if ("m" %in% anno) {
      elements[['labs']][['caption']] <- "Statistics computed with Adonis (1000 permutations)."
      elements[['theme']][['plot.caption']] <- element_text(face = "italic")
    }
    
    if (is.null(facet.by)) {
      elements[['labs']][['subtitle']] <- ggdata[ggdata[['.src']] == "adonis", '.label']
      
    } else {
      elements[['geom_text']] <- list(
        data    = stats_data,
        mapping = aes(label = .label),
        x       = -Inf, 
        y       =  Inf,
        hjust   = -0.04, 
        vjust   = 1.5
      )
      elements[['scale_y_continuous']][['expand']] <- expansion(mult = c(0, 0.08))
    }
  }
  
  
  #-----------------------------------------------
  # Colors, shapes, and sizes.
  #-----------------------------------------------
  
  if (!is.null(color.by)) {
    colors <- assign_colors(colors, ggdata[[color.by]])
    elements[['scale_color_manual']] <- list('values' = colors)
    elements[['scale_fill_manual']]  <- list('values' = colors)
  }
  
  if (!is.null(shape.by)) {
    shapes <- assign_shapes(shapes, ggdata[[shape.by]])
    elements[['scale_shape_manual']] <- list('values' = shapes)
  }
  
  if (any(c("m", "t") %in% names(layers)))
    elements[['scale_size_continuous']] <- list(
      'range'  = c(3, ifelse("m" %in% names(layers), 15, 5)),
      'name'   = "Taxa Abundance",
      'labels' = function (x) paste0(x * 100, "%") )
  
  
  
  #-----------------------------------------------
  # Load all layer datasets into ggplot
  #-----------------------------------------------
  p <- ggplot2::ggplot(ggdata) + theme_bw()
  
  
  
  #--------------------------------------------------------------
  # Facet arguments
  #--------------------------------------------------------------
  if (length(facet.by) == 1) {
    elements[['facet_wrap']] <- list(
      'facets' = backtick(facet.by),
      'scales' = "free")
    
  } else if (length(facet.by) > 1) {
    elements[['facet_grid']] <- list(
      'rows' = as.formula(paste(
        backtick(facet.by[[2]]), "~", backtick(facet.by[[1]]) )),
      'scales' = "free")
  }
  
  
  
  #-----------------------------------------------
  # Add customizable layers - points, lines, ...
  #-----------------------------------------------
  for (layer in names(layers)) {
    fun <- layers[[layer]][['fun']]
    
    # Unprefixed dot arguments, e.g. 'size'=2
    #--------------------------------------------------------------
    for (i in intersect(names(dots), formalArgs(fun)))
      layers[[layer]][['args']][[i]] <- dots[[i]]
    
    # Prefixed dot arguments, e.g. 'point.size'=2
    #--------------------------------------------------------------
    layer.short <- layer
    layer.title <- names(layerlist[layerlist == layer])
    for (qry in paste0(c(layer.short, layer.title), "."))
      for (i in names(dots)[startsWith(names(dots), qry)])
        layers[[layer]][['args']][[sub(qry, "", i)]] <- dots[[i]]
    
    p <- p + do.call(fun, layers[[layer]][['args']])
    
    remove("fun", "layer.short", "layer.title", "qry")
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

