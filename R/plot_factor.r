

plot_factor <- function (
  biom, x, y, layers = "rls", 
  color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, 
  colors   = NULL, patterns   = NULL, shapes   = NULL, 
  p.min = 0.05, p.adj = "fdr", se = "ci95", xlab.angle = 'auto', ...) {
  
  dots <- list(...)
  mode <- attr(y, 'mode', exact = TRUE)
  
  
  #-----------------------------------------------
  # Metadata needed for plotting / stat grouping
  #-----------------------------------------------
  opvars    <- c()
  groupvars <- c()
  group.by  <- NULL
  
  params <- c("color.by", "shape.by", "pattern.by", "facet.by")
  if (isTRUE(mode == "rank")) { x <- ".taxa"
  } else                      { params <- c("x", params) }
  
  for (param in params) {
    arg <- get(param)
    
    if (is.null(arg)) next
    if (is.null(attr(arg, 'mode', exact = TRUE))) {
      arg <- validate_metrics(biom, arg, mode='meta')
      assign(param, arg)
    }
    
    if (arg %in% groupvars) next
    groupvars <- c(groupvars, as.vector(arg))
    opvars    <- c(opvars, paste0(attr(arg, 'op', exact = TRUE), as.vector(arg)))
  }
  
  if (length(groupvars) > 1) {
    group.by <- ".group"
    biom$metadata[[group.by]] <- interaction(
      lapply(groupvars, function (i) { biom$metadata[[i]] }) )
    groupvars <- c(groupvars, group.by)
    opvars    <- c(opvars,    group.by)
  } else {
    group.by <- groupvars
  }
  
  for (i in groupvars)
    biom$metadata[[i]] <- as.factor(biom$metadata[[i]])
  
  
  #--------------------------------------------------------------
  # Settings for the ggplot object
  #--------------------------------------------------------------
  elements <- list()
  
  
  #--------------------------------------------------------------
  # Convert biom object to a data.frame
  #--------------------------------------------------------------
  df <- distill(biom = biom, metric = y, md = opvars, safe=TRUE)
  
  
  
  #--------------------------------------------------------------
  # Control the values displayed on the y axis
  #--------------------------------------------------------------
  if (is.rarefied(biom) && mode %in% c('rank', 'taxon')) {
    df[['.value']] <- df[['.value']] / depth(biom)
    breaks <- base::pretty(df[['.value']])
    labels <- paste0(breaks * 100, "%")
    
  } else {
    breaks <- base::pretty(df[['.value']])
    labels <- breaks
  }
  
  elements[['scale_y']] <- list(
    'limits'       = c(0,max(df[['.value']]) * 1.1),
    'breaks'       = breaks, 
    'labels'       = labels,
    'minor_breaks' = NULL )
  remove("breaks", "labels")
  
  
  #--------------------------------------------------------------
  # y-axis title
  #--------------------------------------------------------------
  ylab <- as.vector(y)
  
  if (mode == "adiv") {
    ylab <- paste(y, "Diversity")
  } else if (mode == "bdiv") {
    ylab <- paste(ifelse(dots[['weighted']], "Weighted", "Unweighted"), y, "Distance")
  } else if (mode %in% c("rank", "taxon")) {
    if (is.rarefied(biom)) { ylab <- paste(y, "\nRelative Abundance")
    } else                 { ylab <- paste(y, "\nRaw Abundance")      }
  }
  
  if (length(trans <- grep("\\.*trans$", names(dots), value=T)) == 1)
    ylab <- paste0(ylab, "\n(", dots[[trans]], " scale)")
  
  elements[['labs']][['y']] <- ylab
  
  
  #--------------------------------------------------------------
  # Special case of "rank ~ ."
  #--------------------------------------------------------------
  if (isTRUE(x == ".taxa")) {
    elements[['labs']] <- c(elements[['labs']], list(x = NULL))
    
    df[[x]]   <- as.factor(df[[x]])
    groupvars <- setdiff(unique(c(groupvars, x)), ".group")
    
    if (length(groupvars) > 1) {
      group.by <- ".group"
      df[[group.by]] <- interaction(
        lapply(groupvars, function (i) { df[[i]] }) )
      groupvars <- c(groupvars, group.by)
    } else {
      group.by <- groupvars
    }
  }
  
  
  #--------------------------------------------------------------
  # Automatically facet by Metric or Taxa
  #--------------------------------------------------------------
  auto_facet <- attr(df, 'facet', exact = TRUE)
  if (!is.null(auto_facet) && isFALSE(x == auto_facet)) {
    if (length(unique(df[[auto_facet]])) > 1) {
      facet.by  <- head(unique(c(auto_facet, facet.by)), 2)
      groupvars <- unique(c(groupvars, facet.by))
      if (is.null(dots[['scales']])) dots[['scales']] <- "free_y"
    }
  }
  remove("auto_facet")
  
  
  #--------------------------------------------------------------
  # Example `layers` argument: "box", c("BOX", "dot"), "V", "vbs"
  #--------------------------------------------------------------
  layers <- tolower(layers)
  layerlist <- c( 'b' = "box", 'v' = "violin",   'e' = "errorbar", 
                  'r' = "bar", 's' = "strip",    'l' = "linerange",
                  'd' = "dot", 'c' = "crossbar", 'p' = "pointrange")
  if (length(layers) == 1 && !layers %in% layerlist) {
    layers <- unname(layerlist[strsplit(layers, "")[[1]]])
  } else {
    layers <- unname(layerlist[pmatch(layers, layerlist)])
  }
  if (!isTRUE(length(layers) > 0)) layers <- "rls"
  remove("layerlist")
  
  
  #--------------------------------------------------------------
  # Interpret custom abbreviations for geom and scale names
  #--------------------------------------------------------------
  layer_geom <- function (layer) {
    
    if (layer == "facet")
      return (ifelse(length(facet.by) == 1, "facet_wrap", "facet_grid"))
    
    short2full <- c(
      'bar'        = "geom_bar",      'linerange'  = "geom_linerange",
      'violin'     = "geom_violin",   'pointrange' = "geom_pointrange",
      'strip'      = "geom_jitter",   'fill'       = "scale_fill_manual",
      'box'        = "geom_boxplot",  'color'      = "scale_color_manual",
      'dot'        = "geom_beeswarm", 'shape'      = "scale_shape_manual",
      'errorbar'   = "geom_errorbar", 'pattern'    = "scale_pattern_type_manual",
      'crossbar'   = "geom_crossbar", 'segment'    = "geom_segment",
      'text'       = "geom_text",     'scale_y'    = "scale_y_continuous")
    
    if (layer %in% names(short2full))
      return (short2full[[layer]])
    
    return (layer)
  }
  
  #--------------------------------------------------------------
  # How to match `dots` parameters to their appropriate layer
  #--------------------------------------------------------------
  layer_regex <- function (layer) {
    short2regex <- c(
      "violin" = "^v(|iolin)\\.", "crossbar"   = "^c(|rossbar)\\.", 
      "box"    = "^b(|ox)\\.",    "errorbar"   = "^e(|rrorbar)\\.",
      "dot"    = "^d(|ot)\\.",    "linerange"  = "^l(|inerange)\\.",
      "strip"  = "^s(|trip)\\.",  "pointrange" = "^p(|ointrange)\\.",
      "theme"  = "^t(|heme)\\.",  "facet"      = "^f(|acet)\\.",
      "bar"    = "^(|ba)r\\.",    "scale_y"    = "^(|scale_)y\\." )
    if (layer %in% names(short2regex))
      return (short2regex[[layer]])
    return (paste0("^", layer, "\\."))
  }
  
  
  #--------------------------------------------------------------
  # Return a ggpattern function only when necessary
  #--------------------------------------------------------------
  layer_func <- function (layer) {
    geom <- layer_geom(layer)
    func <- list('ggplot2', geom)
    if (layer == "dot")     func <- list("ggbeeswarm", geom)
    if (layer == "pattern") func <- list("ggpattern",  geom)
    if (!is.null(pattern.by)) {
      if (layer %in% c("bar", "box", "crossbar", "violin"))
        func <- list('ggpattern', paste0(geom, "_pattern"))
      if (layer %in% c("fill", "color", "shape"))
        func <- list('ggpattern', paste0("scale_pattern_", layer, "_manual"))
    }
    return (do.call(`::`, func))
  }
  
  
  #--------------------------------------------------------------
  # Help multiple plots overlay well together
  #--------------------------------------------------------------
  
  dodged    <- isTRUE(length(unique(c(x, color.by, pattern.by))) > 1)
  patterned <- isTRUE(!is.null(pattern.by))
  dotted    <- isTRUE(any(c("dot", "strip") %in% layers))
  
  dodge  <- ggplot2::position_dodge(width = 0.8)
  jitter <- ggplot2::position_jitter(width = 0.25, height = 0, seed = 0)
  jdodge <- ggplot2::position_jitterdodge(
    dodge.width = 0.8, jitter.width = 0.05, jitter.height = 0, seed = 0)
  
  # Violin arguments
  if ("violin" %in% layers) {
    elements[['violin']] <- list(color = "black")
    if (dodged)            elements[['violin']][['position']] <- dodge
    if (patterned)         elements[['violin']][['pattern']]  <- "magick"
    if ("bar" %in% layers) elements[['violin']][['alpha']]    <- 0.5
  }
  
  # Bar plot arguments
  if ("bar" %in% layers) {
    elements[['bar']] <- list(stat = "summary", fun = "mean")
    if (dodged)             elements[['bar']][['position']] <- dodge
    if (dodged)             elements[['bar']][['width']]    <- 0.7
    if (patterned)          elements[['bar']][['pattern']]  <- "magick"
    if (length(layers) > 1) elements[['bar']][['alpha']]    <- 0.6
    
    # Bar charts need extra help on non-linear axes
    if (length(trans <- grep("\\.*trans$", names(dots), value=T)) == 1) {
      fun <- if (dots[[trans]] == 'sqrt')  { function (y) sqrt(mean(y * y))
      } else if (dots[[trans]] == 'log1p') { function (y) log1p(mean(exp(y) - 1))
      } else { stop('Bar charts can only be re-scaled using sqrt or log1p.') }
      elements[['bar']][['fun']] <- fun
      remove("fun")
    }
    remove("trans")
  }
  
  # Boxplot arguments
  if ("box" %in% layers) {
    elements[['box']] <- list(color = "black", width = 0.7)
    
    if ("violin" %in% layers) {
      elements[['box']][['fill']]  <- "white"
      elements[['box']][['width']] <- 0.1
      
    } else if (patterned)     {
      elements[['box']][['pattern']] <- "magick"
      
    } else {
      elements[['box']][['alpha']] <- 0.6
    }
    
    if (dodged) elements[['box']][['position']]      <- dodge
    if (dotted) elements[['box']][['outlier.shape']] <- NA
  }
  
  # Dotplot arguments for ggbeeswarm::geom_beeswarm
  if ("dot" %in% layers) {
    args <- list()
    if (dodged) args[['dodge.width']] <- 0.8
    if (any(c("violin", "box") %in% layers)) args <- c(args, color = "black", fill = "black")
    if (all(c("violin", "box") %in% layers)) args <- c(args, cex = 1, size = 0.8)
    elements[['dot']] <- args
  }
  
  # Stripchart arguments
  if ("strip" %in% layers) {
    elements[['strip']][['position']] <- if (dodged) jdodge else jitter
    if (any(c("violin", "box", "bar") %in% layers)) {
      elements[['strip']][['color']] = "black"
      elements[['strip']][['size']]  = 0.1
    }
  }
  
  
  #--------------------------------------------------------------
  # Vertical line arguments
  #--------------------------------------------------------------
  if (any(c("crossbar", "errorbar", "linerange", "pointrange") %in% layers)) {
    
    if (substr(se, 1, 2) == "ci") {
      cl <- ifelse(se == "ci", "95", substr(se, 3, nchar(se)))
      cl <- as.numeric(cl) / 100
      vlineFn <- function (vals) {
        tt <- try(t.test(vals, conf.level = cl), silent = TRUE)
        if (!is(tt, "htest")) tt <- list(estimate = NA, conf.int = c(NA, NA))
        data.frame(
          y    = ifelse(isTRUE(unname(tt$estimate) >= 0), unname(tt$estimate), 0),
          ymin = ifelse(isTRUE(tt$conf.int[1]      >= 0), tt$conf.int[1],      0),
          ymax = ifelse(isTRUE(tt$conf.int[2]      >= 0), tt$conf.int[2],      0)
        )
      }
    } else if (se == "range") {
      vlineFn <- function (vals) {
        data.frame(y = mean(vals), ymin = min(vals), ymax = max(vals))
      }
    } else if (se == "mad") {
      vlineFn <- function (vals) {
        med <- median(vals); dev <- mad(vals, med)
        data.frame(y = med, ymin = med - dev, ymax = med + dev)
      }
    } else if (se == "sd") {
      vlineFn <- function (vals) {
        avg <- mean(vals); dev <- sd(vals)
        data.frame(y = avg, ymin = avg - dev, ymax = avg + dev)
      }
    } else if (se == "se") {
      vlineFn <- function (vals) {
        avg <- mean(vals); dev <- sqrt(var(vals)/length(vals))
        data.frame(y = avg, ymin = avg - dev, ymax = avg + dev)
      }
    } else {
      stop("`se` must be one of 'ci95', 'mad', 'sd', 'se', or 'range'.")
    }
    
    
    plyby <- NULL
    for (i in groupvars)
      plyby <- c(plyr::as.quoted(as.name(i)), plyby)
    
    args <- list(
      data     = plyr::ddply(df, plyby, function (v) { vlineFn(v[[".value"]]) }),
      mapping  = aes(y=y, ymin=pmax(ymin, 0), ymax=ymax)
    )
    if (dodged) args[['position']] <- dodge
    
    if (any(c("violin", "box") %in% layers))
      args <- c(args, color = "black")
    
    if ("errorbar"   %in% layers) elements[['errorbar']]   <- c(args, width = 0.5)
    if ("linerange"  %in% layers) elements[['linerange']]  <- args
    if ("pointrange" %in% layers) elements[['pointrange']] <- args
    if ("crossbar"   %in% layers) {
      if (any(c("violin", "box") %in% layers)) { args[['fill']]    <- NA
      } else if (!is.null(pattern.by))         { args[['pattern']] <- "magick"
      } else                                   { args[['fill']]    <- "white" }
      elements[['crossbar']] <- c(args, width = 0.5)
    }
  }
  
  
  #--------------------------------------------------------------
  # Facet arguments
  #--------------------------------------------------------------
  if (length(facet.by) == 1) {
    elements[['facet']] <- list('facets' = backtick(facet.by))
    
  } else if (length(facet.by) > 1) {
    elements[['facet']] <- list(
      'rows' = as.formula(paste(
        backtick(facet.by[[2]]), "~", backtick(facet.by[[1]])
      )))
  }
  
  
  #--------------------------------------------------------------
  # Theme arguments
  #--------------------------------------------------------------
  args <- list()
  
  if (identical(as.vector(x), '.all')) {
    
    args[['axis.title.x']] <- element_blank()
    args[['axis.text.x']]  <- element_blank()
    args[['axis.ticks.x']] <- element_blank()
    
  } else {
    
    # x-axis Text Angle
    if (xlab.angle == 'auto') {
      if (!is.null(x))
        if (sum(nchar(as.character(unique(df[[x]])))) > 40)
          xlab.angle <- 30
    }
    
    if (xlab.angle == 90) {
      args[['axis.text.x']] <- ggplot2::element_text(angle=-90, vjust=0.3, hjust=0)
      
    } else if (xlab.angle == 30) {
      args[['axis.text.x']] <- ggplot2::element_text(angle=-30, vjust=1, hjust=0)
      
      # Ensure long x-axis labels don't get truncated at the figure's edge
      rpad <- strwidth(tail(levels(df[[x]]), 1), units="inches")
      rpad <- rpad * 0.8660254 # sin((90-30) * pi / 180) / sin(90 * pi / 180)
      if (!is.null(color.by))
        rpad <- rpad - max(c(
          strwidth(color.by,               units="inches", cex=1.2) + .20,
          strwidth(levels(df[[color.by]]), units="inches")          + .52 ))
      args[['plot.margin']] <- ggplot2::unit(x=c(.1, max(.1, rpad), .1, .1), units="inches")
    }
  }
  
  elements[['theme']]   <- args
  
  
  
  #--------------------------------------------------------------
  # Set up the global aesthetics
  #--------------------------------------------------------------
  aes_args <- list(x = as.name(x), y = as.name(".value"))
  
  if (!is.null(group.by))
    aes_args[['group']] <- as.name(group.by)
  
  if (!is.null(pattern.by)) {
    aes_args[['pattern_type']] <- as.name(pattern.by)
    
    if (!is.null(color.by))
      aes_args[['pattern_fill']] <- as.name(color.by)
    
  } else if (!is.null(color.by)) {
    aes_args[['color']] <- as.name(color.by)
    aes_args[['fill']]  <- as.name(color.by)
  }
  
  if (!is.null(shape.by))
    aes_args[['shape']] <- as.name(shape.by)
  
  
  
  #--------------------------------------------------------------
  # Define which color, pattern, and shape names to use
  #--------------------------------------------------------------
  if (!is.null(color.by)) {
    colors <- assign_colors(colors, df[[color.by]])
    elements[['fill']][['values']]  <- colors
    elements[['color']][['values']] <- colors
  }
  
  if (!is.null(pattern.by)) {
    patterns <- assign_patterns(patterns, df[[pattern.by]])
    elements[['pattern']][['values']] <- patterns
  }
  
  if (!is.null(shape.by)) {
    shapes <- assign_shapes(shapes, df[[shape.by]])
    elements[['shape']][['values']] <- shapes
  }
  
  
  #--------------------------------------------------------------
  # Stat brackets
  #--------------------------------------------------------------
  stats <- NA
  
  y.pos <- ifelse(isTRUE("voilin" %in% layers), "violin", "max")
  
  if (length(unique(c(x, color.by, pattern.by, shape.by))) == 1) {
    if (isFALSE(x == ".taxa")) {
      
      #--------------------------------------------------------------
      # Brackets between x-values
      #--------------------------------------------------------------
      stats <- stats.table(df, x, ".value", by = facet.by, pairwise = TRUE, adj = p.adj, y.pos = y.pos)
      pvals <- stats[stats[['adj.p']] <= max(p.min),]
      
      
      if (nrow(pvals) > 0) {
        
        xpos <- as.factor(df[[x]])
        xpos <- setNames(seq_along(levels(xpos)), levels(xpos))
        pvals[['.xmin']] <- as.numeric(xpos[pvals[['Group1']]])
        pvals[['.xmax']] <- as.numeric(xpos[pvals[['Group2']]])
        remove("xpos")
        
        
        if (length(p.min) == 1) {
          asterisks <- FALSE
          pvals[['.label']] <- paste("italic(p)==", pvals[['adj.p']])
        } else {
          asterisks <- TRUE
          pvals[['.label']] <- ""
          for (i in p.min)
            pvals[['.label']] <- pvals[['.label']] %>%
            paste0(ifelse(pvals[['adj.p']] <= i, "*", ""))
        }
        
        
        if (is.null(facet.by)) {
          pvals[['.step']] <- pvals[['y.pos']] * .12
          pvals[['y.pos']] <- pvals[['y.pos']] + (pvals[['.step']] * seq_len(nrow(pvals)))
          
        } else {
          pvals <- plyr::ddply(
            .data      = pvals, 
            .variables = plyr::as.quoted(as.name(facet.by)), 
            .fun       = function (z) {
              z[['.step']] <- z[['y.pos']] * .12
              z[['y.pos']] <- z[['y.pos']] + (z[['.step']] * seq_len(nrow(z)))
              return (z)
            })
        }
        pvals[['.tick']] <- pvals[['y.pos']] - pvals[['.step']] / 4
        
        
        inherited_aes_cols <- data.frame(matrix(nrow=nrow(pvals), ncol=0))
        for (i in groupvars)
          inherited_aes_cols[[i]] <- if (is.null(pvals[[i]])) df[1,i] else pvals[[i]]
        
        
        elements[['text']] <- list(
          mapping = aes(x=.x, y=.y, label=.label),
          data    = data.frame(
            check.names = FALSE, 
            inherited_aes_cols,
            .x     = (pvals[['.xmin']] + pvals[['.xmax']]) / 2,
            .y     = pvals[['y.pos']] + ifelse(asterisks, 0, pvals[['.step']] / 5),
            .label = pvals[['.label']] ),
          color   = 'black',
          size    = 3,
          vjust   = 0,
          parse   = !asterisks )
        
        elements[['segment']] <- list(
          mapping = aes(x=.x, xend=.xend, y=.y, yend=.yend),
          data    = data.frame(
            check.names = FALSE, 
            inherited_aes_cols,
            .x    = c(pvals[['Group1']], pvals[['Group1']], pvals[['Group2']]),
            .xend = c(pvals[['Group2']], pvals[['Group1']], pvals[['Group2']]),
            .y    = c(pvals[['y.pos']],  pvals[['y.pos']],  pvals[['y.pos']]),
            .yend = c(pvals[['y.pos']],  pvals[['.tick']],  pvals[['.tick']]) ),
          color   = 'black' )
        
        elements[['scale_y']][['limits']][[2]] <- max(pvals[['y.pos']]) * 1.1
        
      }
    }
    
    
  } else {
    #--------------------------------------------------------------
    # P-value for each x position
    #--------------------------------------------------------------
    statcol <- setdiff(c(color.by, pattern.by, shape.by), x)
    if (length(statcol) > 1) {
      df[['.stat']] <- interaction(lapply(statcol, function (i) { df[[i]] }))
      statcol <- ".stat"
    }
    stats <- stats.table(df, statcol, ".value", by = c(x, facet.by), adj = p.adj, y.pos = y.pos)
    
    
    inherited_aes_cols <- data.frame(matrix(nrow=nrow(stats), ncol=0))
    for (i in groupvars)
      inherited_aes_cols[[i]] <- if (is.null(stats[[i]])) df[1,i] else stats[[i]]
    
    elements[['segment']] <- list(
      mapping = aes(x=.x, xend=.xend, y=.y, yend=.yend),
      data    = data.frame(
        check.names = FALSE, 
        inherited_aes_cols,
        .x    = as.numeric(stats[[x]]) - .4,
        .xend = as.numeric(stats[[x]]) + .4,
        .y    = stats[['y.pos']] * 1.08,
        .yend = stats[['y.pos']] * 1.08 ),
      color   = 'black' )
    
    elements[['text']] <- list(
      mapping = aes(x=.x, y=.y, label=.label),
      data    = data.frame(
        check.names = FALSE, 
        inherited_aes_cols,
        .x     = as.numeric(stats[[x]]),
        .y     = stats[['y.pos']] * 1.10,
        .label = paste("italic(p)==", stats[['adj.p']]) ),
      color   = 'black',
      size    = 3,
      vjust   = 0,
      parse   = TRUE )
    
    
    elements[['scale_y']][['limits']][[2]] <- max(stats[['y.pos']]) * 1.2
    
    remove("statcol")
  }
  
  
  #--------------------------------------------------------------
  # Create the plot and add each layer with its arguments
  #--------------------------------------------------------------
  p <- ggplot(df, do.call(aes_string, aes_args, quote = TRUE)) +
    theme_bw()
  
  
  for (layer in names(elements)) {
    func <- layer_func(layer)
    regx <- layer_regex(layer)
    
    # Unprefixed dot arguments, e.g. 'scales'="free_x"
    #--------------------------------------------------------------
    for (i in intersect(names(dots), formalArgs(func)))
      elements[[layer]][[i]] <- dots[[i]]
    
    # Prefixed dot arguments, e.g. 'facet.scales'="free_x"
    #--------------------------------------------------------------
    for (i in grep(regx, names(dots), value = TRUE))
      elements[[layer]][[sub(regx, "", i, perl = TRUE)]] <- dots[[i]]
      
    p <- p + do.call(func, elements[[layer]])
  }
  
  
  attr(p, 'data')  <- df
  attr(p, 'stats') <- stats
  
  return (p)
  
}
