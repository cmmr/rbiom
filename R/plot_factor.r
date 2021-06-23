

plot_factor <- function (biom, x, y, layers = "box", color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, colors = NULL, patterns = NULL, shapes = NULL, p.min = 0.05, p.adj = "fdr", vline = "ci95", xlab.angle = 'auto', ...) {
  
  dots <- list(...)
  
  #-----------------------------------------------
  # Metadata needed for plotting / stat grouping
  #-----------------------------------------------
  opvars    <- c()
  groupvars <- c()
  for (param in c("x", "color.by", "shape.by", "pattern.by", "facet.by")) {
    arg <- get(param)
    
    if (is.null(arg)) next
    arg <- validate_metrics(biom, arg, mode='meta')
    assign(param, arg)
    
    if (arg %in% groupvars) next
    groupvars <- c(groupvars, as.vector(arg))
    opvars    <- c(opvars, paste0(attr(arg, 'op'), as.vector(arg)))
  }
  
  if (length(groupvars) > 0) {
    biom$metadata[['.group']] <- interaction(
      lapply(groupvars, function (i) { biom$metadata[[i]] }) )
    groupvars <- c(groupvars, '.group')
    opvars    <- c(opvars,    '.group')
    for (i in groupvars)
      biom$metadata[[i]] <- as.factor(biom$metadata[[i]])
  }
  
  
  #--------------------------------------------------------------
  # Convert biom object to a data.frame
  #--------------------------------------------------------------
  df <- distill(biom = biom, metric = y, md = opvars, safe=TRUE)
  y  <- attr(df, 'response')
  
  
  #--------------------------------------------------------------
  # Automatically facet by Metric or Taxa
  #--------------------------------------------------------------
  if (!is.null(attr(df, 'facet'))) {
    facet.by  <- head(unique(c(attr(df, 'facet'), facet.by)), 2)
    groupvars <- unique(c(groupvars, facet.by))
    if (is.null(dots[['scales']])) dots[['scales']] <- "free_y"
  }
  
  
  #--------------------------------------------------------------
  # Example `layers` argument: "box", c("BOX", "dot"), "V", "vbs"
  #--------------------------------------------------------------
  layers <- tolower(layers)
  layerlist <- c( 'b' = "box", 'v' = "violin",   'e' = "errorbar", 
                  'r' = "bar", 's' = "strip",    'l' = "linerange",
                  'd' = "dot", 'c' = "crossbar", 'p' = "pointrange" )
  if (length(layers) == 1 && !layers %in% layerlist) {
    layers <- unname(layerlist[strsplit(layers, "")[[1]]])
  } else {
    layers <- unname(layerlist[pmatch(layers, layerlist)])
  }
  if (!isTRUE(length(layers) > 0)) layers <- "box"
  remove("layerlist")
  
  
  #--------------------------------------------------------------
  # The layers that we're generally working with
  #--------------------------------------------------------------
  layer_geom <- function (layer) {
    
    if (layer == "facet")
      return (ifelse(length(facet.by) == 1, "facet_wrap", "facet_grid"))
    
    c("bar"        = "geom_bar",      "linerange"  = "geom_linerange",
      "violin"     = "geom_violin",   "pointrange" = "geom_pointrange",
      "strip"      = "geom_jitter",   "fill"       = "scale_fill_manual",
      "box"        = "geom_boxplot",  "color"      = "scale_color_manual",
      "dot"        = "geom_beeswarm", "shape"      = "scale_shape_manual",
      "errorbar"   = "geom_errorbar", "pattern"    = "scale_pattern_type_manual",
      "crossbar"   = "geom_crossbar", "theme"      = "theme" )[[layer]]
  }
  
  #--------------------------------------------------------------
  # How to match `dots` parameters to their appropriate layer
  #--------------------------------------------------------------
  layer_regex <- function (layer) c(
    "violin" = "^v(|iolin)\\.", "crossbar"   = "^c(|rossbar)\\.", 
    "box"    = "^b(|ox)\\.",    "errorbar"   = "^e(|rrorbar)\\.",
    "dot"    = "^d(|ot)\\.",    "linerange"  = "^l(|inerange)\\.",
    "strip"  = "^s(|trip)\\.",  "pointrange" = "^p(|ointrange)\\.",
    "theme"  = "^t(|heme)\\.",  "facet"      = "^f(|acet)\\.",
    "bar"    = "^(|ba)r\\." )[[layer]]
  
  
  #--------------------------------------------------------------
  # Return a ggpattern function only when necessary
  #--------------------------------------------------------------
  layer_func <- function (layer) {
    geom <- layer_geom(layer)
    x    <- list('ggplot2', geom)
    if (layer == "dot")     x <- list("ggbeeswarm", geom)
    if (layer == "pattern") x <- list("ggpattern",  geom)
    if (!is.null(pattern.by)) {
      if (layer %in% c("bar", "box", "crossbar", "violin"))
        x <- list('ggpattern', paste0(geom, "_pattern"))
      if (layer %in% c("fill", "color", "shape"))
        x <- list('ggpattern', paste0("scale_pattern_", layer, "_manual"))
    }
    return(do.call(`::`, x))
  }
  
  
  #--------------------------------------------------------------
  # Help multiple plots overlay well together
  #--------------------------------------------------------------
  elements <- list()
  
  dodged    <- isTRUE(length(unique(c(x, color.by, pattern.by))) > 1)
  patterned <- isTRUE(!is.null(pattern.by))
  dotted    <- isTRUE(any(c("dot", "strip") %in% layers))
  
  dodge  <- ggplot2::position_dodge(width = 0.8)
  jitter <- ggplot2::position_jitter(width = 0.25, seed=0)
  jdodge <- ggplot2::position_jitterdodge(dodge.width = 0.8, jitter.width = 0.05, seed = 0)
  
  # Violin arguments
  if ("violin" %in% layers) {
    elements[['violin']] <- list(color = "black")
    if (dodged)            elements[['violin']][['position']] = dodge
    if (patterned)         elements[['violin']][['pattern']]  = "magick"
    if ("bar" %in% layers) elements[['violin']][['alpha']]    = 0.5
  }
  
  # Bar plot arguments
  if ("bar" %in% layers) {
    elements[['bar']] <- list(stat = "summary", fun = "mean")
    if (dodged)             elements[['bar']][['position']] <- dodge
    if (dodged)             elements[['bar']][['width']]    <- 0.7
    if (length(layers) > 1) elements[['bar']][['alpha']]    <- 0.6
  }
  
  # Boxplot arguments
  if ("box" %in% layers) {
    elements[['box']] <- list(color = "black", width = 0.7)
    
    if ("violin" %in% layers) {
      elements[['box']][['fill']]  <- "white"
      elements[['box']][['width']] <- 0.1
      
    } else if (patterned)     {
      elements[['box']][['pattern']] <- "magick"
    }
    
    if (dodged) elements[['box']][['position']] <- dodge
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
    elements[['strip']] <- list()
    elements[['strip']][['position']] <- if (dodged) jdodge else jitter
    if (any(c("violin", "box") %in% layers))
      elements[['strip']][['color']] = "black"
  }
  
  
  #--------------------------------------------------------------
  # Vertical line arguments
  #--------------------------------------------------------------
  if (any(c("crossbar", "errorbar", "linerange", "pointrange") %in% layers)) {
    
    if (substr(vline, 1, 2) == "ci") {
      cl <- ifelse(vline == "ci", "95", substr(vline, 3, nchar(vline)))
      cl <- as.numeric(cl) / 100
      vlineFn <- function (vals) {
        tt <- try(t.test(vals, conf.level = cl), silent = TRUE)
        if (!is(tt, "htest")) tt <- list(estimate = NA, conf.int = c(NA, NA))
        data.frame(y = unname(tt$estimate), ymin = tt$conf.int[1], ymax = tt$conf.int[2])
      }
    } else if (vline == "range") {
      vlineFn <- function (vals) {
        data.frame(y = mean(vals), ymin = min(vals), ymax = max(vals))
      }
    } else if (vline == "mad") {
      vlineFn <- function (vals) {
        med <- median(vals); dev <- mad(vals, med)
        data.frame(y = med, ymin = med - dev, ymax = med + dev)
      }
    } else if (vline == "sd") {
      vlineFn <- function (vals) {
        avg <- mean(vals); dev <- sd(vals)
        data.frame(y = avg, ymin = avg - dev, ymax = avg + dev)
      }
    } else if (vline == "se") {
      vlineFn <- function (vals) {
        avg <- mean(vals); dev <- sqrt(var(vals)/length(vals))
        data.frame(y = avg, ymin = avg - dev, ymax = avg + dev)
      }
    } else {
      stop("vline must be one of 'ci95', 'mad', 'sd', or 'se'.")
    }
    
    
    plyby <- NULL
    for (i in groupvars)
      plyby <- c(plyr::as.quoted(as.name(i)), plyby)
    
    args <- list(
      data     = plyr::ddply(df, plyby, function (v) { vlineFn(v[[y]]) }),
      mapping  = aes(y=y, ymin=ymin, ymax=ymax)
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
  
  elements[['theme']] <- args
  
  
  
  #--------------------------------------------------------------
  # Set up the global aesthetics
  #--------------------------------------------------------------
  aes_args <- list(x = as.name(x), y = as.name(y), group = as.name(".group"))
  
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
  
  
  #--------------------------------------------------------------
  # Define which color, pattern, and shape names to use
  #--------------------------------------------------------------
  if (!is.null(color.by)) {
    colors <- assign_colors(colors, df[[color.by]])
    p <- p + layer_func("fill")(values = colors)
    p <- p + layer_func("color")(values = colors)
  }
  
  if (!is.null(pattern.by)) {
    patterns <- assign_patterns(patterns, df[[pattern.by]])
    p        <- p + layer_func("pattern")(values = patterns)
  }
  
  if (!is.null(shape.by)) {
    shapes <- assign_shapes(shapes, df[[shape.by]])
    p <- p + layer_func("shape")(values = shapes)
  }
  
  
  #--------------------------------------------------------------
  # Stat brackets
  #--------------------------------------------------------------
  stats <- NA
  
  y.pos <- ifelse(isTRUE("voilin" %in% layers), "violin", "max")
  
  if (length(unique(c(x, color.by, pattern.by, shape.by))) == 1) {
    #--------------------------------------------------------------
    # Brackets between x-values
    #--------------------------------------------------------------
    stats <- stats.table(df, x, y, by = facet.by, pairwise = TRUE, adj = p.adj, y.pos = y.pos)
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
      
      inherited_aes_cols <- data.frame(.group=rep(df[1,'.group'], nrow(pvals)))
      for (i in setdiff(groupvars, '.group'))
        inherited_aes_cols[[i]] <- if (is.null(pvals[[i]])) df[1,i] else pvals[[i]]
      
      
      p <- p + ggplot2::geom_text(
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
      
      p <- p + ggplot2::geom_segment(
        mapping = aes(x=.x, xend=.xend, y=.y, yend=.yend),
        data    = data.frame(
          check.names = FALSE, 
          inherited_aes_cols,
          .x    = pvals[['Group1']],
          .xend = pvals[['Group2']],
          .y    = pvals[['y.pos']],
          .yend = pvals[['y.pos']] ),
        color   = 'black' ) + 
        
        ggplot2::geom_segment(
          mapping = aes(x=.x, xend=.xend, y=.y, yend=.yend),
          data    = data.frame(
            check.names = FALSE, 
            inherited_aes_cols,
            .x    = pvals[['Group1']],
            .xend = pvals[['Group1']],
            .y    = pvals[['y.pos']],
            .yend = pvals[['y.pos']] - pvals[['.step']] / 4 ),
          color   = 'black' ) + 
        
        ggplot2::geom_segment(
          mapping = aes(x=.x, xend=.xend, y=.y, yend=.yend),
          data    = data.frame(
            check.names = FALSE, 
            inherited_aes_cols,
            .x    = pvals[['Group2']],
            .xend = pvals[['Group2']],
            .y    = pvals[['y.pos']],
            .yend = pvals[['y.pos']] - pvals[['.step']] / 4 ),
          color   = 'black' )
      
      
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
    stats <- stats.table(df, statcol, y, by = c(x, facet.by), adj = p.adj, y.pos = y.pos)
    
    inherited_aes_cols <- data.frame(.group=rep(df[1,'.group'], nrow(stats)))
    for (i in setdiff(groupvars, '.group'))
      inherited_aes_cols[[i]] <- if (is.null(stats[[i]])) df[1,i] else stats[[i]]
    
    p <- p + ggplot2::geom_segment(
      mapping = aes(x=.x, xend=.xend, y=.y, yend=.yend),
      data    = data.frame(
        check.names = FALSE, 
        inherited_aes_cols,
        .x    = as.numeric(stats[[x]]) - .4,
        .xend = as.numeric(stats[[x]]) + .4,
        .y    = stats[['y.pos']] * 1.08,
        .yend = stats[['y.pos']] * 1.08 ),
      color   = 'black' )
    
    p <- p + ggplot2::geom_text(
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
    
    p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .12)))
    
    remove("statcol")
  }
  
  
  
  
  attr(p, 'data')  <- df
  attr(p, 'stats') <- stats
  
  return (p)
  
}
