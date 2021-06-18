

plot_factor <- function (biom, x, y, layers = "box", color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, colors = NULL, patterns = NULL, shapes = NULL, p.min = 0.05, p.adj = "fdr", vline = "ci95", xlab.angle = 'auto', ...) {
  
  dots <- list(...)
  
  #-----------------------------------------------
  # Which metadata columns will this plot need?
  #-----------------------------------------------
  md <- c()
  for (i in c("x", "color.by", "facet.by", "pattern.by", "shape.by")) {
    if (is.null(get(i))) next
    assign(i, validate_metrics(biom, get(i), mode='meta'))
    md <- c(md, get(i))
  }
  
  
  #--------------------------------------------------------------
  # Grouping variables for stats, etc.
  #--------------------------------------------------------------
  groupvars <- unique(c(x, color.by, pattern.by, shape.by, facet.by))
  
  if (length(groupvars) > 0) {
    biom$metadata[['.group']] <- interaction(
      lapply(groupvars, function (i) { biom$metadata[[i]] })
    )
    groupvars <- c('.group', groupvars)
    for (i in groupvars)
      biom$metadata[[i]] <- as.factor(biom$metadata[[i]])
  }
  
  
  #--------------------------------------------------------------
  # Convert biom object to a data.frame
  #--------------------------------------------------------------
  df <- distill(biom = biom, metric = y, md = unique(c(x, groupvars)))
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
  layerlist <- c(
    "violin", "box", "dot", "strip", 
    "crossbar", "errorbar", "linerange", "pointrange" )
  if (length(layers) == 1) layers <- strsplit(layers, "")[[1]]
  layers <- layerlist[pmatch(tolower(layers), layerlist)]
  
  
  #--------------------------------------------------------------
  # Merge sets of arguments. Preference given in following order:
  #  1. Prefixed `plot` arguments: plot(`box.fill`="red")
  #  2. `plot` arguments matching formalArgs: plot(`fill`="red")
  #  3. `withdots` arguments: withdots("box", `fill`="red")
  #--------------------------------------------------------------
  full2short <- c("facet_grid" = "facet", "facet_wrap" = "facet")
  
  short2full <- c("violin" = "geom_violin",   "crossbar"   = "geom_crossbar", 
                  "box"    = "geom_boxplot",  "errorbar"   = "geom_errorbar",
                  "dot"    = "geom_dotplot",  "linerange"  = "geom_linerange",
                  "strip"  = "geom_jitter",   "pointrange" = "geom_pointrange")
  
  short2regx <- c("violin" = "^v(|iolin)\\.", "crossbar"   = "^c(|rossbar)\\.", 
                  "box"    = "^b(|ox)\\.",    "errorbar"   = "^e(|rrorbar)\\.",
                  "dot"    = "^d(|ot)\\.",    "linerange"  = "^l(|inerange)\\.",
                  "strip"  = "^s(|trip)\\.",  "pointrange" = "^p(|ointrange)\\.",
                  "theme"  = "^t(|heme)\\.",  "facet"      = "^f(|acet)\\.")
  
  alldots <- function (fn, ...) {
    
    defaults <- list(...)
    
    full  <- ifelse(fn %in% names(short2full), short2full[[fn]], fn)
    short <- ifelse(fn %in% names(full2short), full2short[[fn]], fn)
    
    if (!is.null(defaults[['pattern']])) {
      func <- do.call(`::`, list("ggpattern", paste0(full, "_pattern")))
    } else if (full == "geom_dotplot") {
      func <- ggbeeswarm::geom_beeswarm
    } else {
      func <- do.call(`::`, list("ggplot2", full))
    }
    
    regx                 <- short2regx[[short]]
    specific_dots        <- dots[grep(regx, names(dots))]
    names(specific_dots) <- sub(regx, "", names(specific_dots), perl = TRUE)
    
    generic_dots <- dots[intersect(names(dots), formalArgs(func))]
    generic_dots <- generic_dots[setdiff(names(generic_dots), names(specific_dots))]
    args         <- c(generic_dots, specific_dots)
    args         <- c(args, defaults[setdiff(names(defaults), names(args))])
    
    do.call(func, args)
  }
  
  
  #--------------------------------------------------------------
  # Help multiple plots overlay well together
  #--------------------------------------------------------------
  dodge  <- ggplot2::position_dodge(width = 0.8)
  jdodge <- ggplot2::position_jitterdodge(dodge.width = 0.8, jitter.width = 0.05)
  
  # Violin arguments
  v_args <- list(fn = "violin", position = dodge, color = "black")
  if (!is.null(pattern.by)) v_args <- c(v_args, pattern = "magick")
  
  # Boxplot arguments
  b_args <- list(fn = "box", position = dodge, color = "black", width = 0.7)
  b_args[['width']] <- ifelse("violin" %in% layers, 0.1, 0.7)
  if ("violin" %in% layers)        {      b_args[['fill']]    <- "white"
  } else if (!is.null(pattern.by)) {      b_args[['pattern']] <- "magick"  }
  if (any(c("dot", "strip") %in% layers)) b_args[['outlier.shape']] <- NA
  
  # Dotplot arguments for ggplot2::geom_dotplot
  # d_args <- list(fn = "dot", position = dodge, binaxis = "y", stackdir = "center")
  # d_args[['binwidth']] <- diff(range(df[[y]])) / 30
  # if (any(c("violin", "box") %in% layers))
  #   d_args <- c(d_args, color = "black", fill = "black", dotsize = 0.7)
  
  # Dotplot arguments for ggbeeswarm::geom_beeswarm
  d_args <- list(fn = "dot", dodge.width = 0.8)
  if (any(c("violin", "box") %in% layers)) d_args <- c(d_args, color = "black", fill = "black")
  if (all(c("violin", "box") %in% layers)) d_args <- c(d_args, cex = 1, size = 0.8)
  
  # Stripchart arguments
  s_args <- list(fn = "strip",  position = jdodge)
  if (any(c("violin", "box") %in% layers))
    s_args <- c(s_args, color = "black")
  
  # Vertical line arguments
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
      mapping  = aes(y=y, ymin=ymin, ymax=ymax),
      position = dodge
    )
    
    if (any(c("violin", "box") %in% layers))
      args <- c(args, color = "black")
    
    c_args <- c(list(fn = "crossbar"),   args, width = 0.5)
    e_args <- c(list(fn = "errorbar"),   args, width = 0.5)
    l_args <- c(list(fn = "linerange"),  args)
    p_args <- c(list(fn = "pointrange"), args)
    remove("args")
    
    # crossbar fill
    if (any(c("violin", "box") %in% layers)) { c_args[['fill']]    <- NA
    } else if (!is.null(pattern.by))         { c_args[['pattern']] <- "magick"
    } else                                   { c_args[['fill']]    <- "white" }
  }
  
  
  # Facet arguments
  if (length(facet.by) == 1) {
    f_args <- list(
      'fn'     = "facet_wrap", 
      'facets' = backtick(facet.by) )
    
  } else if (length(facet.by) > 1) {
    f_args <- list(
      'fn'   = "facet_grid", 
      'rows' = as.formula(paste(
        backtick(facet.by[[2]]), "~", backtick(facet.by[[1]])
      )))
  }
  
  
  
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
  p <- ggplot(df, do.call(aes_string, aes_args, quote = TRUE))
  
  for (layer in intersect(layerlist, layers)) {
    args <- get(paste0(substr(layer, 1, 1), "_args"))
    p    <- p + do.call("alldots", args)
  }
  
  
  #--------------------------------------------------------------
  # Define which color, pattern, and shape names to use
  #--------------------------------------------------------------
  if (!is.null(color.by)) {
    colors <- assign_colors(colors, df[[color.by]])
    
    if (is.null(pattern.by)) {
      p <- p + ggplot2::scale_fill_manual(values = colors)
      p <- p + ggplot2::scale_color_manual(values = colors)
    } else {
      p <- p + ggpattern::scale_pattern_fill_manual(values = colors)
      p <- p + ggpattern::scale_pattern_color_manual(values = colors)
    }
  }
  
  if (!is.null(pattern.by)) {
    patterns <- assign_patterns(patterns, df[[pattern.by]])
    p        <- p + ggpattern::scale_pattern_type_manual(values = patterns)
  }
  
  if (!is.null(shape.by)) {
    shapes <- assign_shapes(shapes, df[[shape.by]])
    
    if (is.null(pattern.by)) {
      p <- p + ggplot2::scale_shape_manual(values = shapes)
    } else {
      p <- p + ggpattern::scale_pattern_shape_manual(values = shapes)
    }
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
      
      group.by <- unique(c(x, color.by, pattern.by, shape.by))
      groups   <- setNames(seq_along(levels(df[[group.by]])), levels(df[[group.by]]))
      pvals[['.xmin']] <- as.numeric(groups[pvals[['Group1']]])
      pvals[['.xmax']] <- as.numeric(groups[pvals[['Group2']]])
      remove("group.by", "groups")
      
      
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
  
  
  
  #--------------------------------------------------------------
  # Facet the plot
  #--------------------------------------------------------------
  if (!is.null(facet.by)) p <- p + do.call("alldots", f_args)
  
  
  #--------------------------------------------------------------
  # Arguments to pass to ggplot2::theme()
  #--------------------------------------------------------------
  t_args <- list(fn = "theme")
  
  # x-axis Text Angle
  if (xlab.angle == 'auto') {
    if (!is.null(x))
      if (sum(nchar(as.character(unique(df[[x]])))) > 40)
        xlab.angle <- 30
  }
  
  if (xlab.angle == 90) {
    t_args[['axis.text.x']] <- ggplot2::element_text(angle=-90, vjust=0.3, hjust=0)
    
  } else if (xlab.angle == 30) {
    t_args[['axis.text.x']] <- ggplot2::element_text(angle=-30, vjust=1, hjust=0)
    
    # Ensure long x-axis labels don't get truncated at the figure's edge
    rpad <- strwidth(tail(levels(df[[x]]), 1), units="inches")
    rpad <- rpad * 0.8660254 # sin((90-30) * pi / 180) / sin(90 * pi / 180)
    if (!is.null(color.by))
      rpad <- rpad - max(c(
        strwidth(color.by,               units="inches", cex=1.2) + .20,
        strwidth(levels(df[[color.by]]), units="inches")          + .52 ))
    t_args[['plot.margin']] <- ggplot2::unit(x=c(.1, max(.1, rpad), .1, .1), units="inches")
  }
  
  
  p <- p + do.call("alldots", t_args)
  
  attr(p, 'data')  <- df
  attr(p, 'stats') <- stats
  
  return (p)
  
}
