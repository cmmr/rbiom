

plot_factor <- function (
  biom, x, y, layers = "rls", xvals = NULL,
  color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, 
  colors   = NULL, patterns   = NULL, shapes   = NULL, facets   = NULL, 
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
  
  for (i in groupvars)
    biom$metadata[[i]] <- as.factor(biom$metadata[[i]])
  
  
  #-----------------------------------------------
  # Layers to include in the plot
  #-----------------------------------------------
  layerlist <- c( 'b' = "box", 'v' = "violin",   'e' = "errorbar", 
                  'r' = "bar", 's' = "strip",    'l' = "linerange",
                  'd' = "dot", 'c' = "crossbar", 'p' = "pointrange")
  layers <- layer_match(layers, layerlist, "rls")
  
  if (!is.null(color.by))   layers %<>% c('fill', 'color')
  if (!is.null(shape.by))   layers %<>% c('shape')
  if (!is.null(facet.by))   layers %<>% c('facet')
  if (!is.null(pattern.by)) layers %<>% c('pattern')
  
  layers <- c('ggplot', layers, 'x_disc', 'y_cont', 'labs', 'theme', 'theme_bw')
  layers <- sapply(layers, function (x) list())
  
  remove("layerlist")
  
  
  #--------------------------------------------------------------
  # Convert biom values to relative abundances.
  #--------------------------------------------------------------
  if (is.rarefied(biom) && mode %in% c('rank', 'taxon'))
    biom[['counts']][['v']] <- biom[['counts']][['v']] / attr(biom, 'rarefaction')
  
  
  #--------------------------------------------------------------
  # Convert biom object to a data.frame
  #--------------------------------------------------------------
  ggdata <- distill(biom = biom, metric = y, md = opvars, safe=TRUE)
  names(ggdata) <- sub(attr(ggdata, 'response'), ".y", names(ggdata))
  ggdata[['.src']] <- "points"
  
  
  #--------------------------------------------------------------
  # Drop values not named in colors, patterns, or shapes.
  # Useful for filtering displayed combinations in beta div BPs.
  #
  # plot(hmp50, Bray ~ `==Sex`, color.by="Body Site", colors=c(
  #    'Saliva' = "green", 
  #    'Saliva vs Stool' = "red", 
  #    'Buccal mucosa vs Saliva' = "blue" ))
  #--------------------------------------------------------------
  ggdata <- local({
    conf <- list(
      "x"          = unname(xvals), 
      "facet.by"   = unname(facets), 
      "color.by"   = names(colors), 
      "pattern.by" = names(patterns), 
      "shape.by"   = names(shapes) )
    for (i in seq_along(conf)) {
      col  <- get(names(conf)[[i]])        %||% next
      keep <- conf[[i]]                    %||% next
      curr <- as.character(ggdata[[col]]) %||% next
      if (isTRUE(length(intersect(keep, curr)) > 0)) {
        ggdata        <- ggdata[curr %in% keep,,drop=F]
        ggdata[[col]] <- factor(curr[curr %in% keep], levels = keep)
      }
    }
    ggdata
  })
  
  
  #--------------------------------------------------------------
  # Add a group column to ensure groups are always separated
  #--------------------------------------------------------------
  if (length(groupvars) > 1) {
    group.by <- ".group"
    ggdata[[group.by]] <- interaction(
      lapply(groupvars, function (i) { ggdata[[i]] }) )
    groupvars <- c(groupvars, group.by)
    opvars    <- c(opvars,    group.by)
  } else {
    group.by <- groupvars
  }
  
  
  #--------------------------------------------------------------
  # Remove and report NA/NaN/Inf values
  #--------------------------------------------------------------
  err <- finite_check(ggdata)
  if (!is.null(err)) {
    ggdata <- ggdata[-err[['bad']],,drop=F]
    errMsg <- err[['msg']]
  }
  remove("err")
  
  
  #--------------------------------------------------------------
  # Add a column of alternate facet labels when `facets`= is used
  #--------------------------------------------------------------
  if (!is.null(facet.by) && !is.null(names(facets))) {
    if (!identical(names(facets), unname(facets))) {
      if (length(intersect(names(facets), as.character(ggdata[[facet.by]]))) > 0) {
        ggdata[['.facet']] <- facets[as.character(ggdata[[facet.by]])] %>%
          factor(levels = unname(facets))
        facet.by <- '.facet'
      }
    }
  }
  
  
  
  #--------------------------------------------------------------
  # Control the values displayed on the y axis
  #--------------------------------------------------------------
  if (is.rarefied(biom) && mode %in% c('rank', 'taxon')) {
    layers[['y_cont']][['labels']] <- ~ paste0(. * 100, "%")
  }
  layers[['y_cont']][['breaks']] <- base::pretty(ggdata[['.y']])
  layers[['y_cont']] %<>% c(list(
    'minor_breaks' = NULL,
    'limits'       = ggdata[['.y']] %>% 
      range(na.rm = TRUE, finite = TRUE) * c(0,1.1) ))
  
  
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
  
  layers[['labs']][['y']] <- ylab
  
  
  #--------------------------------------------------------------
  # Special case of "rank ~ ."
  #--------------------------------------------------------------
  if (isTRUE(x == ".taxa")) {
    layers[['labs']] %<>% c(list(x = NULL))
    
    ggdata[[x]] <- as.factor(ggdata[[x]])
    groupvars   <- setdiff(unique(c(groupvars, x)), ".group")
    
    if (length(groupvars) > 1) {
      group.by <- ".group"
      ggdata[[group.by]] <- interaction(
        lapply(groupvars, function (i) { ggdata[[i]] }) )
      groupvars <- c(groupvars, group.by)
    } else {
      group.by <- groupvars
    }
  }
  
  
  #--------------------------------------------------------------
  # Automatically facet by Metric or Taxa
  #--------------------------------------------------------------
  auto_facet <- attr(ggdata, 'facet', exact = TRUE)
  if (!is.null(auto_facet) && isFALSE(x == auto_facet)) {
    if (length(unique(ggdata[[auto_facet]])) > 1) {
      facet.by  <- head(unique(c(auto_facet, facet.by)), 2)
      groupvars <- unique(c(groupvars, facet.by))
      if (is.null(dots[['scales']])) dots[['scales']] <- "free_y"
    }
  }
  remove("auto_facet")
  
  
  #--------------------------------------------------------------
  # Help multiple plots overlay well together
  #--------------------------------------------------------------
  
  dodged    <- isTRUE(length(unique(c(x, color.by, pattern.by))) > 1)
  patterned <- isTRUE(!is.null(pattern.by))
  dotted    <- isTRUE(any(c("dot", "strip") %in% names(layers)))
  
  dodge  <- as.cmd(position_dodge(width = 0.8))
  jitter <- as.cmd(position_jitter(width = 0.25, height = 0, seed = 0))
  jdodge <- as.cmd(position_jitterdodge(
    dodge.width = 0.8, jitter.width = 0.05, jitter.height = 0, seed = 0) )
  
  
  
  # Violin arguments
  if ("violin" %in% names(layers)) {
    layers[['violin']] %<>% c(list(color = "black"))
    if (dodged)                   layers[['violin']][['position']] <- dodge
    if ("bar" %in% names(layers)) layers[['violin']][['alpha']]    <- 0.5
    if (patterned) {
      layers[['violin']][['pattern']] <- "magick"
      layers[['violin']][['fill']]    <- "white"
      layers[['violin']][['color']]   <- "black"
    }
  }
  
  # Bar plot arguments
  if ("bar" %in% names(layers)) {
    layers[['bar']] %<>% c(list(stat = "summary", fun = "mean"))
    if (dodged)             layers[['bar']][['position']] <- dodge
    if (dodged)             layers[['bar']][['width']]    <- 0.7
    if (length(layers) > 1) layers[['bar']][['alpha']]    <- 0.6
    if (patterned) {
      layers[['bar']][['pattern']] <- "magick"
      layers[['bar']][['fill']]    <- "white"
      layers[['bar']][['color']]   <- "black"
    }
    layers[['y_cont']][['breaks']] <- base::pretty(c(0, ggdata[['.y']]))
    
    # Bar charts need extra help on non-linear axes
    if (length(trans <- grep("\\.*trans$", names(dots), value=T)) == 1) {
      fun <- if (dots[[trans]] == 'sqrt')  { function (y) sqrt(mean(y * y))
      } else if (dots[[trans]] == 'log1p') { function (y) log1p(mean(exp(y) - 1))
      } else { stop('Bar charts can only be re-scaled using sqrt or log1p.') }
      layers[['bar']][['fun']] <- fun
      remove("fun")
    }
    remove("trans")
  }
  
  # Boxplot arguments
  if ("box" %in% names(layers)) {
    layers[['box']] %<>% c(list(color = "black", width = 0.7))
    
    if ("violin" %in% names(layers)) {
      layers[['box']][['fill']]  <- "white"
      layers[['box']][['width']] <- 0.1
      
    } else if (patterned)     {
      layers[['box']][['pattern']] <- "magick"
      layers[['box']][['fill']]    <- "white"
      layers[['box']][['color']]   <- "black"
      
    } else {
      layers[['box']][['alpha']] <- 0.6
    }
    
    if (dodged) layers[['box']][['position']]      <- dodge
    if (dotted) layers[['box']][['outlier.shape']] <- NA
  }
  
  # Dotplot arguments for ggbeeswarm::geom_beeswarm
  if ("dot" %in% names(layers)) {
    args <- list()
    if (dodged) args[['dodge.width']] <- 0.8
    if (any(c("violin", "box") %in% names(layers))) args <- c(args, color = "black", fill = "black")
    if (all(c("violin", "box") %in% names(layers))) args <- c(args, cex = 1, size = 0.8)
    layers[['dot']] %<>% c(args)
  }
  
  # Stripchart arguments for ggbeeswarm::geom_quasirandom
  if ("strip" %in% names(layers)) {
    args <- list('size' = 0.5)
    if (dodged) args[['dodge.width']] <- 0.8
    if (any(c("violin", "box", "bar") %in% names(layers)))
      args %<>% c(color = "black", fill = "black")
    layers[['strip']] %<>% c(args)
  }
  
  
  #--------------------------------------------------------------
  # Vertical line arguments
  #--------------------------------------------------------------
  if (any(c("crossbar", "errorbar", "linerange", "pointrange") %in% names(layers))) {
    
    if (substr(se, 1, 2) == "ci") {
      cl <- ifelse(se == "ci", "95", substr(se, 3, nchar(se)))
      cl <- as.numeric(cl) / 100
      vlineFn <- function (vals) {
        tt <- try(t.test(vals, conf.level = cl), silent = TRUE)
        if (!is(tt, "htest")) tt <- list(estimate = NA, conf.int = c(NA, NA))
        data.frame(
          .y    = ifelse(isTRUE(unname(tt$estimate) >= 0), unname(tt$estimate), 0),
          .ymin = ifelse(isTRUE(tt$conf.int[1]      >= 0), tt$conf.int[1],      0),
          .ymax = ifelse(isTRUE(tt$conf.int[2]      >= 0), tt$conf.int[2],      0)
        )
      }
    } else if (se == "range") {
      vlineFn <- function (vals) {
        data.frame(.y = mean(vals), .ymin = min(vals), .ymax = max(vals))
      }
    } else if (se == "mad") {
      vlineFn <- function (vals) {
        med <- median(vals); dev <- mad(vals, med)
        data.frame(.y = med, .ymin = med - dev, .ymax = med + dev)
      }
    } else if (se == "sd") {
      vlineFn <- function (vals) {
        avg <- mean(vals); dev <- sd(vals)
        data.frame(.y = avg, .ymin = avg - dev, .ymax = avg + dev)
      }
    } else if (se == "se") {
      vlineFn <- function (vals) {
        avg <- mean(vals); dev <- sqrt(var(vals)/length(vals))
        data.frame(.y = avg, .ymin = avg - dev, .ymax = avg + dev)
      }
    } else {
      stop("`se` must be one of 'ci95', 'mad', 'sd', 'se', or 'range'.")
    }
    
    
    #--------------------------------------------------------------
    # Add the ggdata rows for the vline function
    #--------------------------------------------------------------
    plyby <- NULL
    for (i in groupvars)
      plyby <- c(plyr::as.quoted(as.name(i)), plyby)
    
    vline_df <- plyr::ddply(ggdata, plyby, function (v) { vlineFn(v[[".y"]]) })
    vline_df[['.src']] <- "vline"
    ggdata %<>% append_df(vline_df)
    
    remove("plyby", "vline_df")
    
    
    #--------------------------------------------------------------
    # Additional arguments for the vline layer. `data` and
    # `mapping` parameters are handled by layers_toPlot().
    #--------------------------------------------------------------
    args <- list()
    
    if (dodged)                                     args[['position']] <- dodge
    if (patterned)                                  args[['color']]    <- "black"
    if (any(c("violin", "box") %in% names(layers))) args[['color']]    <- "black"
    
    if ("errorbar"   %in% names(layers)) layers[['errorbar']]   %<>% c(args)
    if ("linerange"  %in% names(layers)) layers[['linerange']]  %<>% c(args)
    if ("pointrange" %in% names(layers)) layers[['pointrange']] %<>% c(args)
    if ("crossbar"   %in% names(layers)) {
      if (any(c("violin", "box", "bar") %in% names(layers))) { args[['fill']]    <- NA
      } else if (patterned)                                  { args[['pattern']] <- "magick"
      } else                                                 { args[['fill']]    <- "white" }
      layers[['crossbar']] %<>% c(args)
    }
    remove("args")
  }
  
  
  #--------------------------------------------------------------
  # Facet arguments
  #--------------------------------------------------------------
  if (length(facet.by) == 1) {
    layers[['facet']]     <- list('facets' = backtick(facet.by))
    attr(dots, 'nfacets') <- length(unique(metadata(biom, facet.by)))
    
  } else if (length(facet.by) > 1) {
    layers[['facet']] <- list(
      'rows' = as.formula(paste(
        backtick(facet.by[[2]]), "~", backtick(facet.by[[1]])
      )))
    attr(dots, 'facet.nrow') <- length(unique(metadata(biom, facet.by[[2]])))
    attr(dots, 'facet.ncol') <- length(unique(metadata(biom, facet.by[[1]])))
  }
  
  
  
  #--------------------------------------------------------------
  # Aesthetic mappings for args not customized above
  #--------------------------------------------------------------
  mappings <- list(
    "box"        = c('x', 'y', 'group', 'color', 'fill', 'pattern_type', 'pattern_fill'),
    "violin"     = c('x', 'y', 'group', 'color', 'fill', 'pattern_type', 'pattern_fill'),
    "bar"        = c('x', 'y', 'group', 'color', 'fill', 'pattern_type', 'pattern_fill'),
    "dot"        = c('x', 'y', 'group', 'color', 'fill', 'shape'),
    "strip"      = c('x', 'y', 'group', 'color', 'fill', 'shape'),
    "crossbar"   = c('x',      'group', 'color', 'fill'),
    "pointrange" = c('x',      'group', 'color', 'fill'),
    "linerange"  = c('x',      'group', 'color'),
    "errorbar"   = c('x',      'group', 'color') )
  
  args <- list('x' = x, 'y' = ".y")
  if (!is.null(group.by))                         args[['group']]        <- group.by
  if (!is.null(color.by))                         args[['color']]        <- color.by
  if (!is.null(color.by))                         args[['fill']]         <- color.by
  if (!is.null(shape.by))                         args[['shape']]        <- shape.by
  if (!is.null(pattern.by))                       args[['pattern_type']] <- pattern.by
  if (!is.null(pattern.by) && !is.null(color.by)) args[['pattern_fill']] <- color.by
  
  for (layer in intersect(names(layers), names(mappings))) {
    for (arg in intersect(mappings[[layer]], names(args))) {
      if (is.null(layers[[layer]][[arg]])) {
        if (!"mapping" %in% names(layers[[layer]]))
          layers[[layer]][['mapping']] <- list()
        val <- args[[arg]]
        if (isTRUE(arg == "group" && val == ".all")) next
        val <- if (isTRUE(val == ".all")) NA else as.name(val)
        layers[[layer]][['mapping']][[arg]] <- val
      }
    }
  }
  remove("args")
  
  
  #--------------------------------------------------------------
  # Define which color, pattern, and shape names to use
  #--------------------------------------------------------------
  if (!is.null(color.by)) {
    ggdata[[color.by]] %<>% as.factor()
    colors <- assign_colors(colors, ggdata[[color.by]])
    layers[['fill']][['values']]  <- colors
    layers[['color']][['values']] <- colors
  }
  
  if (!is.null(pattern.by)) {
    ggdata[[pattern.by]] %<>% as.factor()
    patterns <- assign_patterns(patterns, ggdata[[pattern.by]])
    layers[['pattern']][['values']] <- patterns
  }
  
  if (!is.null(shape.by)) {
    ggdata[[shape.by]] %<>% as.factor()
    shapes <- assign_shapes(shapes, ggdata[[shape.by]])
    layers[['shape']][['values']] <- shapes
  }
  
  
  #--------------------------------------------------------------
  # Theme arguments
  #--------------------------------------------------------------
  args <- list()
  args[['text']] <- as.cmd(element_text('size' = 14))
  
  if (identical(as.vector(x), '.all')) {
    
    args[['axis.title.x']] <- as.cmd(element_blank())
    args[['axis.text.x']]  <- as.cmd(element_blank())
    args[['axis.ticks.x']] <- as.cmd(element_blank())
    
  } else {
    
    # x-axis Text Angle
    if (tolower(xlab.angle) == 'auto') {
      if (!is.null(x))
        if (sum(nchar(as.character(unique(ggdata[[x]])))) > 40)
          xlab.angle <- 30
    }
    
    if (xlab.angle == 90 || tolower(xlab.angle) == "vertical") {
      args[['axis.text.x']] <- as.cmd(element_text(angle=-90, vjust=0.3, hjust=0))
      
    } else if (xlab.angle == 30 || tolower(xlab.angle) == "angled") {
      args[['axis.text.x']] <- as.cmd(element_text(angle=-30, vjust=1, hjust=0))
      
      # Ensure long x-axis labels don't get truncated at the figure's edge
      rpad <- strwidth(tail(levels(ggdata[[x]]), 1), units="inches")
      rpad <- rpad * 0.8660254 # sin((90-30) * pi / 180) / sin(90 * pi / 180)
      if (!is.null(color.by))
        rpad <- rpad - max(c(
          strwidth(color.by, units="inches", cex=1.2) + .20,
          strwidth(levels(ggdata[[color.by]]), units="inches") + .52 ))
      rpad <- max(.1, signif(rpad, digits = 3))
      args[['plot.margin']] <- as.cmd(unit(x=c(.1, rpad, .1, .1), units='inches'), list(rpad=rpad))
    }
  }
  
  layers[['theme']] %<>% c(args)
  
  
  #--------------------------------------------------------------
  # Stat brackets
  #--------------------------------------------------------------
  y.pos <- ifelse(isTRUE("voilin" %in% names(layers)), "violin", "max")
  
  if (isTRUE(is.numeric(p.min) && p.min >= 0)) {
    if (length(unique(c(x, color.by, pattern.by, shape.by))) == 1) {
      if (isFALSE(x == ".taxa")) {
        
        #--------------------------------------------------------------
        # Brackets between x-values
        #--------------------------------------------------------------
        stats <- stats.table(ggdata, x, ".y", by = facet.by, pairwise = TRUE, adj = p.adj, y.pos = y.pos)
        stats[['Group1']] %<>% factor(levels = levels(ggdata[[x]]))
        stats[['Group2']] %<>% factor(levels = levels(ggdata[[x]]))
        
        pvals <- stats[stats[['adj.p']] <= max(p.min),]
        
        
        if (nrow(pvals) > 0) {
          
          xpos <- as.factor(ggdata[[x]])
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
          
          pvals_df <- pvals[,intersect(names(pvals), groupvars),drop=FALSE]
          
          ggdata %<>% append_df(data.frame(
            check.names = FALSE, 
            pvals_df,
            .src   = "stats",
            .x     = (pvals[['.xmin']] + pvals[['.xmax']]) / 2,
            .y     = pvals[['y.pos']] + ifelse(asterisks, 0, pvals[['.step']] / 5),
            .label = pvals[['.label']] ))
          
          
          if (nrow(pvals_df) == 0 || ncol(pvals_df) == 0) {
            pvals_df <- data.frame(row.names = seq_len(nrow(pvals) * 3))
          } else {
            pvals_df <- pvals_df %>% rbind(pvals_df) %>% rbind(pvals_df)
          }
          
          ggdata %<>% append_df(data.frame(
            check.names = FALSE, 
            pvals_df,
            .src  = "brackets",
            .x    = c(pvals[['Group1']], pvals[['Group1']], pvals[['Group2']]) %>% as.numeric(),
            .xend = c(pvals[['Group2']], pvals[['Group1']], pvals[['Group2']]) %>% as.numeric(),
            .y    = c(pvals[['y.pos']],  pvals[['y.pos']],  pvals[['y.pos']]),
            .yend = c(pvals[['y.pos']],  pvals[['.tick']],  pvals[['.tick']]) ))
          
          layers[['stats_text']] <- list('size' = 3, 'vjust' = 0, 'parse' = !asterisks)
          layers[['brackets']]   <- list()
          layers[['y_cont']][['limits']][[2]] <- max(pvals[['y.pos']]) * 1.1
          
        }
      }
      
      
    } else {
      
      #--------------------------------------------------------------
      # P-value for each x position
      #--------------------------------------------------------------
      stats <- stats.table(
        biom  = ggdata, 
        x     = setdiff(c(color.by, pattern.by, shape.by), x), 
        y     = ".y", 
        by    = c(x, facet.by), 
        adj   = p.adj, 
        y.pos = y.pos )
      
      stats[[x]] <- factor(stats[[x]], levels = levels(ggdata[[x]]))
      stats_df <- stats[,intersect(names(stats), groupvars),drop=FALSE]
      
      
      ggdata %<>% append_df(data.frame(
        check.names = FALSE, 
        stats_df,
        .src   = "stats",
        .x     = as.numeric(stats[[x]]),
        .y     = stats[['y.pos']] * 1.10,
        .label = paste("italic(p)==", stats[['adj.p']]) ))
      
      ggdata %<>% append_df(data.frame(
        check.names = FALSE, 
        stats_df,
        .src  = "brackets",
        .x    = as.numeric(stats[[x]]) - .4,
        .xend = as.numeric(stats[[x]]) + .4,
        .y    = stats[['y.pos']] * 1.08,
        .yend = stats[['y.pos']] * 1.08 ))
      
      
      layers[['stats_text']] <- list(size = 3, vjust = 0, parse = TRUE)
      layers[['brackets']]   <- list()
      layers[['y_cont']][['limits']][[2]] <- max(stats[['y.pos']]) * 1.2
      
    }
  }
  
  
  
  #--------------------------------------------------------------
  # One-off tweaks
  #--------------------------------------------------------------
  layers[['labs']][['y']] %<>% sub(
    fixed       = TRUE,
    pattern     = "OTUs Diversity", 
    replacement = "Observed OTUs" )
  
  
  #--------------------------------------------------------------
  # Create the plot and add each layer with its arguments.
  # Also attach a human-readable version of the plot command.
  #--------------------------------------------------------------
  layers[['ggplot']] %<>% c(list(data = ggdata))
  p <- layers_toPlot(layers, dots)
  
  
  #--------------------------------------------------------------
  # Remove extra columns that the user won't need
  #--------------------------------------------------------------
  attr(p, 'data') <- local({
    dropcols <- c(".depth", ".metric", ".all")
    if (isTRUE(length(unique(ggdata[['.src']])) == 1)) dropcols %<>% c(".src")
    ggdata[,setdiff(names(ggdata), dropcols),drop=FALSE]
  })
  
  if ("stats" %in% ls())
    attr(p, 'stats') <- local({
      dropcols <- c(".id", ".all", "y.pos")
      stats[,setdiff(names(stats), dropcols),drop=FALSE]
    })
  
  if ("errMsg" %in% ls())
    attr(p, 'err') <- errMsg
  
  
  return (p)
  
}
