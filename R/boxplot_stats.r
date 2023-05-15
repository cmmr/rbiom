
#________________________________________________________
# Computes p-values for categorical differences
#________________________________________________________
boxplot_stats <- function (layers) {
  
  params      <- attr(layers, 'params', exact = TRUE)
  xcol        <- attr(layers, 'xcol', exact = TRUE)
  ycol        <- attr(layers, 'ycol', exact = TRUE)
  
  color.by    <- names(params[['color.by']])
  shape.by    <- names(params[['shape.by']])
  pattern.by  <- names(params[['pattern.by']])
  facet.by    <- params[['facet.by']]
  
  p.label     <- params[['p.label']]
  p.adj       <- params[['p.adj']]
  
  if (identical(p.label, TRUE))  p.label <- 0.05
  if (identical(p.label, FALSE)) p.label <- -Inf
  
  free_x <- free_y <- NULL
  if (hasLayer("facet")) {
    free_x <- isTRUE(layers[['facet']][['scales']] %in% c("free", "free_x"))
    free_y <- isTRUE(layers[['facet']][['scales']] %in% c("free", "free_y"))
  }
  
  
  y.pos <- unique(c(
    if (hasLayer("violin"))     "violin" else NULL,
    if (hasLayer("dot"))        "max"    else NULL,
    if (hasLayer("strip"))      "max"    else NULL,
    if (hasLayer("box"))        "box"    else NULL,
    if (hasLayer("errorbar"))   "vline"  else NULL,
    if (hasLayer("pointrange")) "vline"  else NULL,
    if (hasLayer("crossbar"))   "vline"  else NULL,
    if (hasLayer("linerange"))  "vline"  else NULL,
    if (hasLayer("bar"))        "mean"   else NULL
  ))
  
  stats  <- NULL
  ggdata <- attr(layers, 'data', exact = TRUE)
  
  if (isTRUE(nrow(ggdata) > 0 && is.numeric(p.label))) {
    
    #________________________________________________________
    # Significance values to put on top of the brackets
    #________________________________________________________
    p_annotations <- function (adj.p) {
      if (length(p.label) == 1)
        return (paste("italic(p)==", formatC(adj.p, format="g", digits = 1)))
      
      labels <- rep_len("", length(adj.p))
      for (i in p.label)
        labels <- paste0(labels, ifelse(adj.p <= i, "*", ""))
      return (labels)
    }
    
    
    if (length(unique(c(xcol, color.by, shape.by, pattern.by))) == 1) {
      if (isFALSE(xcol == ".taxa")) {
        
        #________________________________________________________
        # Brackets between x-values
        #________________________________________________________
        stats <- stats_table(
          biom        = ggdata, 
          x           = xcol, 
          y           = ycol, 
          by          = facet.by, 
          pairwise    = TRUE, 
          adj         = p.adj, 
          y.pos       = y.pos, 
          y.pos.facet = ".metric" ) # y.pos.facet )
        stats[['Group1']] %<>% factor(levels = levels(ggdata[[xcol]]))
        stats[['Group2']] %<>% factor(levels = levels(ggdata[[xcol]]))
        
        
        #________________________________________________________
        # Use p.top to retain only most significant taxa.
        #________________________________________________________
        if (isTRUE(is.finite(params[['p.top']]))) {
          
          taxa <- if (params[['p.top']] >= 1) {
            rank(plyr::daply(stats, '.taxa', function (s) {
              min(c(Inf, s[['.p.val']]), na.rm = TRUE) }), ties.method = "first")
            
          } else {
            plyr::daply(stats, '.taxa', function (s) {
              min(c(Inf, s[['.adj.p']]), na.rm = TRUE) })
          }
          
          taxa   <- names(taxa[taxa <= params[['p.top']]])
          stats  <- stats[stats[['.taxa']]   %in% taxa,,drop=FALSE]
          ggdata <- ggdata[ggdata[['.taxa']] %in% taxa,,drop=FALSE]
          stats[['.taxa']]  %<>% factor()
          ggdata[['.taxa']] %<>% factor()
          
          for (i in names(attributes(ggdata))) {
            obj <- attr(ggdata, i, exact = TRUE)
            if (!is.data.frame(obj))         next
            if (!".taxa" %in% colnames(obj)) next
            obj <- obj[obj[['.taxa']] %in% taxa,,drop=FALSE]
            obj[['.taxa']] %<>% factor()
            attr(ggdata, i) <- obj
          }
          
          remove("taxa")
        }
        
        
        #________________________________________________________
        # Don't display insignificant p-values on the plot
        #________________________________________________________
        pvals <- stats[which(stats[['.adj.p']] <= max(p.label)),,drop=FALSE]
        if (nrow(pvals) > 0) {
          
          xpos <- as.factor(ggdata[[xcol]])
          xpos <- setNames(seq_along(levels(xpos)), levels(xpos))
          pvals[['.xmin']]  <- as.numeric(xpos[as.character(pvals[['Group1']])])
          pvals[['.xmax']]  <- as.numeric(xpos[as.character(pvals[['Group2']])])
          pvals[['.label']] <- p_annotations(pvals[['.adj.p']])
          
          
          if (is_null(facet.by)) {
            pvals[['.step']] <- pvals[['y.pos']] * .13
            pvals[['y.pos']] <- pvals[['y.pos']] + (pvals[['.step']] * seq_len(nrow(pvals)))
            
          } else {
            
            # When some x-values are absent from some facets, we'll need
            # to re-number the remaining x-positions.
            
            plyby <- ply_cols(rev(facet.by))
            pvals <- plyr::ddply(pvals, plyby, function (z) {
                
              # Position bracket height on a per-facet basis
              z[['.step']] <- z[['y.pos']] * .13
              z[['y.pos']] <- z[['y.pos']] + (z[['.step']] * seq_len(nrow(z)))
              
              # Drop x categories that are absent from this facet
              if (free_x) {
                ggdata2      <- plyr::match_df(ggdata, z, on=facet.by)
                xpos2        <- intersect(names(xpos), ggdata2[[xcol]])
                xpos2        <- setNames(seq_along(xpos2), xpos2)
                z[['.xmin']] <- as.numeric(xpos2[as.character(z[['Group1']])])
                z[['.xmax']] <- as.numeric(xpos2[as.character(z[['Group2']])])
              }
              
              return (z)
            })
          }
          pvals[['.tick']] <- pvals[['y.pos']] - pvals[['.step']] / 4
          remove("xpos")
          
          pvals_df <- pvals[,intersect(names(pvals), names(ggdata)),drop=FALSE]
          
          attr(ggdata, 'stat_labels') <- data.frame(
            check.names = FALSE, 
            pvals_df,
            .x     = (pvals[['.xmin']] + pvals[['.xmax']]) / 2,
            .y     = pvals[['y.pos']] + pvals[['.step']] / 15,
            .label = pvals[['.label']] )
          
          
          if (nrow(pvals_df) == 0 || ncol(pvals_df) == 0) {
            pvals_df <- data.frame(row.names = seq_len(nrow(pvals) * 3))
          } else {
            pvals_df <- pvals_df %>% rbind(pvals_df) %>% rbind(pvals_df)
          }
          
          attr(ggdata, 'stat_brackets') <- data.frame(
            check.names = FALSE, 
            pvals_df,
            .x    = c(pvals[['.xmin']], pvals[['.xmin']], pvals[['.xmax']]),
            .xend = c(pvals[['.xmax']], pvals[['.xmin']], pvals[['.xmax']]),
            .y    = c(pvals[['y.pos']], pvals[['y.pos']], pvals[['y.pos']]),
            .yend = c(pvals[['y.pos']], pvals[['.tick']], pvals[['.tick']]) )
          
          # layers[['yaxis']][['limits']][[2]] <- max(pvals[['y.pos']]) * 1.1
          setLayer(
            'layer'        = "brackets",
            'mapping|x'    = ".x",
            'mapping|xend' = ".xend",
            'mapping|y'    = ".y",
            'mapping|yend' = ".yend" )
          
          setLayer(
            'layer'         = "stats_text",
            'size'          = 3,
            'vjust'         = ifelse(isTRUE(params[['flip']]), 0.5, 0),
            'parse'         = length(p.label) == 1,
            'mapping|x'     = ".x",
            'mapping|y'     = ".y",
            'mapping|label' = ".label" )
          
        }
      }
      
      
    } else {
      
      #________________________________________________________
      # P-value for each x position
      #________________________________________________________
      stats <- stats_table(
        biom        = ggdata, 
        x           = setdiff(c(color.by, pattern.by, shape.by), xcol), 
        y           = ycol, 
        by          = c(xcol, facet.by), 
        adj         = p.adj, 
        y.pos       = y.pos,
        y.pos.facet = ".metric" ) # y.pos.facet )
      
      
      #________________________________________________________
      # Use p.top to retain only most significant taxa.
      #________________________________________________________
      if (isTRUE(is.finite(params[['p.top']]))) {
        
        taxa <- if (params[['p.top']] >= 1) {
          rank(plyr::daply(stats, '.taxa', function (s) {
            min(c(Inf, s[['.p.val']]), na.rm = TRUE) }), ties.method = "first")
          
        } else {
          plyr::daply(stats, '.taxa', function (s) {
            min(c(Inf, s[['.adj.p']]), na.rm = TRUE) })
        }
        
        taxa   <- names(taxa[taxa <= params[['p.top']]])
        stats  <- stats[stats[['.taxa']]   %in% taxa,,drop=FALSE]
        ggdata <- ggdata[ggdata[['.taxa']] %in% taxa,,drop=FALSE]
        stats[['.taxa']]  %<>% factor()
        ggdata[['.taxa']] %<>% factor()
        
        for (i in names(attributes(ggdata))) {
          obj <- attr(ggdata, i, exact = TRUE)
          if (!is.data.frame(obj))         next
          if (!".taxa" %in% colnames(obj)) next
          obj <- obj[obj[['.taxa']] %in% taxa,,drop=FALSE]
          obj[['.taxa']] %<>% factor()
          attr(ggdata, i) <- obj
        }
        
        remove("taxa")
      }
      
      
      #________________________________________________________
      # Don't display insignificant p-values on the plot
      #________________________________________________________
      pvals <- stats[which(stats[['.adj.p']] <= max(p.label)),,drop=FALSE]
      if (nrow(pvals) > 0) {
      
        pvals[[xcol]] <- factor(pvals[[xcol]], levels = levels(ggdata[[xcol]]))
      
      
        # When some x-values are absent from some facets, we'll need
        # to re-number the remaining x-positions.
        
        if (!is_null(facet.by) && isTRUE(free_x)) {
          
          all_x <- levels(pvals[[xcol]])
          
          pvals <- plyr::ddply(
            .data      = pvals, 
            .variables = ply_cols(facet.by), 
            .fun       = function (z) {
                
              facet_x <- plyr::match_df(ggdata, z, on=facet.by)[[xcol]] %>%
                as.character() %>%
                unique()
              
              map_x <- sapply(all_x, function (x) which(facet_x == x) %>% ifelse(is_null(.), NA, .))
              
              z[['x.pos']] <- map_x[as.numeric(z[[xcol]])]
              
              
              return (z)
            })
          
        } else {
          pvals[['x.pos']] <- as.numeric(pvals[[xcol]])
        }
        
        
        attr(ggdata, 'stat_labels') <- data.frame(
          check.names = FALSE, 
          pvals,
          .x     = pvals[['x.pos']],
          .y     = pvals[['y.pos']] * 1.10,
          .label = p_annotations(pvals[['.adj.p']]) )
        
        if (!isTRUE(params[['flip']]))
          attr(ggdata, 'stat_brackets') <- data.frame(
            check.names = FALSE, 
            pvals,
            .x    = pvals[['x.pos']] - .4,
            .xend = pvals[['x.pos']] + .4,
            .y    = pvals[['y.pos']] * 1.08,
            .yend = pvals[['y.pos']] * 1.08 )
        
        if (!isTRUE(params[['flip']]))
          setLayer(
            'layer'        = "brackets",
            'mapping|x'    = ".x",
            'mapping|xend' = ".xend",
            'mapping|y'    = ".y",
            'mapping|yend' = ".yend" )
        
        setLayer(
          'layer'         = "stats_text",
          'size'          = 3,
          'hjust'         = ifelse(isTRUE(params[['flip']]), 0.1, 0.5),
          'vjust'         = ifelse(isTRUE(params[['flip']]), 0.5, 0),
          'parse'         = length(p.label) == 1,
          'mapping|x'     = ".x",
          'mapping|y'     = ".y",
          'mapping|label' = ".label" )
      
      
        setLayer("theme", "plot.caption" = element_text(face = "italic", color = "gray"))
        setLayer("labs",  "caption"      = local({
          
          test <- glue_collapse(
            x    = na.omit(unique(stats[['Test']])), 
            sep  = ", ", 
            last = " and " )
          
          meth <- switch(
            EXPR = p.adj,
            "holm"       = "FDR-adjusted as per Holm (1979)", 
            "hochberg"   = "FDR-adjusted as per Hochberg (1988)",
            "hommel"     = "FDR-adjusted as per Hommel (1988)",
            "bonferroni" = "FDR-adjusted using the Bonferroni correction",
            "BH"         = "FDR-adjusted as per Benjamini & Hochberg (1995)",
            "fdr"        = "FDR-adjusted as per Benjamini & Hochberg (1995)",
            "BY"         = "FDR-adjusted as per Benjamini & Yekutieli (2001)",
            "not corrected for multiple comparisons" )
          
          return(glue("{test} p-values {meth}."))
          
        }))
        
      }
    }
  }

  if (is_null(stats)) {
    setLayer("yaxis", expand = c(0.02, 0, 0.02, 0) )
    
  } else {
    
    if (isTRUE(params[['flip']])) {
      setLayer("yaxis", expand = c(0.02, 0, 0.15, 0) )
    } else {
      setLayer("yaxis", expand = c(0.02, 0, 0.08, 0) )
    }
    
    
    #________________________________________________________
    # Don't label the y-axis beyond the data range.
    #________________________________________________________
    if (hasLayer('facet')) {
      
      setLayer("yaxis", limits = c(0, NA))
      
    } else {
      
      setLayer(
        layer = "yaxis",
        breaks = local({
          ymax   <- min(stats[["y.pos"]])
          breaks <- labeling::extended(0, ymax, 5, only.loose = TRUE)
          return (breaks[breaks <= ymax])
      }))
      
      setLayer(
        'layer'   = "stats_bg",
        'color'   = NA,
        'fill'    = "white",
        'mapping' = aes_string(
          xmin = -Inf, 
          xmax = Inf, 
          ymin = signif(min(stats[['y.pos']]), 3), 
          ymax = Inf ))
    }
    
    
    dropcols <- c(".id", ".all", "x.pos", "y.pos")
    keepcols <- setdiff(names(stats), dropcols)
    stats    <- stats[, keepcols, drop=FALSE]
    attr(layers, "stats") <- stats
  }
  
  
  attr(layers, "data")  <- ggdata
  
  return (layers)
}
