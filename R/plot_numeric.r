
# See also ggpmisc package
# https://docs.r4photobiology.info/ggpmisc/articles/model-based-annotations.html


plot_numeric <- function (
  biom, x, y, layers = "pr", 
  color.by = NULL, shape.by = NULL, facet.by = NULL, 
  colors   = NULL, shapes   = NULL, 
  se = "ci95", model = y ~ x, regr = "lm", 
  p.adj = "fdr", anno = "mpr", ...) {
  
  dots <- list(...)
  
  
  
  #-----------------------------------------------
  # ggplot2 layers
  #-----------------------------------------------
  layerspecs <- list(
      'point'      = list('.fun' = ggplot2::geom_point,         '.regex' = "^p(|oint)\\."),
      'regression' = list('.fun' = ggplot2::stat_smooth,        '.regex' = "^r(|egression)\\."),
      'xaxis'      = list('.fun' = ggplot2::scale_x_continuous, '.regex' = "^x(|axis)\\."),
      'yaxis'      = list('.fun' = ggplot2::scale_y_continuous, '.regex' = "^y(|axis)\\."),
      'labs'       = list('.fun' = ggplot2::labs,               '.regex' = "^labs\\."),
      'theme'      = list('.fun' = ggplot2::theme,              '.regex' = "^t(|heme)\\."),
      'fill'       = list('.fun' = ggplot2::scale_fill_manual,  '.regex' = "^fill\\."),
      'color'      = list('.fun' = ggplot2::scale_color_manual, '.regex' = "^color\\."),
      'shape'      = list('.fun' = ggplot2::scale_shape_manual, '.regex' = "^shape\\."),
      'facet'      = list('.fun' = ggplot2::facet_wrap,         '.regex' = "^f(|acet)\\.") )
  
  
  layers <- layer_match(layers, c("p" = "point", "r" = "regression"), "pr")
  layers <- c(layers, 'xaxis', 'yaxis', 'labs', 'theme')
  
  if (!is.null(color.by)) layers <- c(layers, 'fill', 'color')
  if (!is.null(shape.by)) layers <- c(layers, 'shape')
  if (!is.null(facet.by)) layers <- c(layers, 'facet')
  
  layers <- layerspecs[layers]
  
  layers[['labs']][['y']] <- as.vector(y)
  
  
  #-----------------------------------------------
  # Add dot arguments to layers
  #-----------------------------------------------
  for (layer in names(layers)) {
    
    func <- layers[[layer]][['.fun']]
    regx <- layers[[layer]][['.regex']]
    
    # Unprefixed dot arguments, e.g. 'scales'="free_x"
    #--------------------------------------------------------------
    for (i in intersect(names(dots), formalArgs(func)))
      layers[[layer]][[i]] <- dots[[i]]
    
    # Prefixed dot arguments, e.g. 'facet.scales'="free_x"
    #--------------------------------------------------------------
    for (i in grep(regx, names(dots), value = TRUE))
      layers[[layer]][[sub(regx, "", i, perl = TRUE)]] <- dots[[i]]
  }
  remove(list = intersect(ls(), c("layer", "func", "regx", "i")))
  
  
  
  
  #-----------------------------------------------
  # The x, y, etc to plot
  #-----------------------------------------------
  for (i in unique(c(color.by, shape.by, facet.by)))
    biom[['metadata']][[i]] %<>% as.factor()
  
  md <- unique(c(x, color.by, shape.by, facet.by))
  df <- distill(biom = biom, metric = y, md = md, safe = TRUE)
  
  
  
  #-----------------------------------------------
  # Trend Line and Statistics
  #-----------------------------------------------
  if ("regression" %in% names(layers)) {
    
    
    # Parse the formula from 'model'
    #-----------------------------------------------
    tryCatch({
      formula <- if (!is(model, 'formula')) as.formula(model) else model
      stopifnot(is(formula, 'formula'))
    }, error = function (e) stop("Invalid value for `formula`. ", e) )
    
    
    # Parse the regression function from 'regr'
    #-----------------------------------------------
    tryCatch({
      
      if (is.character(regr)) {
        if (any(grepl("::", regr, fixed = TRUE))) { # regr = "stats::lm"
          method <- do.call(`::`, as.list(strsplit(regr, "::")[[1]]))
          
        } else { # regr = "lm"
          method <- get(regr)
        }
        
      } else { # regr = lm  or  regr = function(f,data) { ... }
        method <- regr
        regr   <- substr(deparse(substitute(regr)), 1, 10)
      }
      stopifnot(is(method, 'function'))
      
      if (!all(c('formula', 'data') %in% formalArgs(method)))
        stop("function does not have parameters 'formula' and 'data'.")
      
    }, error = function (e) stop("Invalid value for `regr`. ", e) )
    
    
    
    # Parameters for stat_smooth()
    #-----------------------------------------------
    layers[['regression']][['formula']] <- formula
    layers[['regression']][['method']]  <- method
    
    if (is.null(se)) {
      layers[['regression']][['se']] <- FALSE
      
    } else {
      se <- tolower(se)
      if (!grepl("^ci\\d+$", se)[[1]])
        stop("`se` must be of the form 'ci95', 'ci99', etc.")
      
      layers[['regression']][['se']]    <- TRUE
      layers[['regression']][['level']] <-  as.numeric(sub("ci", "", se)) / 100
    }
    
    
    # Run Statistics
    #-----------------------------------------------
    
    # change `y ~ x` into `Shannon ~ BMI`
    formula <- do.call(substitute, list(expr = model, env = list(x=as.name(x), y=as.name(y))))
    
    ply_by <- unique(intersect(names(df), c(color.by, shape.by, facet.by)))
    ply_df <- df[,intersect(names(df), c(".value", all.vars(formula), ply_by)),drop=F]
    names(ply_df) <- sub(".value", y, names(ply_df))
    
    # broomable <- methods(broom::glance) %>% as.character() %>% sub("glance.", "", .)
    
    glances <- plyr::ddply(ply_df, backtick(ply_by), function(x) {
      tryCatch(
        error   = function (e) return (NULL),
        warning = function (w) return (NULL),
        {
          method(formula, data = x) %>% broom::glance() %>% signif(3)
        })
    })
    if ('p.value' %in% names(glances))
      glances[['adj.p.value']] <- p.adjust(glances[['p.value']], method = p.adj)
    
    
    # Add Annotations
    #-----------------------------------------------
    anno <- layer_match(anno, c(names(glances), "method"), "mpr")
    
    if ('method' %in% anno)
      layers[['labs']][['caption']] <- paste0("Regression formula: ", regr, "(", format(formula), ")" )
    
    anno <- setdiff(anno, "method")
    if (isTRUE(length(anno) > 0 && nrow(glances) > 0)) {
      glances[['.label']] <- apply(glances[,anno,drop=F], 1L, function (x) paste(names(x), sep=" = ", x, collapse="; "))
      
      if (is.null(facet.by)) {
        layers[['labs']][['subtitle']] <- glances[['.label']]
        
      } else {
        layers[['stats_anno']] <- list(
          .fun    = ggplot2::geom_text,
          data    = glances[,c(ply_by, '.label')],
          mapping = aes(label = .label),
          x       = -Inf, 
          y       =  Inf,
          hjust   = -0.04, 
          vjust   = 1.5
        )
        layers[['yaxis']][['expand']] <- expansion(mult = c(0, 0.08))
      }
    }
  }
  
  
  #--------------------------------------------------------------
  # Facet arguments
  #--------------------------------------------------------------
  if (length(facet.by) == 1) {
    layers[['facet']][['facets']] = backtick(facet.by)
    
  } else if (length(facet.by) > 1) {
    layers[['facet']][['.fun']] <- ggplot2::facet_grid
    rows <- rows %>% head(2) %>% rev() %>% backtick()
    rows <- as.formula(paste(rows, collapse = " ~ "))
    layers[['facet']][['rows']] <- rows
    remove("rows")
  }
  
  
  #--------------------------------------------------------------
  # Set up the global aesthetics
  #--------------------------------------------------------------
  aes_args <- list(x = as.name(x), y = ".value")
  
  
  
  #--------------------------------------------------------------
  # Define which color and shape names to use
  #--------------------------------------------------------------
  if (!is.null(color.by)) {
    colors <- assign_colors(colors, df[[color.by]])
    aes_args[['color']] <- as.name(color.by)
    aes_args[['fill']]  <- as.name(color.by)
    layers[['fill']][['values']]  <- colors
    layers[['color']][['values']] <- colors
  }
  
  if (!is.null(shape.by)) {
    shapes <- assign_shapes(shapes, df[[shape.by]])
    aes_args[['shape']] <- as.name(shape.by)
    layers[['shape']][['values']] <- shapes
  }
  
  
  
  #--------------------------------------------------------------
  # Create the plot and add each layer with its arguments
  #--------------------------------------------------------------
  p <- ggplot(df, do.call(aes_string, aes_args, quote = TRUE)) + 
    theme_bw()
  
  for (layer in names(layers)) {
    
    func <- layers[[layer]][['.fun']]
    args <- layers[[layer]]
    args <- args[grep("^\\.", names(args), invert = TRUE)]
    
    p <- p + do.call(func, args)
    
  }
  
  
  # Attach table of stats
  if ("glances" %in% ls())
    attr(p, 'stats') <- glances
  
  
  return (p)
}
