
#--------------------------------------------------------------
# Create the plot and add each layer with its arguments.
# Also attach a human-readable version of the plot command.
#--------------------------------------------------------------
plot_build <- function (layers) {
  
  stopifnot('ggplot' %in% names(layers))
  
  ggdata <- attr(layers, 'data',   exact = TRUE)
  params <- attr(layers, 'params', exact = TRUE)
  stopifnot(is.data.frame(ggdata))
  
  
  #--------------------------------------------------------------
  # Return a placeholder plot if no data is present
  #--------------------------------------------------------------
  if (nrow(ggdata) == 0) {
    p <- ggplot() + 
      geom_text(mapping = aes(x=1,y=1), label="No Data to Display") + 
      theme_void()
    p$plot_env <- emptyenv()
    attr(p, 'facet.nrow') <- 1
    attr(p, 'facet.ncol') <- 1
    return (p)
  }
  
  
  #--------------------------------------------------------------
  # Standardize the list order: ggplot() first, theme() last, etc
  #--------------------------------------------------------------
  layer_order <- c(
    'ggplot', 'violin', 'point', 'bar', 'box', 'spider', 'dot', 'ellipse', 'strip', 
    'regression', 'name', 'crossbar', 'linerange', 'errorbar', 'pointrange', 
    'mean', 'arrow', 'taxon', 'brackets', 'stats_text', 'stack', 'labs', 
    'color', 'fill', 'shape', 'pattern', 'size', 'continuous_scale', 'scale_size', 
    'facet', 'xaxis', 'yaxis', 'theme_bw', 'theme' )
  layer_order <- c(
    intersect(layer_order, names(layers)),
    setdiff(names(layers), layer_order) )
  
  
  #--------------------------------------------------------------
  # Suppress x-axis labels when they're identical to facet labels
  #--------------------------------------------------------------
  if(attr(layers, 'xcol', exact = TRUE) %in% params[['facet.by']]) {
    layers[['theme']][['axis.text.x']]        <- element_blank()
    layers[['theme']][['axis.ticks.x']]       <- element_blank()
    layers[['theme']][['axis.title.x']]       <- element_blank()
    layers[['theme']][['panel.grid.major.x']] <- element_blank()
  }
  
  
  #--------------------------------------------------------------
  # Note any transformation in the y-axis label
  #--------------------------------------------------------------
  if (!is.null(layers[['labs']][['y']]))
    if (!is.null(layers[['yaxis']][['trans']]))
      layers[['labs']][['y']] %<>% paste0(" (", layers[['yaxis']][['trans']], " scale)")
  
  
  
  #--------------------------------------------------------------
  # See if there's a better group column than '.group'
  #--------------------------------------------------------------
  gcol <- '.group'
  if (hasName(ggdata, '.group')) {
    
    gvals <- ggdata[['.group']] %>% as.character() %>% factor() %>% as.numeric()
    gvals[is.na(gvals)] <- 0
    
    prefer <- params[grep("\\.by$", names(params))] %>% unlist() %>% unname()
    mcols  <- grep("^\\.", names(ggdata), invert = TRUE, value = TRUE)
    mcols  <- mcols[order(!mcols %in% prefer)]
    
    for (i in mcols) {
      if (!is.factor(ggdata[[i]])) next
      
      vals <- ggdata[[i]] %>% as.character() %>% factor() %>% as.numeric()
      vals[is.na(vals)] <- 0
      
      if (all(vals == gvals)) {
        gcol <- i
        break
      }
    }
    
    remove("gvals", "prefer", "mcols")
  }
  
  
  # #--------------------------------------------------------------
  # # When a layer only has a single row in ggdata, inline it.
  # #--------------------------------------------------------------
  # for (layer in layer_order) {
  #   
  #   args <- layers[[layer]]
  #   src  <- attr(args, 'src', exact = TRUE) %||% next
  #   
  #   row <- ggdata[ggdata[['.src']] == src,,drop=FALSE]
  #   if (nrow(row) != 1) next
  #   if (!'mapping' %in% names(args)) next
  #   if (!all(unname(args[['mapping']]) %in% names(row))) next
  #   
  #   for (i in seq_along(args[['mapping']])) {
  #     arg <- names(args[['mapping']][i])
  #     col <- unname(args[['mapping']][[i]])
  #     args[['mapping']][[arg]] <- row[[col]]
  #     
  #     remove("arg", "col")
  #   }
  #   
  #   ggdata <- ggdata[ggdata[['.src']] != src,,drop=FALSE]
  #   
  #   args[['data']]    <- as.cmd(data.frame())
  #   attr(args, 'src') <- NULL
  #   layers[[layer]]   <- args
  #   
  #   remove("args", "src", "row")
  # }
  
  
  #--------------------------------------------------------------
  # Track which columns in `data` the plot command uses
  #--------------------------------------------------------------
  mapped_cols <- params[['facet.by']]
  
  
  #--------------------------------------------------------------
  # Clean up args; remove vestigial columns from ggdata
  #--------------------------------------------------------------
  
  for (layer in layer_order) {
    
    args <- layers[[layer]]
    fn   <- attr(args, 'fn',       exact = TRUE)
    fun  <- attr(args, 'function', exact = TRUE)
    src  <- attr(args, 'src',      exact = TRUE)
    
    
    if (is.null(src)) { data_cols <- colnames(ggdata)
    } else            { data_cols <- colnames(attr(ggdata, src, exact = TRUE)) }
      
    
    
    # Create the aes object for `mapping`=
    #--------------------------------------------------------------
    if ('mapping' %in% names(args)) {
      
      aes_args <- args[['mapping']]
      aes_args <- aes_args[!sapply(aes_args, is.null)]
      
      # Put group last, so we can see if it's optional when it comes up
      aes_cols <- c()
      if ('group' %in% names(aes_args))
        aes_args <- aes_args[c(setdiff(names(aes_args), 'group'), 'group')]
      
      
      for (arg in names(aes_args)) {
        
        val <- tryCatch(expr = { as.vector(aes_args[[arg]]) }, error = function (e) browser() )
        # val <- as.vector(aes_args[[arg]])
        
        if (identical(val, ".all"))   val <- NA
        if (identical(val, ".group")) val <- gcol
        
        
        if (is.character(val) && length(val) == 1) {
          
          if (val %in% data_cols) {
            
            if (arg == "group" && val %in% aes_cols) {
              
              # Avoid this: color = Sex, group = Sex
              val <- NULL
              
            } else {
            
              aes_cols <- unique(c(aes_cols, val))
              
              if (is.null(src))
                mapped_cols <- unique(c(mapped_cols, val))
              
              # Add backticks to column names
              val <- capture.output(as.name(val))
            }
          
          } else {
            # Add double-quotes to other strings
            val <- glue::double_quote(val)
          }
        }
        
        aes_args[[arg]] <- val
      }
      
      if (isTRUE(length(aes_args) > 0)) {
        args[['mapping']] <- do.call(
          what  = aes_string, 
          args  = aes_args, 
          quote = TRUE )
        
      } else {
        args <- args[names(args) != 'mapping']
      }
      
      remove("aes_cols")
    }
    
    
    # Don't specify an argument if it's already the default
    #--------------------------------------------------------------
    defaults <- tryCatch(expr = { formals(fun) }, error = function (e) browser() )
    # defaults <- formals(fun)
    for (i in intersect(names(defaults), names(args)))
      if (identical(defaults[[i]], args[[i]]))
        args[[i]] <- NULL
    
    
    layers[[layer]] <- list(fn = fn, fun = fun, args = args)
  }
  
  # Drop columns while keeping attributes
  ggdata <- within(ggdata, rm(list = setdiff(names(ggdata), mapped_cols)))
  remove("mapped_cols")
  
  
  #--------------------------------------------------------------
  # Assemble the ggplot object (p)
  #--------------------------------------------------------------
  p    <- NULL
  cmds <- c("library(ggplot2)")
  
  for (layer in layer_order) {
    
    fn   = layers[[layer]][['fn']]
    fun  = layers[[layer]][['fun']]
    args = layers[[layer]][['args']]
    
    
    # Skip theme() unless it has arguments
    #--------------------------------------------------------------
    if (layer == "theme" && length(args) == 0)
      next
    
    
    # Force sqrt scale to display zero tick mark.
    #--------------------------------------------------------------
    if (layer == "yaxis" && identical(args[['trans']], "sqrt"))
      args[['trans']] <- as.cmd(scales::trans_new("sqrt0", base::sqrt, function(y) ifelse(y<0, 0, y^2)))
    
    
    # Show ggplot() layer as "ggplot(data)", rest more verbosely
    #--------------------------------------------------------------
    if (fn == "ggplot") {
      args[['data']] <- structure(ggdata, display = "data")
      p    <- do.call(fun, args) 
      cmds <- sprintf("%s(%s)", fn, as.args(args, fun=fun))
      
    } else {
      p    <- p + do.call(fun, args)
      cmds <- c(cmds, sprintf("%s(%s)", fn, as.args(args, indent = 4, fun=fun)))
    }
  }
  
  
  #--------------------------------------------------------------
  # Attach the number of facet rows and cols as plot attributes
  #--------------------------------------------------------------
  attr(p, 'facet.nrow') <- attr(layers, 'facet.nrow', exact = TRUE)
  attr(p, 'facet.ncol') <- attr(layers, 'facet.ncol', exact = TRUE)
  
  
  attr(p, 'cmd')   <- paste(collapse=" +\n  ", cmds)
  attr(p, 'data')  <- ggdata
  attr(p, 'stats') <- attr(layers, "stats", exact = TRUE)
  
  
  #--------------------------------------------------------------
  # Prevent p from including all of this environment's variables.
  #--------------------------------------------------------------
  p$plot_env <- emptyenv()
  
  
  return (p)
}