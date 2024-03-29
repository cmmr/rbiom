
#______________________________________________________________
# Create the plot and add each layer with its arguments.
# Also attach a human-readable version of the plot command.
#______________________________________________________________
plot_build <- function (params) {
  
  
  ggdata <- params$.ggdata
  layers <- params$layers
  
  stopifnot(has_layer(params, 'ggplot'))
  stopifnot(is(ggdata, 'data.frame'))
  
  
  #______________________________________________________________
  # Return a placeholder plot if no data is present
  #______________________________________________________________
  if (nrow(ggdata) == 0) {
    p <- ggplot() + 
      geom_text(mapping = aes(x=1,y=1), label="No Data to Display") + 
      theme_void()
    p$plot_env <- emptyenv()
    attr(p, 'facet.nrow')  <- 1
    attr(p, 'facet.ncol')  <- 1
    attr(p, 'facet.count') <- 1
    return (p)
  }
  
  
  #________________________________________________________
  # Add in coord_flip when flip=TRUE
  #________________________________________________________
  if (isTRUE(params$flip)) add_layer(params, 'flip')
  
  
  #______________________________________________________________
  # Suppress vertical gridlines (horizontal if flipped).
  #______________________________________________________________
  if (isTRUE(params$flip)) {
    set_layer(params, 'theme', panel.grid.major.y = element_blank())
  } else {
    set_layer(params, 'theme', panel.grid.major.x = element_blank())
  }
  
  
  #______________________________________________________________
  # Standardize the list order: ggplot() first, theme() last, etc
  #______________________________________________________________
  layer_order <- c(
    'ggplot', 'stats_bg', 'stripe', 'violin', 'point', 'smooth', 
    'scatter', 'bar', 'box', 'spider', 'dot', 'ellipse', 'strip', 
    'name', 'crossbar', 'linerange', 'rect', 'errorbar', 'pointrange', 'mean', 
    'arrow', 'taxon', 'brackets', 'stats_text', 'stack', 'hline', 'vline', 
    'labs', 'color', 'fill', 'shape', 'pattern', 'size', 'continuous_scale', 
    'scale_size', 'facet', # 'free_y', 
    'xaxis', 'yaxis', 'flip', 'theme_bw', 'theme' )
  layer_order <- c(
    intersect(layer_order, env_names(layers)),
    setdiff(env_names(layers), layer_order) )
  
  
  #______________________________________________________________
  # Suppress x-axis labels when they're identical to facet labels
  #______________________________________________________________
  if (params$.xcol %in% params$facet.by) {
    
    if (isTRUE(params$flip)) {
      set_layer(
        params     = params, 
        layer      = 'theme', 
        .overwrite = TRUE,
        'axis.text.x'        = element_blank(),
        'axis.ticks.x'       = element_blank(),
        'axis.title.x'       = element_blank(),
        'panel.grid.major.x' = element_blank() )
      
    } else {
      set_layer(
        params     = params, 
        layer      = 'theme', 
        .overwrite = TRUE,
        'axis.text.y'        = element_blank(),
        'axis.ticks.y'       = element_blank(),
        'axis.title.y'       = element_blank(),
        'panel.grid.major.y' = element_blank() )
    }
  }
  
  
  
  #______________________________________________________________
  # x-axis and y-axis transformations
  #______________________________________________________________
  
  if (is_string(params$trans, c("rank", "log", "log1p", "sqrt")))
    params$.ylab %<>% sprintf("%s(%s)", params$trans, .)
  
  if (eq(params$trans, 'percent'))
    set_layer(params, 'yaxis', labels = scales::percent)
  
  
  for (axis in c('x', 'y')) {
    
    layer <- sprintf("%saxis", axis)
    label <- sprintf(".%slab", axis)
    trans <- params$layers[[layer]][['trans']]
    
    if (is.null(trans)) next
    
    
    # Set axis label to, e.g., `.ylab <- "Abundance (log scale)"`
    #______________________________________________________________
    params[[label]] %<>% sprintf("%s (%s scale)", ., trans)
    
    
    
    # Change 10000 to 10k for log scales.
    #______________________________________________________________
    if (trans %in% c("log", "log10", "log1p"))
      set_layer(params, layer, labels = label_number(scale_cut = cut_si("")))
    
    
    
    # Force sqrt scale to display zero tick mark.
    # Handle sqrt- and log1p- transforming of Inf and -Inf
    #______________________________________________________________
    if (eq(trans, "log1p") && isTRUE(params$stripe)) {
      
      params$layers[[layer]][['trans']] <- trans_new(
        name      = paste0(axis, "_log1p"), 
        transform = function (y) { y[is.finite(y)] <- base::log1p(y[is.finite(y)]); return (y); }, 
        inverse   = as.cmd(base::expm1) )
    
    } else if (eq(trans, "sqrt")) {
      
      if (isTRUE(params$stripe)) {
        
        params$layers[[layer]][['trans']] <- trans_new(
          name      = paste0(axis, "_sqrt"), 
          transform = function (y) { y[is.finite(y)] <- base::sqrt(y[is.finite(y)]); return (y); }, 
          inverse   = function (y) ifelse(y<0, 0, y^2) )
        
      } else {
        
        params$layers[[layer]][['trans']] <- trans_new(
          name      = paste0(axis, "_sqrt"), 
          transform = as.cmd(base::sqrt), 
          inverse   = function (y) ifelse(y<0, 0, y^2) )
      }  
      
    }
  }
  
  
  
  #______________________________________________________________
  # x-axis and y-axis labels
  #______________________________________________________________
  
  if (hasName(params, '.xlab')) set_layer(params, 'labs', x = params$.xlab)
  if (hasName(params, '.ylab')) set_layer(params, 'labs', y = params$.ylab)
  
  
  
  
  #______________________________________________________________
  # See if there's a better group column than '.group'
  #______________________________________________________________
  gcol <- '.group'
  if (hasName(ggdata, '.group')) {
    
    gvals <- ggdata[['.group']] %>% as.character() %>% factor() %>% as.numeric()
    gvals[is.na(gvals)] <- 0
    
    
    prefer <- lapply(env_names(params), function (i) {
      if (grepl("^(x|.+\\.by)$", i)) params[[i]] else NULL
    }) %>% unlist() %>% unname()
    mcols  <- grep("^\\.", colnames(ggdata), invert = TRUE, value = TRUE)
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
  
  
  #______________________________________________________________
  # Track which columns in `data` the plot command uses
  #______________________________________________________________
  mapped_cols <- params$facet.by
  
  
  #______________________________________________________________
  # Clean up args; remove vestigial columns from ggdata
  #______________________________________________________________
  
  for (layer in layer_order) {
    
    args <- layers[[layer]]
    fun  <- attr(args, 'function', exact = TRUE)
    src  <- attr(args, 'src',      exact = TRUE)
    
    
    if (is_null(src)) { data_cols <- colnames(ggdata)
    } else            { data_cols <- colnames(attr(ggdata, src, exact = TRUE)) }
    
    
    # Create the aes object for `mapping`=
    #______________________________________________________________
    if (hasName(args, 'mapping')) {
      
      if (is(args[['mapping']], "uneval")) {
        
        # mapping is already an aes() object.
        mapped_cols <- unique(c(mapped_cols, colnames(ggdata)))
        
      } else {
      
        aes_args <- args[['mapping']]
        aes_args <- aes_args[!sapply(aes_args, is.null)]
        
        # Put group last, so we can see if it's optional when it comes up
        aes_cols <- c()
        if ('group' %in% names(aes_args))
          aes_args <- aes_args[c(setdiff(names(aes_args), 'group'), 'group')]
        
        
        for (arg in names(aes_args)) {
          
          val <- as.vector(aes_args[[arg]])
          
          # if (eq(val, ".all"))   val <- NA
          if (eq(val, ".group")) val <- gcol
          
          
          if (is.character(val) && length(val) == 1) {
            
            if (val %in% data_cols) {
              
              if (arg == "group" && val %in% aes_cols) {
                
                # Avoid this: color = Sex, group = Sex
                val <- NULL
                
              } else {
              
                aes_cols <- unique(c(aes_cols, val))
                
                if (is_null(src))
                  mapped_cols <- unique(c(mapped_cols, val))
                
                # Add backticks to column names
                val <- as.name(val)
                #val <- sym(val)
                #val <- enquo(val)
              }
            
            }
            
          }
          
          aes_args[[arg]] <- val
        }
        
        
        if (isTRUE(length(aes_args) > 0)) {
          args[['mapping']] <- do.call(aes, aes_args)
          
        } else {
          args <- args[names(args) != 'mapping']
        }
        
        remove("aes_args", "aes_cols")
      }
    }
    
    
    # Don't specify an argument if it's already the default
    #______________________________________________________________
    defaults <- formals(fun)
    for (i in intersect(names(defaults), names(args)))
      if (eq(defaults[[i]], args[[i]]))
        args[[i]] <- NULL
    
    
    layers[[layer]] <- list(fun = fun, args = args)
  }
  
  # Drop columns while keeping attributes
  ggdata <- within(ggdata, rm(list = setdiff(names(ggdata), mapped_cols)))
  remove("mapped_cols")
  
  
  #______________________________________________________________
  # Assemble the ggplot object (p)
  #______________________________________________________________
  p    <- NULL
  cmds <- c("library(ggplot2)")
  
  for (layer in layer_order) {
    
    fun  <- layers[[layer]][['fun']]
    args <- layers[[layer]][['args']]
    fn   <- attr(fun, 'fn', exact = TRUE)
    
    
    # Skip theme() unless it has arguments
    #______________________________________________________________
    if (layer == 'theme' && length(args) == 0)
      next
    
    
    
    # Show ggplot() layer as "ggplot(data)", rest more verbosely
    #______________________________________________________________
    if (fn == "ggplot") {
      
      args[['data']] <- structure(ggdata, display = "data")
      p    <- do.call(fun, args) 
      cmds <- attr(p, 'display')
      
    } else {
      args[['.indent']] <- 4
      p    <- p + (q <- do.call(fun, args))
      cmds <- c(cmds, attr(q, 'display'))
    }
    
  }
  
  
  #______________________________________________________________
  # Attach code, facet info, etc
  #______________________________________________________________
  for (i in names(params$.plot_attrs))
    attr(p, i) <- params$.plot_attrs[[i]]
  
  attr(p, 'code') <- cmds %>% 
    paste(collapse=" +\n  ") %>% 
    add_class('rbiom_code')
  
  
  #________________________________________________________
  # Enable accessing attributes with `$`.
  #________________________________________________________
  p %<>% add_class('rbiom_plot')
  
  
  #______________________________________________________________
  # Prevent p from including all of this environment's variables.
  #______________________________________________________________
  p$plot_env <- emptyenv()
  
  
  return (p)
}
