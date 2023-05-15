#' Visualize BIOM data with boxplots.
#' 
#' @name bdiv_boxplot
#' 
#' @inherit adiv_boxplot params sections return
#' 
#' @family plotting
#' 
#' @param x   A categorical metadata column name. Prefix the column name with 
#'        \code{==} or \code{!=} to limit comparisons to within or between
#'        groups, respectively. The default, \code{NULL} groups all distances 
#'        into a single column.
#' 
#' @param metric   Beta diversity metric(s) to use. Options are 
#'        \code{"Manhattan"}, \code{"Euclidean"}, \code{"Bray-Curtis"}, 
#'        \code{"Jaccard"}, and \code{"UniFrac"}. UniFrac requires a 
#'        phylogenetic tree. Default: \code{"Bray-Curtis"}.
#'        
#' @param color.by,pattern.by,shape.by,facet.by,limit.by   Metadata columns to 
#'        use for data partitioning. Prefix the column name with 
#'        \code{==} or \code{!=} to limit comparisons to within or between
#'        groups, respectively. Default: \code{NULL}
#'        
#' @param weighted   When employing a beta diversity metric, use the weighted
#'        version. Default: \code{TRUE}.
#'        
#' @param tree   A phylogenetic tree for use in calculating UniFrac distance.
#'        The default, \code{NULL}, will use the BIOM object's tree.
#'        
#' @param filter.fun   A function that takes a \code{bdiv_table()} as input and
#'        returns a \code{bdiv_table()} as output. This modified table will be
#'        used for plotting. Default: \code{NULL}.
#'        
#' @param ...   Parameters are matched to formal arguments of ggplot2
#'        functions. Prefixing parameter names with a layer name ensures that
#'        a particular parameter is passed to, and only to, that layer. For
#'        instance, \code{dot.size = 2} or \code{d.size = 2} ensures only the 
#'        dotplot layer has its size set to \code{2}.
#' 
#' 
#' @export
#' @seealso \code{\link{stats_table}}
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     bdiv_boxplot(biom, x="==Body Site", metric="UniFrac", color.by="Body Site")
#'     
#'
bdiv_boxplot <- function (
  biom, x = NULL, metric = "Bray-Curtis", layers = "lsb",
  color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, limit.by = NULL, 
  p.adj = "fdr", p.label = TRUE, ci = 95, xlab.angle = 'auto', 
  weighted = TRUE, tree = NULL, ...) {
  
  with_cache(local({
  
    
    #________________________________________________________
    # Record the function call in a human-readable format.
    #________________________________________________________
    params <- c(as.list(environment()), list(...))
    params[['...']] <- NULL
    history <- attr(biom, 'history')
    history %<>% c(sprintf("bdiv_boxplot(%s)", as.args(params, fun = bdiv_boxplot)))
    remove(list = setdiff(ls(), c("params", "history")))
    
    
    #________________________________________________________
    # Sanity checks. x and *.by are checked by boxplot_build.
    #________________________________________________________
    params %<>% within({
      if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
      metric %<>% validate_arg(biom, 'metric', 'bdiv', n = c(1,Inf), tree = tree)
    })
    
    
    
    #________________________________________________________
    # Defaults
    #________________________________________________________
    if (is_null(params[['x']])) {
      params[['x']] <- ".all"
      params[['biom']][['metadata']][[".all"]] <- factor("all")
    }
    params[['weighted']] %<>% as.logical()
    params[['.cmp']]    <- TRUE
    
    
    
    #________________________________________________________
    # Use the generalized boxplot function to make the plot
    #________________________________________________________
    p <- boxplot_build(params, bdiv_boxplot, bdiv_boxplot_data, bdiv_boxplot_layers)
    
    attr(p, 'history') <- history
    
    
    return (p)
    
  }))
}


#______________________________________________________________
# Convert biom object to a data.frame
#______________________________________________________________
bdiv_boxplot_data <- function (params) {
  
  
  # Compute each bdiv metric separately - bray, etc
  #________________________________________________________
  ggdata <- plyr::ddply(
    .data      = expand.grid(
      '.metric' = params[['metric']], 
      '.weight' = params[['weighted']] ), 
    .variables = c(".metric", ".weight"), 
    .fun       =  function (x) {
      
      tbl <- bdiv_table(
        biom     = params[['biom']], 
        method   = x[1, '.metric'], 
        weighted = x[1, '.weight'],
        md       = paste0(params[['.md.cmps']], names(params[['.md.cmps']])),
        tree     = params[['tree']],
        safe     = TRUE )
      
      if (!is_null(params[['filter.fun']]))
        tbl <- params[['filter.fun']](tbl)
      
      return (tbl)
  })
  
  ggdata <- within(
    data = ggdata, 
    expr = {
      .y      <- .value
      .metric <- paste(ifelse(.weight, "Weighted", "Unweighted"), .metric)
      rm(".value", ".weight")
  })
  
  
  # Allow multiple `metric` beta div metrics
  #________________________________________________________
  if (length(params[['metric']]) > 1 || length(params[['weighted']]) > 1) {
    params[['facet.by']] %<>% c(".metric")
    ggdata[['.metric']]  %<>% factor(levels = unique(ggdata[['.metric']]))
  }
  
  
  attr(ggdata, 'params') <- params
  attr(ggdata, 'xcol')   <- params[['x']]
  attr(ggdata, 'ycol')   <- ".y"
  
  return (ggdata)
}



#______________________________________________________________
# Make bdiv-specific layer tweaks
#______________________________________________________________
bdiv_boxplot_layers <- function (layers) {
  
  # y-axis title
  #________________________________________________________
  yraw <- attr(layers, 'params', exact = TRUE)[['metric']]
  
  if (length(yraw) > 1) {
    ylab <- structure(
      .Data   = "Dissimilarity (\u03B2)", 
      display = '"Dissimilarity (\\u03B2)"' )
    
  } else {
    ylab <- paste(yraw, "Distance")
  }
  
  setLayer("labs", y = ylab)
  
  
  return (layers)
}


