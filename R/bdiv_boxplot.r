#' Visualize BIOM data with boxplots.
#' 
#' @name bdiv_boxplot
#' 
#' @inherit adiv_boxplot params sections return
#' @inherit bdiv_distmat params
#' 
#' @family visualization
#' 
#' @param x   A categorical metadata column name. Prefix the column name with 
#'        \code{==} or \code{!=} to limit comparisons to within or between
#'        groups, respectively. The default, \code{NULL} groups all distances 
#'        into a single column.
#'        
#' @param color.by,pattern.by,shape.by,facet.by,limit.by   Metadata columns to 
#'        use for data partitioning. Prefix the column name with 
#'        \code{==} or \code{!=} to limit comparisons to within or between
#'        groups, respectively. Default: \code{NULL}
#'        
#' @param filter.fun   A function that takes a \code{bdiv_table()} as input and
#'        returns a \code{bdiv_table()} as output. This modified table will be
#'        used for plotting. Default: \code{NULL}.
#'        
#' @param ...   Parameters are matched to formal arguments of ggplot2
#'        functions. Prefixing parameter names with a layer name ensures that
#'        a particular parameter is passed to, and only to, that layer. For
#'        instance, \code{dot.size = 2} or \code{d.size = 2} ensures only the 
#'        dotplot layer has its size set to \code{2}. The special prefix
#'        \code{pt.} will control both the dot and strip layers.
#' 
#' 
#' @export
#' @seealso [biom_stats()]
#' @examples
#'     library(rbiom)
#'     
#'     biom <- sample_rarefy(hmp50)
#'     bdiv_boxplot(biom, x="==Body Site", bdiv="UniFrac", color.by="Body Site")
#'     
#'
bdiv_boxplot <- function (
  biom, x = NULL, bdiv = "Bray-Curtis", layers = "lsb",
  color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, limit.by = NULL, 
  flip = FALSE, stripe = flip, p.adj = "fdr", p.label = TRUE, ci = 95, outliers = NULL, 
  xlab.angle = 'auto', weighted = TRUE, tree = NULL, ...) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(c(as.list(environment()), list(...)), eval)
  cache_file <- get_cache_file("bdiv_boxplot", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  arg_str <- as.args(params, fun = bdiv_boxplot, indent = 2)
  history <- paste0(collapse = "\n", c(
    attr(biom, 'history', exact = TRUE),
    sprintf("fig  <- bdiv_boxplot(%s)", arg_str) ))
  remove(list = setdiff(ls(), c("params", "history", "cache_file")))
  
  
  #________________________________________________________
  # Sanity checks. x and *.by are checked by boxplot_build.
  #________________________________________________________
  params %<>% within({
    if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
    bdiv %<>% validate_arg(biom, 'bdiv', 'bdiv', n = c(1,Inf), tree = tree)
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
  
  
  set_cache_value(cache_file, p)
  return (p)
}


#______________________________________________________________
# Convert biom object to a data.frame
#______________________________________________________________
bdiv_boxplot_data <- function (params) {
  
  
  # Compute each bdiv distance metric separately: bray, etc
  #________________________________________________________
  ggdata <- plyr::ddply(
    .data      = expand.grid(
      '.bdiv'     = params[['bdiv']], 
      '.weighted' = params[['weighted']] ), 
    .variables = c(".bdiv", ".weighted"), 
    .fun       =  function (x) {
      
      tbl <- bdiv_table(
        biom     = params[['biom']], 
        bdiv     = as.character(x[1, '.bdiv']), 
        weighted = x[1, '.weighted'],
        md       = paste0(params[['.md.cmps']], names(params[['.md.cmps']])),
        tree     = params[['tree']] )
      
      if (!is_null(params[['filter.fun']]))
        tbl <- params[['filter.fun']](tbl)
      
      return (tbl)
  })
  
  ggdata %<>% within({
    .bdiv <- paste(ifelse(.weighted, "Weighted", "Unweighted"), .bdiv) })
  
  
  # Allow multiple bdiv metrics
  #________________________________________________________
  if (length(params[['bdiv']]) > 1 || length(params[['weighted']]) > 1) {
    params[['facet.by']] %<>% c(".bdiv")
    ggdata[['.bdiv']]  %<>% factor(levels = unique(ggdata[['.bdiv']]))
  }
  
  
  attr(ggdata, 'params')   <- params
  attr(ggdata, 'xcol')     <- params[['x']]
  attr(ggdata, 'response') <- ".distance"
  attr(ggdata, 'ycol')     <- ".distance"
  
  return (ggdata)
}



#______________________________________________________________
# Make bdiv-specific layer tweaks
#______________________________________________________________
bdiv_boxplot_layers <- function (layers) {
  
  # y-axis title
  #________________________________________________________
  yraw <- attr(layers, 'params', exact = TRUE)[['bdiv']]
  
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


