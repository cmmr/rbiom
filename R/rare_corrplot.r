#' Visualize rarefaction curves with scatterplots and trendlines.
#' 
#' @name rare_corrplot
#' 
#' @inherit adiv_corrplot params sections return
#' 
#' @family plotting
#' 
#' 
#' @param metric   Alpha diversity metric(s) to use. Options are: 
#'        \code{"OTUs"}, \code{"Shannon"}, \code{"Chao1"}, \code{"Simpson"}, 
#'        and/or \code{"InvSimpson"}. Default: \code{"OTUs"}.
#' 
#' @param depths   Rarefaction depths to show in the plot. Passed to
#'        \link{adiv_table}. The default, \code{"multi_even"}, uses a heuristic
#'        to pick 10 evenly spaced depths.
#'        
#' @param rline   Where to draw a horizontal line on the plot, intended to show
#'        a particular rarefaction depth. Set to \code{TRUE} to show an 
#'        auto-selected rarefaction depth or \code{NULL} to not show a line.
#'        Default: \code{NULL}.
#'        
#' @param caption   Display information about the method used for trendline
#'        fitting beneath the plot. Default: \code{TRUE}.
#'        
#' @param ...   Additional parameters to pass along to ggplot2 
#'        functions. Prefix a parameter name with either \code{p.}, 
#'        \code{r.}/\code{v.}, or \code{s.} to ensure it gets passed to (and
#'        only to)  \link[ggplot2]{geom_point}, \link[ggplot2]{vline}, or 
#'        \link[ggplot2]{geom_smooth}, respectively. For instance, 
#'        \code{p.size = 2} ensures only the points have their size set to 
#'        \code{2}.
#'        
#' @return A \code{ggplot2} plot. The computed data points and statistics will 
#'         be attached as \code{attr(p, 'data')} and \code{attr(p, 'stats')}, 
#'         respectively.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     rare_corrplot(hmp50, color.by="Body Site", metric=c("shannon", "otus"), facet.by="Sex") 
#'     

rare_corrplot <- function (
    biom, metric = "OTUs", depths = NULL, points = FALSE,
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    ci = 95, rline = NULL, caption = TRUE, ...) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(c(as.list(environment()), list(...)), eval)
  cache_file <- get_cache_file("rare_corrplot", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  arg_str <- as.args(params, fun = rare_corrplot, indent = 2)
  history <- paste0(collapse = "\n", c(
    attr(biom, 'history', exact = TRUE),
    sprintf("fig  <- rare_corrplot(%s)", arg_str) ))
  remove(list = setdiff(ls(), c("params", "history", "cache_file")))
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  params %<>% within({
    if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
    metric %<>% validate_arg(biom, 'metric', 'adiv', n = c(1,Inf))
    stopifnot(is_null(depths) || is.numeric(depths))
    stopifnot(is.logical(rline) || is_null(rline) || is.numeric(rline))
  })
  
  
  #________________________________________________________
  # Subset biom by requested metadata and aes.
  #________________________________________________________
  params %<>% metadata_params(contraints = list(
    color.by = list(n = 0:1),
    facet.by = list(col_type = "cat"),
    limit.by = list() ))
  
  
  #________________________________________________________
  # Default rarefaction depth.
  #________________________________________________________
  params[['rline']] <- with(params, {
    
    if (isFALSE(rline)) return (NULL)  # rline = FALSE
    if (!isTRUE(rline)) return (rline) # rline = 1380
    
    # rline = TRUE (auto-select rarefaction depth)
    ss    <- as.vector(sample_sums(biom))
    rline <- (sum(ss) * .1) / length(ss)
    rline <- min(ss[ss >= rline])
    return (rline)
  })
  
  
  #________________________________________________________
  # Pull rarefaction stats for each depth.
  #________________________________________________________
  if (is_null(params[['depths']]))
    params[['depths']] <- "multi_even"
  
  data <- adiv_table(
    biom    = params[['biom']],
    rarefy  = params[['depths']],
    metrics = params[['metric']],
    long    = TRUE, 
    md      = unique(c(
      params[['color.by']] %>% names(), 
      params[['facet.by']])) )
  
  
  #________________________________________________________
  # Facet by metric when there's more than one.
  #________________________________________________________
  if (length(params[['metric']]) > 1)
    params[['facet.by']] %<>% c(".metric")
  
  
  #________________________________________________________
  # Tell corrplot_build to use lm(y ~ log(x))
  #________________________________________________________
  params[['model']] <- "logarithmic"
  
  
  #________________________________________________________
  # Initialize the `layers` object.
  #________________________________________________________
  layers <- structure(
    list(),
    'data'     = data,
    'params'   = params,
    'function' = rare_corrplot,
    'xcol'     = ".depth",
    'ycol'     = attr(data, 'response', exact = TRUE),
    'xmode'    = "numeric" )
  
  initLayer(c('ggplot', 'smooth', 'labs', 'theme_bw'))
  
  if (!is_null(params[['color.by']])) initLayer("color")
  if (isTRUE(params[['points']]))     initLayer("point")
  
  
  #________________________________________________________
  # Add default layer parameters.
  #________________________________________________________
  setLayer("ggplot", mapping = list(x = ".depth", y = attr(data, 'response', exact = TRUE)))
  setLayer("xaxis",  labels = si_units)
  setLayer("labs",   x = "Rarefaction Depth", y = params[['metric']])
  
  if (length(params[['metric']]) > 1) {
    setLayer("labs", y = "Diversity")
    setLayer("facet", scales = "free_y")
  }
  
  if (!is_null(params[['rline']]))
    setLayer("vline", xintercept = params[['rline']], color = "red", linetype="dashed")
  
  
  
  # Convert layer definitions into a plot.
  #________________________________________________________
  p <- corrplot_build(layers)
  
  attr(p, 'history') <- history
  
  
  set_cache_value(cache_file, p)
  return (p)
}






