#' Visualize rarefaction curves with scatterplots and trendlines.
#' 
#' @name rare_corrplot
#' 
#' @inherit adiv_corrplot params sections return
#' 
#' @family visualization
#' 
#' 
#' @param depths   Rarefaction depths to show in the plot. Passed to
#'        [adiv_table()]. The default, \code{"multi_even"}, uses a heuristic
#'        to pick 10 evenly spaced depths.
#'        
#' @param rline   Where to draw a horizontal line on the plot, intended to show
#'        a particular rarefaction depth. Set to \code{TRUE} to show an 
#'        auto-selected rarefaction depth or \code{NULL} to not show a line.
#'        Default: \code{NULL}.
#'        
#' @param model   What type of trendline to fit to the data. Options are: 
#'        \code{"lm"} (linear), \code{"log"} (logarithmic), or \code{"gam"} 
#'        (generalized additive). You can alternatively provide a list of 
#'        length two where the first element is a character vector of length 1 
#'        naming a function, and the second element is a list of arguments to 
#'        pass to that function. One of the function's arguments must be named 
#'        'formula'. For example, 
#'        \code{model = list("stats::lm", list(formula = y ~ x))}.
#'        Default: \code{"log"}.
#'        
#' @param ...   Additional parameters to pass along to ggplot2 
#'        functions. Prefix a parameter name with either \code{p.}/\code{pt.}, 
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
#'     rare_corrplot(hmp50, color.by="Body Site", adiv=c("shannon", "otus"), facet.by="Sex") 
#'     

rare_corrplot <- function (
    biom, adiv = "Shannon", depths = NULL, layers = "t", rline = TRUE,
    color.by = NULL, facet.by = NULL, limit.by = NULL, model = "log", 
    stats = "emtrends", p.adj = "fdr", ci = 95, caption = FALSE, ...) {
  
  
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
    adiv %<>% validate_arg(biom, 'adiv', 'adiv', n = c(1,Inf))
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
    biom   = params[['biom']],
    rarefy = params[['depths']],
    adiv   = params[['adiv']],
    long   = TRUE, 
    md     = unique(c(
      params[['color.by']] %>% names(), 
      params[['facet.by']])) )
  
  
  #________________________________________________________
  # Facet by adiv metric when there's more than one.
  #________________________________________________________
  if (length(params[['adiv']]) > 1)
    params[['facet.by']] %<>% c(".adiv")
  
  
  #________________________________________________________
  # Tell corrplot_build to use lm(y ~ log(x))
  #________________________________________________________
  params[['model']] <- "log"
  
  
  #________________________________________________________
  # Convert user's `layers` spec to layer names.
  #________________________________________________________
  layer_names <- local({
    layerlist <- c('t' = "trend", 's' = "scatter")
    
    layer_match(params[['layers']], choices = layerlist, default = "t") %>%
      c('ggplot', ., 'labs', 'theme_bw')
  })
  
  if (!'scatter' %in% layer_names) params[['shape.by']] <- NULL
  
  
  
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
  
  initLayer(layer_names)
  
  if (!is_null(params[['color.by']])) initLayer("color")
  if (!is_null(params[['shape.by']])) initLayer("shape")
  
  
  #________________________________________________________
  # Add default layer parameters.
  #________________________________________________________
  setLayer("ggplot", mapping = list(x = ".depth", y = attr(data, 'response', exact = TRUE)))
  setLayer("xaxis",  labels = si_units)
  setLayer("labs",   x = "Rarefaction Depth", y = params[['adiv']])
  
  if (length(params[['adiv']]) > 1) {
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






