#' Visualize rarefaction curves with scatterplots and trendlines.
#' 
#' @name rare_corrplot
#' 
#' @family plotting
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#'           
#' @param layers   \code{"point"}, and/or \code{"smooth"}. Single letter 
#'        abbreviations are also accepted. For instance, 
#'        \code{c("point", "smooth")} is equivalent to \code{c("p", "s")} and 
#'        \code{"ps"}. Default: \code{"ps"}.
#' 
#' @param metric   Alpha diversity metric(s) to use. Options are: 
#'        \code{"OTUs"}, \code{"Shannon"}, \code{"Chao1"}, \code{"Simpson"}, 
#'        and/or \code{"InvSimpson"}. Default: \code{"OTUs"}.
#' 
#' @param depths   Rarefaction depths to show in the plot. Passed to
#'        \link{adiv_table}. The default, \code{"multi_even"}, uses a heuristic
#'        to pick 10 evenly spaced depths.
#'        
#' @param color.by,facet.by   Metadata column to color and/or facet by. If that
#'        column is a \code{factor}, the ordering of levels will be maintained 
#'        in the plot. Default: \code{color.by = NULL, facet.by = NULL}
#'        
#' @param colors,facets   Names of the colors and/or facets values to use in
#'        the plot. Available values for colors are given by \code{colors()}. 
#'        Use a named character vector to map them to specific factor levels in
#'        the metadata. \code{facets} are coerced to unnamed character vectors.
#'        If the length of these vectors is less than the values present in 
#'        their corresponding metadata column, then the data set will be 
#'        subseted accordingly. Default: \code{colors = NULL, facets = NULL}
#'        
#' @param method,formula   Passed to \link[ggplot2]{stat_smooth} to
#'        define the fitted curve.
#'        Default: \code{method = "lm", formula = y ~ log(x)}.
#' 
#' @param ci   The confidence interval to display around the fitted curve. Set
#'        to \code{FALSE} to hide the confidence interval. Default: \code{95}.
#'        
#' @param vline   Where to draw a vertical line on the plot, intended to show
#'        a particular rarefaction depth. Default: \code{NULL}.
#'        
#' @param ...   Additional parameters to pass along to ggplot2 
#'        functions. Prefix a parameter name with either \code{p.}, \code{s.}, 
#'        or \code{v.} to ensure it gets passed to (and only to) 
#'        \link[ggplot2]{geom_point}, \link[ggplot2]{vline}, or 
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
    biom, layers = "ps", metric = "OTUs", depths = NULL,
    color.by = NULL, facet.by = NULL, colors = NULL, facets = NULL,
    method = "lm", formula = y ~ log(x), ci = 95, vline = NULL, ...) {
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  metric <- as.vector(validate_metrics(biom, metric, "adiv", multi=TRUE))
  stopifnot(is.null(depths) || is.numeric(depths))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  params <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  history <- sprintf("rare_corrplot(%s)", as.args(params, fun = rare_corrplot))
  
  
  #________________________________________________________
  # Subset the biom file.
  #________________________________________________________
  for (i in c("color", "facet")) {
    key <- get(paste0(i, ".by"))
    val <- get(paste0(i, "s"))
    if (!is.null(key) && !is.null(val))
      biom <- subset(
        biom = biom,
        expr = a %in% b,
        env  = list(a = as.name(key), b = val),
        fast = TRUE )
  }
  
  
  #________________________________________________________
  # Pull rarefaction stats for each depth.
  #________________________________________________________
  data <- adiv_table(
    biom    = biom,
    rarefy  = ifelse(is.null(depths), "multi_even", depths),
    metrics = metric,
    long    = TRUE,
    md      = unique(c(color.by, facet.by)), 
    safe    = TRUE )
  
  
  #________________________________________________________
  # Facet by metric when there's more than one.
  #________________________________________________________
  if (length(metric) > 1)
    params[['facet.by']] <- unique(c(facet.by, ".metric"))
  
  
  #________________________________________________________
  # Initialize the `layers` object.
  #________________________________________________________
  layers <- structure(
    list(),
    'data'     = data,
    'params'   = params,
    'function' = rare_corrplot,
    'xcol'     = ".depth",
    'ycol'     = ".value",
    'xmode'    = "numeric" )
  
  layer_names <- params[['layers']] %>%
    layer_match(c('p' = "point", 's' = "smooth"), "ps") %>%
    c('ggplot', ., 'labs', 'theme_bw')
  initLayer(layer_names)
  remove("layer_names")
  
  
  #________________________________________________________
  # Add default layer parameters.
  #________________________________________________________
  setLayer("ggplot", mapping = list(x = ".depth", y = ".value"))
  setLayer("labs",   x = "Rarefaction Depth")
  setLayer("labs",   y = ifelse(length(metric) == 1, metric, "Diversity"))
  if (length(metric) > 1) setLayer("facet", scales = "free_y")
  
  if (!is.null(vline))
    setLayer("vline", xintercept = vline, color = "red", linetype="dashed")
  
  
  
  # Convert layer definitions into a plot.
  #________________________________________________________
  p <- corrplot_build(layers)
  
  
  # Attach history of biom modifications and this call.
  #________________________________________________________
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  return (p)
  
  
}






