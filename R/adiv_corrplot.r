#' Visualize alpha diversity with scatterplots and trendlines.
#' 
#' @name adiv_corrplot
#' 
#' @family plotting
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#' 
#' @param x   A numeric metadata column name to use for the x-axis. Required.
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
#'        Default: \code{method = "lm", formula = y ~ x}.
#' 
#' @param ci   The confidence interval to display around the fitted curve. Set
#'        to \code{FALSE} to hide the confidence interval. Default: \code{95}.
#'        
#' @param ...   Additional parameters to pass along to ggplot2
#'        functions. Prefix a parameter name with either \code{p.} or \code{s.}
#'        to ensure it gets passed to (and only to) 
#'        \link[ggplot2]{geom_point} or 
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
#'     adiv_corrplot(rarefy(hmp50), "Age", color.by="Body Site", metric=c("shannon", "otus"), facet.by = "Sex", ci = 90)
#'     

adiv_corrplot <- function (
    biom, x, layers = "ps", metric = "OTUs",
    color.by = NULL, facet.by = NULL, colors = NULL, facets = NULL,
    method = "lm", formula = y ~ x, ci = 95, ...) {
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  metric <- as.vector(validate_metrics(biom, metric, "adiv", multi=TRUE))
  x <- validate_metrics(biom, x, "meta")
  stopifnot(identical(attr(x, 'mode', exact = TRUE), "numeric"))
  x <- as.vector(x)
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  params <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  history <- sprintf("adiv_corrplot(%s)", as.args(params, fun = adiv_corrplot))
  
  
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
  # Compute alpha diversity values.
  #________________________________________________________
  data <- adiv_table(
    biom    = biom,
    metrics = metric,
    long    = TRUE,
    md      = unique(c(x, color.by, facet.by)), 
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
    'function' = adiv_corrplot,
    'xcol'     = x,
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
  setLayer("ggplot", mapping = list(x = x, y = ".value"))
  setLayer("labs", x = x, y = local({
    ylab <- ifelse(length(metric) == 1, metric, "Diversity")
    if (!is_rarefied(biom)) ylab %<>% paste("[UNRAREFIED]")
    return (ylab) }))
  if (length(metric) > 1) setLayer("facet", scales = "free_y")
  
  
  
  
  # Convert layer definitions into a plot.
  #________________________________________________________
  p <- corrplot_build(layers)
  
  
  # Attach history of biom modifications and this call.
  #________________________________________________________
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  return (p)
  
  
}






