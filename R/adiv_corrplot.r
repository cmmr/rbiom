
# See also ggpmisc package
# https://docs.r4photobiology.info/ggpmisc/articles/model-based-annotations.html

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
#' @param metric   Alpha diversity metric(s) to use. Options are: 
#'        \code{"OTUs"}, \code{"Shannon"}, \code{"Chao1"}, \code{"Simpson"}, 
#'        and/or \code{"InvSimpson"}. Default: \code{"OTUs"}.
#'           
#' @param layers   \code{"trend"}, \code{"scatter"}. Single letter 
#'        abbreviations are also accepted. For instance, 
#'        \code{c("trend", "scatter")} is equivalent to \code{c("t", "s")} and 
#'        \code{"ts"}. See \code{vignette("corrplots")} for examples of each. 
#'        Default: \code{"t"}.
#'        
#' @param model   What type of trendline to fit to the data. Options are: 
#'        \code{"linear"}, \code{"logarithmic"}, or \code{"local"}. You can
#'        alternatively provide a list of length two where the first element is
#'        a character vector of length 1 naming a function, and the second 
#'        element is a list of arguments to pass to that function. One of the 
#'        function's arguments must be named 'formula'. For example, 
#'        \code{model = list("stats::lm", list(formula = y ~ x))}.
#'        Default: \code{"linear"}.
#' 
#' @param ci   The confidence interval to display around the fitted curve. Set
#'        to \code{FALSE} to hide the confidence interval. Default: \code{95}.
#'        
#' @param color.by,facet.by,limit.by   Metadata columns to use for aesthetics 
#'        and partitioning. See below for details. Default: \code{NULL}
#'        
#' @param ...   Additional parameters to pass along to ggplot2
#'        functions. Prefix a parameter name with either \code{t.} or \code{s.}
#'        to ensure it gets passed to (and only to) 
#'        \link[ggplot2]{geom_point} or \link[ggplot2]{geom_smooth}, 
#'        respectively. For instance, \code{s.size = 2} ensures only the 
#'        scatterplot points have their size set to \code{2}.
#'        
#' @return A \code{ggplot2} plot. The computed data points and statistics will 
#'         be attached as \code{attr(p, 'data')} and \code{attr(p, 'stats')}, 
#'         respectively.
#' 
#' 
#' @section Aesthetics and Partitions:
#' 
#' Metadata can be used to flexibly subset, partition, and apply aesthetics 
#' when creating a plot. Common use cases are provided below. More thorough 
#' documentation is available at \url{https://cmmr.github.io/rbiom}.
#' 
#' \preformatted{  ## Colors ----------------------------
#'   color.by = "Body Site"
#'   color.by = list('Body Site' = "bright")
#'   color.by = list('Body Site' = c("Stool", "Saliva"))
#'   color.by = list('Body Site' = list('values' = c("Stool", "Saliva"), 'colors' = "bright"))
#'   color.by = list('Body Site' = c('Stool' = "blue", 'Saliva' = "green"))
#'   
#'   ## Facets ----------------------------
#'   facet.by = "Body Site"
#'   facet.by = c("Body Site", "Sex")
#'   facet.by = list('Body Site' = c("Stool", "Saliva"), "Sex")
#'   
#'   ## Limits ----------------------------
#'   limit.by = list('Sex' = "Male", 'Age' = c(20,40))
#'   limit.by = list('Body Site' = c("Saliva", "Anterior nares"), 'Age' = c(NA,35))
#' }
#' 
#' \itemize{
#'   \item{\code{color.by} - }{Any metadata column. (Max 1)}
#'   \item{\code{facet.by} - }{Only categorical metadata column(s).}
#'   \item{\code{limit.by} - }{Any metadata column(s).}
#' }
#' 
#' All built-in color palettes are colorblind-friendly.
#' 
#' The available categorical palette names are: \code{"okabe"}, \code{"carto"}, 
#' \code{"r4"}, \code{"polychrome"}, \code{"tol"}, \code{"bright"}, 
#' \code{"light"}, \code{"muted"}, \code{"vibrant"}, \code{"tableau"}, 
#' \code{"classic"}, \code{"alphabet"}, \code{"tableau20"}, \code{"kelly"}, 
#' and \code{"fishy"}.
#' 
#' The available numeric palette names are: \code{"reds"}, \code{"oranges"}, 
#' \code{"greens"}, \code{"purples"}, \code{"grays"}, \code{"acton"}, 
#' \code{"bamako"}, \code{"batlow"}, \code{"bilbao"}, \code{"buda"}, 
#' \code{"davos"}, \code{"devon"}, \code{"grayC"}, \code{"hawaii"}, 
#' \code{"imola"}, \code{"lajolla"}, \code{"lapaz"}, \code{"nuuk"}, 
#' \code{"oslo"}, \code{"tokyo"}, \code{"turku"}, \code{"bam"}, 
#' \code{"berlin"}, \code{"broc"}, \code{"cork"}, \code{"lisbon"}, 
#' \code{"roma"}, \code{"tofino"}, \code{"vanimo"}, \code{"vik"}
#'  
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     adiv_corrplot(rarefy(hmp50), "Age", color.by="Body Site", metric=c("shannon", "otus"), facet.by = "Sex", ci = 90) 
#'     
adiv_corrplot <- function (
    biom, x, metric = "OTUs", layers = "t", 
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    model = "linear", p.adj = "fdr", ci = 95, ...) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(c(as.list(environment()), list(...)), eval)
  cache_file <- get_cache_file("adiv_corrplot", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  arg_str <- as.args(params, fun = adiv_corrplot, indent = 2)
  history <- paste0(collapse = "\n", c(
    attr(biom, 'history', exact = TRUE),
    sprintf("fig  <- adiv_corrplot(%s)", arg_str) ))
  remove(list = setdiff(ls(), c("params", "history", "cache_file")))
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  params %<>% within({
    if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
    metric %<>% validate_arg(biom, 'metric', 'adiv', n = c(1,Inf))
  })
  
  
  #________________________________________________________
  # Subset biom by requested metadata and aes.
  #________________________________________________________
  params %<>% metadata_params(contraints = list(
    x        = list(n = 1, col_type = "num"),
    color.by = list(n = c(0, 1)),
    facet.by = list(n = c(0, Inf), col_type = "cat"),
    limit.by = list(n = c(0, Inf)) ))
  
  
  #________________________________________________________
  # Compute alpha diversity values.
  #________________________________________________________
  data <- adiv_table(
    biom    = params[['biom']],
    metrics = params[['metric']],
    long    = TRUE,
    md      = unique(c(
      params[['x']], 
      names(params[['color.by']]), 
      params[['facet.by']] )))
  
  
  #________________________________________________________
  # Facet by metric when there's more than one.
  #________________________________________________________
  ylab <- NULL
  ycol <- attr(data, 'response', exact = TRUE)
  hist <- attr(data, 'history',  exact = TRUE)
  
  if (length(unique(data[['.metric']])) > 1) {
    ylab <- "Diversity (multiple metrics)"
    data <- drop_cols(data, '.depth')
    params[['facet.by']] <- unique(c(params[['facet.by']], ".metric"))
    
  } else {
    data <- rename_response(data, paste0(".", unique(data[['.metric']])))
    ycol <- attr(data, 'response', exact = TRUE)
    ylab <- sprintf("Diversity (%s)", unique(data[['.metric']]))
    data <- drop_cols(data, '.depth', '.metric')
  }
  
  if (!is_rarefied(params[['biom']]))
    ylab %<>% paste("\nWARNING: DATA NOT RAREFIED")
  
  attr(data, 'response') <- ycol
  attr(data, 'history')  <- hist
  remove("ycol", "hist")
  
  
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
    'function' = adiv_corrplot,
    'xcol'     = params[['x']],
    'ycol'     = attr(data, 'response', exact = TRUE),
    'xmode'    = "numeric" )
  
  initLayer(layer_names)
  
  if (!is_null(params[['color.by']])) initLayer("color")
  if (!is_null(params[['shape.by']])) initLayer("shape")
  
  
  
  #________________________________________________________
  # Add default layer parameters.
  #________________________________________________________
  setLayer("ggplot", mapping = list(x = params[['x']], y = attr(data, 'response', exact = TRUE)))
  setLayer("labs", x = params[['x']], y = ylab)
  if (length(params[['metric']]) > 1) setLayer("facet", scales = "free_y")
  
  
  
  
  #________________________________________________________
  # Convert layer definitions into a plot.
  #________________________________________________________
  p <- corrplot_build(layers)
  
  attr(p, 'history') <- history
  
  
  set_cache_value(cache_file, p)
  return (p)
}






