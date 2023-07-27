#' Visualize taxa abundance with scatterplots and trendlines.
#' 
#' @name taxa_corrplot
#' 
#' @inherit adiv_corrplot params sections return
#' @inherit taxa_boxplot params
#' 
#' @family plotting
#' 
#'        
#' @param rank   What rank of taxa to display. E.g. "Phylum", "Genus", etc. 
#'        Run \code{taxa_ranks()} to see all options for a given BIOM object. 
#'        The default, \code{NULL}, selects the lowest level.
#'        
#' @param taxa   Which taxa to display. An integer value will show the top n
#'        most abundant taxa. A value 0 <= n < 1 will show any taxa with that 
#'        mean abundance or greater (e.g. 0.1). A character vector of
#'        taxon names will show only those taxa. Default: \code{6}.
#'        
#' @param stats   Which statistic to display on the plot. Options are: 
#'        \describe{
#'          \item{\code{"fit"}}{ How well does the model fit all the data? }
#'          \item{\code{"emtrends"}}{ Does any trendline have a non-zero slope? }
#'          \item{\code{"emt_pairs"}}{ Are the slopes of any trendlines different? }
#'          \item{\code{"emm_pairs"}}{ Are the means of any trendlines different? }
#'          \item{\code{"hide"}}{ Don't show stats on the plot, but still compute them. }
#'          \item{\code{"none"}}{ Do not compute or show statistics. }
#'        }
#'        Statistical tests are run separately on each facet (except 
#'        \code{"fit"}). P-values are adjusted for multiple comparisons by 
#'        considering all facets together. Unless \code{stats = "none"}, all 
#'        stats are attached to the plot as \code{attr(,'stats')}. 
#'        Default: \code{"emtrends"}
#'        
#' @param ...   Additional parameters to pass along to ggplot2
#'        functions. Prefix a parameter name with either \code{t.} or 
#'        \code{s.}/\code{pt.} to ensure it gets passed to (and only to) 
#'        \link[ggplot2]{geom_smooth} or \link[ggplot2]{geom_point}, 
#'        respectively. For instance, \code{s.size = 2} ensures only the 
#'        scatterplot points have their size set to \code{2}.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     taxa_corrplot(rarefy(hmp50), "BMI", color.by="Body Site") 
#'     

taxa_corrplot <- function (
    biom, x, rank = NULL, taxa = 6, layers = "t",
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    model = "linear", stats = "emtrends", 
    p.top = Inf, p.adj = "fdr", ci = 95, ...) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(c(as.list(environment()), list(...)), eval)
  cache_file <- get_cache_file("taxa_corrplot", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  arg_str <- as.args(params, fun = taxa_corrplot, indent = 2)
  history <- paste0(collapse = "\n", c(
    attr(biom, 'history', exact = TRUE),
    sprintf("fig  <- taxa_corrplot(%s)", arg_str) ))
  remove(list = setdiff(ls(), c("params", "history", "cache_file")))
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  params %<>% within({
    if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
    if (is_scalar_logical(stats)) stats %<>% ifelse("nz", "none")
    stopifnot(is_string(stats, c("fit", "emtrends", "emt_pairs", "emm_pairs", "hide", "none")))
    rank %<>% validate_arg(biom, 'rank', n = 1, default = tail(c('OTU', taxa_ranks(biom)), 1))
  })
  
  
  #________________________________________________________
  # Subset biom by requested metadata and aes.
  #________________________________________________________
  params %<>% metadata_params(contraints = list(
    x        = list(n = 1, col_type = "num"),
    color.by = list(n = 0:1),
    facet.by = list(col_type = "cat"),
    limit.by = list() ))
  
  
  #________________________________________________________
  # Compute taxa abundance values.
  #________________________________________________________
  
  if (is_rarefied(params[['biom']]))
    params[['biom']] %<>% as_percent()
  
  data <- taxa_table(
    biom = params[['biom']],
    rank = params[['rank']],
    taxa = params[['taxa']],
    md   = unique(c(
      params[['x']], 
      names(params[['color.by']]), 
      params[['facet.by']] )))
  
  params[['facet.by']] %<>% c(".taxa")
  
  
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
    'function' = taxa_corrplot,
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
  
  setLayer("labs", x = params[['x']], y = local({
    ylab <- paste(params[['rank']], "Abundance")
    if (!is_rarefied(params[['biom']])) ylab %<>% paste("[UNRAREFIED]")
    return (ylab) }))
  
  setLayer("facet", scales = "free_y")
  
  
  
  
  # Convert layer definitions into a plot.
  #________________________________________________________
  p <- corrplot_build(layers)
  
  attr(p, 'history') <- history
  
  
  set_cache_value(cache_file, p)
  return (p)
}






