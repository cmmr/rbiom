
# See also ggpmisc package
# https://docs.r4photobiology.info/ggpmisc/articles/model-based-annotations.html

#' Visualize alpha diversity with scatterplots and trendlines.
#' 
#' @inherit adiv_boxplot params
#' @family visualization
#' 
#' @param biom   A BIOM object, as returned from [read_biom()].
#' 
#' @param x   A numeric metadata column name to use for the x-axis. Required.
#'           
#' @param layers   \code{"trend"}, \code{"scatter"}. Single letter 
#'        abbreviations are also accepted. For instance, 
#'        \code{c("trend", "scatter")} is equivalent to \code{c("t", "s")} and 
#'        \code{"ts"}. See \code{vignette("corrplots")} for examples of each. 
#'        Default: \code{"t"}.
#'        
#' @param color.by,facet.by,limit.by   Metadata columns to use for aesthetics 
#'        and partitioning. See below for details. Default: \code{NULL}
#'        
#' @param model   What type of trendline to fit to the data. Options are: 
#'        \itemize{
#'          \item{\code{"lm"} - }{  Linear model: \code{stats::lm(formula = y ~ x)}.) }
#'          \item{\code{"log"} - }{ Logarithmic model: \code{stats::lm(formula = y ~ log(x))}. }
#'          \item{\code{"gam"} - }{ Generalized additive model: \code{mgcv::gam(formula = y ~ s(x, bs = "cs"), method = "REML")}. }
#'        }
#'        Default: \code{"lm"} \cr\cr
#'        You can alternatively provide a list of length two where the first 
#'        element is a character vector of length 1 naming a function, and the 
#'        second element is a list of arguments to pass to that function. One 
#'        of the function's arguments must be named 'formula'. 
#'        For example, \code{model = list("stats::lm", list(formula = y ~ x))}.
#'        
#' @param stats   Which statistic to display on the plot. Options are: 
#'        \itemize{
#'          \item{\code{"fit"} - }{ How well does the model fit the data? }
#'          \item{\code{"terms"} - }{ How strongly does 'x' influence 'y'? }
#'          \item{\code{"emmeans"} - }{ Is the average 'y' value non-zero? }
#'          \item{\code{"emtrends"} - }{ Does any trendline have a non-zero slope? }
#'          \item{\code{"emm_pairs"} - }{ Are the means of any trendlines different? }
#'          \item{\code{"emt_pairs"} - }{ Are the slopes of any trendlines different? }
#'          \item{\code{"hide"} - }{ Don't show stats on the plot, but still compute them. }
#'          \item{\code{"none"} - }{ Do not compute or show statistics. }
#'        }
#'        Default: \code{"emtrends"} \cr\cr
#'        Note: \code{"emm_pairs"} and \code{"emt_pairs"} can only be calculated
#'        when using a \code{color.by} metadata column with more than one level. \cr\cr
#'        Statistical tests are run separately on each facet. P-values are 
#'        adjusted for multiple comparisons by considering all facets together. 
#'        Unless \code{stats = "none"}, all stats are attached to the plot as 
#'        \code{attr(,'stats')}.
#'
#' @param p.adj   Method to use for multiple comparisons adjustment of p-values.
#'        Run \code{p.adjust.methods} for a list of available options.
#'        Default: \code{"fdr"}.
#' 
#' @param ci   The confidence interval to display around the trendline. 
#'        Default: \code{95}.
#'        
#' @param caption   Display information about the method used for trendline
#'        fitting beneath the plot. Default: \code{FALSE}.
#'        
#' @param ...   Additional parameters to pass along to ggplot2
#'        functions. Prefix a parameter name with either \code{t.} or 
#'        \code{s.}/\code{pt.} to ensure it gets passed to (and only to) 
#'        \link[ggplot2]{geom_smooth} or \link[ggplot2]{geom_point}, 
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
#'     biom <- sample_rarefy(hmp50)
#'     adiv_corrplot(biom, "Age", color.by="Body Site", adiv=c("shannon", "otus"), facet.by = "Sex", ci = 90) 
#'     
adiv_corrplot <- function (
    biom, x, adiv = "Shannon", layers = "t", 
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    model = "lm", stats = "emtrends", 
    p.adj = "fdr", ci = 95, caption = FALSE, ...) {
  
  
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
    adiv %<>% validate_arg(biom, 'adiv', 'adiv', n = c(1,Inf))
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
    biom = params[['biom']],
    adiv = params[['adiv']],
    long = TRUE,
    md   = unique(c(
      params[['x']], 
      names(params[['color.by']]), 
      params[['facet.by']] )))
  
  
  #________________________________________________________
  # Facet by adiv metric when there's more than one.
  #________________________________________________________
  ylab <- NULL
  ycol <- attr(data, 'response', exact = TRUE)
  hist <- attr(data, 'history',  exact = TRUE)
  
  if (length(unique(data[['.adiv']])) > 1) {
    ylab <- "Diversity (multiple metrics)"
    data <- drop_cols(data, '.depth')
    params[['facet.by']] <- unique(c(params[['facet.by']], ".adiv"))
    
  } else {
    data <- rename_response(data, paste0(".", unique(data[['.adiv']])))
    ycol <- attr(data, 'response', exact = TRUE)
    ylab <- sprintf("Diversity (%s)", unique(data[['.adiv']]))
    data <- drop_cols(data, '.depth', '.adiv')
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
    layerlist <- c(
      't' = "trend", 'c' = "confidence", 
      's' = "scatter", 'n' = "name", 'r' = "residual")
    
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
  
  initLayer(setdiff(layer_names, c('trend', 'confidence')))
  
  if (!is_null(params[['color.by']])) initLayer("color")
  if (!is_null(params[['shape.by']])) initLayer("shape")
  
  
  #________________________________________________________
  # Merge trend and confidence into a single layer.
  #________________________________________________________
  if (any(c('trend', 'confidence') %in% layer_names)) {
    initLayer("trend")
    if (!'trend'      %in% layer_names) setLayer("trend", color = NA)
    if (!'confidence' %in% layer_names) setLayer("trend", se = FALSE)
  }
  
  
  #________________________________________________________
  # Add default layer parameters.
  #________________________________________________________
  setLayer("ggplot", mapping = list(x = params[['x']], y = attr(data, 'response', exact = TRUE)))
  setLayer("labs", x = params[['x']], y = ylab)
  if (length(params[['adiv']]) > 1) setLayer("facet", scales = "free_y")
  
  
  
  
  #________________________________________________________
  # Convert layer definitions into a plot.
  #________________________________________________________
  p <- corrplot_build(layers)
  
  attr(p, 'history') <- history
  
  
  set_cache_value(cache_file, p)
  return (p)
}






