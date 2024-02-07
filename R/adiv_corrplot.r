
# See also ggpmisc package
# https://docs.r4photobiology.info/ggpmisc/articles/model-based-annotations.html

#' Visualize alpha diversity with scatterplots and trendlines.
#' 
#' @inherit documentation_corrplot
#' @inherit documentation_default
#' 
#' @family alpha_diversity
#' @family visualization
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     adiv_corrplot(biom, "age", color.by = list(Body = c("a", "sa")))
#'     adiv_corrplot(
#'       biom     = biom, 
#'       x        = "Age", 
#'       adiv     = c("OTUs", "Shannon", "Simpson"), 
#'       color.by = "Body Site", 
#'       facet.by = "Sex", 
#'       layers   = "trend" )
#'     
adiv_corrplot <- function (
    biom, x, adiv = "Shannon", layers = "tc", 
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    formula = y ~ x, engine = "lm", level = 0.95, 
    trans = "none", caption = TRUE, ...) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment(), ...)
  cmd    <- sprintf("adiv_corrplot(%s)", as.args(params, 2, adiv_corrplot))
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  with(params, {
    
    validate_adiv(max = Inf)
    validate_meta_aes('x', col_type = "num")
    validate_meta_aes('color.by', col_type = "cat", null_ok = TRUE)
    validate_meta_aes('facet.by', col_type = "cat", null_ok = TRUE, max = Inf)
    validate_meta_aes('limit.by',                   null_ok = TRUE, max = Inf)
    
    sync_metadata()
  })
  
  
  #________________________________________________________
  # Compute alpha diversity values.
  #________________________________________________________
  with(params, {
    
    
    # Compute each adiv metric separately: shannon, etc
    #________________________________________________________
    .ggdata <- adiv_table(
      biom  = biom,
      adiv  = adiv,
      md    = unique(c(x, names(color.by), facet.by)),
      trans = trans )
    
    
    # axis titles
    #________________________________________________________
    .xlab <- x
    if (length(adiv) == 1)
      .ylab <- switch(
        EXPR    = adiv,
        'OTUs'  = "Observed OTUs",
        'Depth' = "Sequencing Depth",
        paste(adiv, "Diversity") )
    
    
    
    # Facet on multiple adiv metrics
    #________________________________________________________
    if (length(adiv) > 1) {
      .free_y <- TRUE
      facet.by %<>% c('.adiv')
      .ylab <- aa("\u03B1 Diversity", display = '"\\u03B1 Diversity"')
    }
    
  })
  
  
  
  #________________________________________________________
  # Create and customize layer definitions.
  #________________________________________________________
  init_corrplot_layers(params)
  
  
  #________________________________________________________
  # Build the plot.
  #________________________________________________________
  fig <- params %>% 
    plot_facets() %>% 
    plot_build()
  
  
  attr(fig, 'cmd') <- cmd
  set_cache_value(cache_file, fig)
  
  return (fig)
}





