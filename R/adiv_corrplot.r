
# See also ggpmisc package
# https://docs.r4photobiology.info/ggpmisc/articles/model-based-annotations.html

#' Visualize alpha diversity with scatterplots and trendlines.
#' 
#' @inherit documentation_test.trends
#' @inherit documentation_model.lm
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
#'     adiv_corrplot(biom, "age", color.by="body", adiv=c("sha", "otu"), facet.by="sex")
#'     
adiv_corrplot <- function (
    biom, x, adiv = "Shannon", layers = "t", 
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    test = "trends", model = "lm", 
    p.adj = "fdr", level = 0.95, caption = TRUE, ...) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment(), ...)
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
        biom = biom,
        adiv = adiv,
        md   = unique(c(x, names(color.by), facet.by)) )
    
    
    # Facet on multiple adiv metrics
    #________________________________________________________
    if (length(adiv) > 1) {
      .free_y <- TRUE
      facet.by %<>% c('.adiv')
    }
    
  })
  
  
  
  #________________________________________________________
  # Create and customize layer definitions.
  #________________________________________________________
  init_corrplot_layers(params)
  
  
  #________________________________________________________
  # y-axis title and scale
  #________________________________________________________
  set_layer(params, 'labs', y = with(params, {
    
    if (length(adiv) > 1) {
      
      "Diversity (\u03B1)" %>% 
        aa(display = '"Diversity (\\u03B1)"')
      
    } else {
      
      switch (
        EXPR    = adiv,
        'OTUs'  = "Observed OTUs",
        'Depth' = "Sequencing Depth",
        paste(adiv, "Diversity") )
    }
    
  }))
  
  
  
  #________________________________________________________
  # Build the plot.
  #________________________________________________________
  fig <- params %>% 
    plot_facets() %>% 
    corrplot_stats() %>%
    plot_build()
  
  
  set_cache_value(cache_file, fig)
  
  return (fig)
}





