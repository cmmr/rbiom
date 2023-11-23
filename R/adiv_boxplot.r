#' Visualize alpha diversity with boxplots.
#' 
#' @inherit documentation_boxplot
#' @inherit documentation_default
#' 
#' @family alpha_diversity
#' @family visualization
#' 
#' @export
#' 
#' @examples
#'     library(rbiom)
#'     
#'     biom <- sample_rarefy(hmp50)
#'     
#'     adiv_boxplot(biom, x="Sex", adiv=c("otu", "shan"), color.by="Body Site")
#'     
#'     adiv_boxplot(biom, x="body", adiv="all", layers="p", color.by="sex", flip=TRUE)
#'     
#'     adiv_boxplot(biom, x="Body Site", color.by="Body Site")
#'     
#'     
#'     # Each plot object includes additional information.
#'     fig <- adiv_boxplot(biom, x="Body Site")
#'     
#'     ## Provenance History ---------------------
#'     fig$history
#'     
#'     ## Computed Data Points -------------------
#'     fig$data
#'     
#'     ## Statistics Table -----------------------
#'     fig$stats
#'     
#'     ## ggplot2 Command ------------------------
#'     fig$code
#'     

adiv_boxplot <- function (
    biom, x = NULL, adiv = "Shannon", layers = "bld",
    color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, limit.by = NULL, 
    flip = FALSE, stripe = NULL, p.adj = "fdr", p.label = 0.05, ci = "ci", level = 0.95, 
    outliers = NULL, xlab.angle = 'auto', ...) {
  
  
  validate_biom(clone = FALSE)
  
  params  <- eval_envir(environment(), ...)
  history <- append_history('fig ', params)
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
    
    validate_meta_aes('x',          col_type = "cat", null_ok = TRUE)
    validate_meta_aes('color.by',   col_type = "cat", null_ok = TRUE)
    validate_meta_aes('pattern.by', col_type = "cat", null_ok = TRUE)
    validate_meta_aes('shape.by',   col_type = "cat", null_ok = TRUE)
    validate_meta_aes('facet.by',   col_type = "cat", null_ok = TRUE, max = Inf)
    validate_meta_aes('limit.by',                     null_ok = TRUE, max = Inf)
    
    sync_metadata()
  })
  
  
  
  #________________________________________________________
  # Compute alpha diversity and set up adiv facet.
  #________________________________________________________
  with(params, {
    
    # Compute each adiv metric separately: shannon, etc
    #________________________________________________________
    .ggdata <- adiv_table(
      biom = biom, 
      adiv = adiv, 
      md   = TRUE )
    
    
    # Facet on multiple adiv metrics
    #________________________________________________________
    if (length(adiv) > 1) {
      facet.by %<>% c('.adiv')
      .free_y <- TRUE
    }
    
  })
  
  
  
  #________________________________________________________
  # Initialize the layers environment.
  #________________________________________________________
  init_boxplot_layers(params)
  
  
  
  #________________________________________________________
  # Default y-axis title
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
      
    }}))
  
  
  
  
  #________________________________________________________
  # Use the generalized boxplot functions to make the plot
  #________________________________________________________
  fig <- params %>%
    boxplot_facets() %>%
    boxplot_stats() %>%
    plot_build()
  
  
  
  attr(fig, 'history') <- history
  set_cache_value(cache_file, fig)
  
  return (fig)
}


