
#' Visualize beta diversity with scatterplots and trendlines.
#' 
#' @inherit documentation_test.trends
#' @inherit documentation_model.lm
#' @inherit documentation_corrplot
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' @family visualization
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     bdiv_corrplot(biom, "Age", color.by="==Body Site", facet.by="==Sex")
#'     
bdiv_corrplot <- function (
    biom, x, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, layers = "t", 
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    within = NULL, between = NULL, 
    test = "trends", model = "lm", trans = "rank", 
    p.adj = "fdr", level = 0.95, caption = TRUE, ...) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment(), ...)
  cmd    <- sprintf("bdiv_corrplot(%s)", as.args(params, 2, bdiv_corrplot))
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
    
    validate_bdiv(max = Inf)
    validate_meta_aes('x', col_type = "num")
    validate_meta_aes('color.by', col_type = "cat", cmp = TRUE, null_ok = TRUE)
    validate_meta_aes('facet.by', col_type = "cat", cmp = TRUE, null_ok = TRUE, max = Inf)
    validate_meta_aes('limit.by',                   cmp = TRUE, null_ok = TRUE, max = Inf)
    
    # Validates and appends to `within` and `between`.
    validate_meta_cmp(c('color.by', 'facet.by', 'limit.by'))
    
    sync_metadata()
  })
  
  
  #________________________________________________________
  # Compute beta diversity values.
  #________________________________________________________
  with(params, {
    
    
    # Compute each bdiv metric separately: jaccard, etc
    #________________________________________________________
    .ggdata <- bdiv_table(
      biom     = biom, 
      bdiv     = bdiv, 
      weighted = weighted, 
      tree     = tree, 
      md       = unique(c(x, names(color.by), facet.by)), 
      within   = within, 
      between  = between,
      trans    = trans )
    
    
    # axis titles
    #________________________________________________________
    .xcol <- x
    .ycol <- '.distance'
    
    .xlab <- aa(paste("\u0394", x), display = paste0('"\\u0394 ', x, '"'))
    .ylab <- paste(bdiv, "Distance")
    
    
    # Facet on multiple bdiv metrics
    #________________________________________________________
    if (length(bdiv) > 1) {
      facet.by %<>% c('.bdiv')
      .ylab <- aa("\u03B2 Dissimilarity", display = '"\\u03B2 Dissimilarity"')
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
    corrplot_stats() %>%
    plot_build()
  
  
  attr(fig, 'cmd') <- cmd
  set_cache_value(cache_file, fig)
  
  return (fig)
}





