#' Visualize taxa abundance with scatterplots and trendlines.
#' 
#' @inherit documentation_test.trends
#' @inherit documentation_model.lm
#' @inherit documentation_corrplot
#' @inherit documentation_default
#' 
#' @family taxa_abundance
#' @family visualization
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- sample_rarefy(hmp50)
#'     taxa_corrplot(biom, "BMI", color.by="Body Site", taxa = 4) 
#'     

taxa_corrplot <- function (
    biom, x, rank = -1, taxa = 6, layers = "t",
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    test = "trends", model = "lm", level = 0.95, 
    p.top = Inf, p.adj = "fdr", caption = FALSE, ...) {
  
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
    validate_rank(max = Inf)
    validate_meta_aes('x', col_type = "num")
    validate_meta_aes('color.by', null_ok = TRUE)
    validate_meta_aes('facet.by', null_ok = TRUE, max = Inf, col_type = "cat")
    validate_meta_aes('limit.by', null_ok = TRUE, max = Inf)
    
    sync_metadata()
  })
  
  
  
  #________________________________________________________
  # Compute taxa abundance values.
  #________________________________________________________
  with(params, {
    
    if (is_rarefied(biom))
      biom %<>% as_percent()
    
    
    # Compute each rank separately: phylum, genus, etc
    #________________________________________________________
    .ggdata <- taxa_table(
      biom = biom,
      rank = rank,
      taxa = taxa,
      md   = unique(c(x, names(color.by), facet.by)) )
    
    
    # Facet on multiple taxa/ranks
    #________________________________________________________
    facet.by %<>% c('.taxa', .)
    .free_y <- TRUE
    
    if (length(rank) > 1)
      facet.by  %<>% c('.rank', .)
    
  })
  
  
  
  #________________________________________________________
  # Create and customize layer definitions.
  #________________________________________________________
  init_corrplot_layers(params)
  
  
  
  
  #________________________________________________________
  # y-axis title and scale
  #________________________________________________________
  set_layer(params, 'labs', y = with(params, {
    
    ifelse(
      test = is_rarefied(biom), 
      yes  = paste(rank, "Relative Abundance"), 
      no   = paste(rank, "Abundance [UNRAREFIED]"))
    
  }))
  
  if (all(params$.ggdata[['.y']] <= 1))
    set_layer(params, 'yaxis', labels = scales::percent)
  
  
  
  #________________________________________________________
  # Build the plot.
  #________________________________________________________
  fig <- params %>% 
    plot_facets() %>% 
    corrplot_stats() %>%
    plot_build()
  
  
  attr(fig, 'history') <- history
  set_cache_value(cache_file, fig)
  
  return (fig)
}



