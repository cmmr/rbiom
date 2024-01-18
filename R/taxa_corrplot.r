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
#'     biom <- rarefy(hmp50)
#'     
#'     taxa_corrplot(biom, "BMI", color.by="Body Site", taxa = 4) 
#'     

taxa_corrplot <- function (
    biom, x, rank = -1, taxa = 6, layers = "t",
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    test = "trends", model = "lm", trans = "rank", 
    p.top = Inf, p.adj = "fdr", level = 0.95, caption = TRUE, ...) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment(), ...)
  cmd    <- sprintf("taxa_corrplot(%s)", as.args(params, 2, taxa_corrplot))
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
    
    
    # Compute each rank separately: phylum, genus, etc
    #________________________________________________________
    .ggdata <- taxa_table(
      biom  = biom, 
      rank  = rank, 
      taxa  = taxa, 
      md    = unique(c(x, names(color.by), facet.by)), 
      trans = trans )
    
    
    # axis titles
    #________________________________________________________
    .ylab <- {
      if      (eq(trans, 'percent')) { "Relative Abundance" }
      else if (is.null(biom$depth))  { "Unrarefied Counts"  }
      else                           { "Rarefied Counts"    }
    }
    
    
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



