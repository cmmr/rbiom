
#' Visualize rarefaction curves with scatterplots and trendlines.
#' 
#' @inherit documentation_test.pw_means
#' @inherit documentation_model.log
#' @inherit documentation_corrplot
#' @inherit documentation_default
#' 
#' @family rarefaction
#' @family visualization
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     rare_corrplot(hmp50, color.by="Body Site", adiv=c("shan", "otus"), facet.by="Sex")
#'     

rare_corrplot <- function (
    biom, adiv = "Shannon", depths = NULL, layers = "t", rline = TRUE,
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    test = "pw_means", model = "log", 
    p.adj = "fdr", level = 0.95, caption = FALSE, ...) {
  
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
  # Sanity Checks
  #________________________________________________________
  with(params, {
    stopifnot(is_null(depths) || is.numeric(depths))
    stopifnot(is_scalar_logical(rline) || is_scalar_integerish(rline))
    stopifnot(!is.na(rline) && !anyNA(depths))
  })
  
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  with(params, {
    validate_adiv(max = Inf)
    validate_meta_aes('color.by', null_ok = TRUE)
    validate_meta_aes('facet.by', null_ok = TRUE, max = Inf, col_type = "cat")
    validate_meta_aes('limit.by', null_ok = TRUE, max = Inf)
    
    sync_metadata()
  })
  
  
  #________________________________________________________
  # Default rarefaction depth.
  #________________________________________________________
  with(params, {
    
    if (isTRUE(rline)) {
      .ss   <- as.vector(sample_sums(biom))
      rline <- (sum(.ss) * .1) / length(.ss)
      rline <- min(.ss[.ss >= rline])
      remove(".ss")
    }
    
    if (isFALSE(rline)) rline <- NULL
  })
  
  
  #________________________________________________________
  # Pull alpha diversity metrics for each depth.
  #________________________________________________________
  with(params, {
    
    if (is_null(depths))
      depths <- "multi_even"
    
    .ggdata <- adiv_table(
      biom   = biom,
      rarefy = depths,
      adiv   = adiv,
      md     = unique(c(names(color.by), facet.by)) )
    
    .xcol <- '.depth'
    
    
    # Facet on multiple adiv metrics
    #________________________________________________________
    if (length(adiv) > 1) {
      .free_y <- TRUE
      facet.by  %<>% c('.adiv')
    }
    
  })
  
  
  
  #________________________________________________________
  # Create and customize layer definitions.
  #________________________________________________________
  init_corrplot_layers(params)
  
  
  #________________________________________________________
  # x-axis title and scale
  #________________________________________________________
  set_layer(params, 'labs',  x = "Rarefaction Depth")
  set_layer(params, 'xaxis', labels = si_units)
  if (!is_null(params$rline))
    set_layer(params, 'vline', xintercept = params$rline, color = "red", linetype="dashed")
  
  
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
  
  
  attr(fig, 'history') <- history
  set_cache_value(cache_file, fig)
  
  return (fig)
}






