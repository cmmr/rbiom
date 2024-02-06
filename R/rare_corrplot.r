
#' Visualize rarefaction curves with scatterplots and trendlines.
#' 
#' @inherit documentation_corrplot
#' @inherit documentation_default
#' 
#' @family rarefaction
#' @family visualization
#' 
#' 
#' @param formula   Relationship between variables. In `rare_corrplot()`, `x` is 
#'        the sequencing depth and `y` is alpha diversity.
#'        Default: `y ~ log(x)`
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     rare_corrplot(hmp50, color.by="Body Site", adiv=c("shan", "otus"), facet.by="Sex")
#'     

rare_corrplot <- function (
    biom, adiv = "Shannon", layers = "t", rline = TRUE,
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    formula = y ~ log(x), engine = "lm", level = 0.95, 
    trans = "none", caption = TRUE, ...) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment(), ...)
  cmd    <- sprintf("rare_corrplot(%s)", as.args(params, 2, rare_corrplot))
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
    stopifnot(is_scalar_logical(rline) || is_scalar_integerish(rline))
    stopifnot(!is.na(rline))
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
    if (isTRUE(rline))  rline <- rare_suggest(biom$counts)
    if (isFALSE(rline)) rline <- NULL
  })
  
  
  #________________________________________________________
  # Select 10 depths and compute adiv metrics at each.
  #________________________________________________________
  with(params, {
    
    .ggdata <- local({
      
      upper <- fivenum(sample_sums(biom))[[4]]
      rLvls <- floor(seq(from = 5, to = upper, length.out = 10))
      
      plyr::ldply(rLvls, .id = ".depth", function (rLvl) {
        adiv_table(
          biom  = rarefy(biom, depth = rLvl),
          adiv  = adiv,
          md    = unique(c(names(color.by), facet.by)),
          trans = trans )
      }) %>% as_tibble()
    })
    
    
    # axis titles
    #________________________________________________________
    .xcol <- '.depth'
    .xlab <- "Rarefaction Depth"
    
    .ycol <- '.diversity'
    if (length(adiv) == 1)
      .ylab <- switch (
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
  
  if (!is_null(params$rline))
    set_layer(
      params     = params, 
      layer      = 'vline', 
      xintercept = params$rline, 
      color      = "red", 
      linetype   = "dashed" )
  
  
  
  
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






