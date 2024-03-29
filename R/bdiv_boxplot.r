#' Visualize BIOM data with boxplots.
#' 
#' @inherit documentation_cmp
#' @inherit documentation_boxplot
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' @family visualization
#' 
#' @export
#' 
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     bdiv_boxplot(biom, x="==Body Site", bdiv="UniFrac", color.by="Body Site")
#'     
#'
bdiv_boxplot <- function (
  biom, x = NULL, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, layers = "x",
  color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, limit.by = NULL, 
  within = NULL, between = NULL, trans = "none", flip = FALSE, stripe = NULL, p.adj = "fdr", 
  p.label = 0.05, ci = "ci", level = 0.95, outliers = NULL, xlab.angle = 'auto', caption = TRUE, ...) {
  
  
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)
  
  params <- eval_envir(environment(), ...)
  cmd    <- sprintf("bdiv_boxplot(%s)", as.args(params, 2, bdiv_boxplot))
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
    validate_bool('weighted', max = Inf)
    
    validate_meta_aes('x',          col_type = "cat", null_ok = TRUE, cmp = TRUE)
    validate_meta_aes('color.by',   col_type = "cat", null_ok = TRUE, cmp = TRUE)
    validate_meta_aes('pattern.by', col_type = "cat", null_ok = TRUE, cmp = TRUE)
    validate_meta_aes('shape.by',   col_type = "cat", null_ok = TRUE, cmp = TRUE)
    validate_meta_aes('facet.by',   col_type = "cat", null_ok = TRUE, cmp = TRUE, max = Inf)
    validate_meta_aes('limit.by',                     null_ok = TRUE, cmp = TRUE, max = Inf)
    
    # Validates and appends to `within` and `between`.
    validate_meta_cmp(c('x', 'color.by', 'pattern.by', 'shape.by', 'facet.by', 'limit.by'))
    
    sync_metadata()
  })
  
  
  
  #________________________________________________________
  # Compute beta diversity and set up bdiv facet.
  #________________________________________________________
  with(params, {
    
    # Compute each bdiv distance metric separately: bray, etc
    #________________________________________________________
    .ggdata <- bdiv_table(
        biom     = biom,
        bdiv     = bdiv,
        weighted = weighted,
        tree     = tree, 
        md       = ".all",
        within   = within,
        between  = between, 
        trans    = trans ) %>%
      within({
        .bdiv <- paste(ifelse(.weighted, 'Weighted', 'Unweighted'), .bdiv)
        .bdiv <- factor(.bdiv, levels = unique(.bdiv))
        remove(".weighted") })
    
    
    # axis titles
    #________________________________________________________
    .ylab <- paste(bdiv, "Distance")
    
    
    # Facet on multiple bdiv metrics
    #________________________________________________________
    if (length(bdiv) > 1 || length(weighted) > 1) {
      facet.by %<>% c(".bdiv")
      .ylab <- aa("\u03B2 Dissimilarity", display = '"\\u03B2 Dissimilarity"')
    }
  })
  
  
  
  #________________________________________________________
  # Initialize the layers environment.
  #________________________________________________________
  init_boxplot_layers(params)
  
  
  #________________________________________________________
  # Note any within/between groups in the plot caption.
  #________________________________________________________
  within  <- paste(collapse = ", ", params$within)
  between <- paste(collapse = ", ", params$between)
  
  caption <- paste(collapse = "\n", c(
    if (has_layer(params, 'labs')) params$layers[['labs']][['caption']] else NULL,
    if (nzchar(within))            paste("Within groups:",  within)     else NULL,
    if (nzchar(between))           paste("Between groups:", between)    else NULL ))
  
  if (nchar(caption))
    set_layer(params, 'labs', caption = caption)
  
  
  
  
  #________________________________________________________
  # Use the generalized boxplot functions to make the plot
  #________________________________________________________
  fig <- params %>%
    boxplot_facets() %>%
    boxplot_stats() %>%
    plot_build()
  
  
  
  attr(fig, 'cmd') <- cmd
  set_cache_value(cache_file, fig)
  
  return (fig)
}
