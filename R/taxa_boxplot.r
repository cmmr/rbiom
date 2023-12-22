
# To-do: add metacoder for overlaying on a phylogenetic tree.

#' Visualize BIOM data with boxplots.
#' 
#' @inherit documentation_boxplot
#' @inherit documentation_default
#' 
#' @family taxa_abundance
#' @family visualization
#' 
#' @param x   A categorical metadata column name to use for the x-axis. The 
#'        default, `NULL` puts taxa along the x-axis.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     taxa_boxplot(biom, rank = c("Phylum", "Genus"), flip = TRUE)
#'     taxa_boxplot(biom, taxa = 3, layers = "ps", color.by = list("Body Site" = c('Saliva' = "blue", 'Stool' = "red")))
#'     

taxa_boxplot <- function (
    biom, x = NULL, rank = -1, taxa = 6, layers = 'bld', unc = 'singly', other = FALSE,
    color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, limit.by = NULL, 
    flip = FALSE, stripe = NULL, p.top = Inf, p.adj = 'fdr', p.label = 0.05, 
    ci = 'ci', level = 0.95, outliers = NULL, xlab.angle = 'auto', y.trans = 'sqrt', 
    caption = TRUE, ...) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment(), ...)
  cmd    <- sprintf("taxa_boxplot(%s)", as.args(params, 2, taxa_boxplot))
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
    validate_var_range('p.top', c(0, Inf))
    
    validate_meta_aes('x',          col_type = "cat", null_ok = TRUE)
    validate_meta_aes('color.by',   col_type = "cat", null_ok = TRUE)
    validate_meta_aes('pattern.by', col_type = "cat", null_ok = TRUE)
    validate_meta_aes('shape.by',   col_type = "cat", null_ok = TRUE)
    validate_meta_aes('facet.by',   col_type = "cat", null_ok = TRUE, max = Inf)
    validate_meta_aes('limit.by',                     null_ok = TRUE, max = Inf)
    
    sync_metadata()
  })
  
  
  #________________________________________________________
  # init_layers ignores formalArgs, so copy into equivalent.
  #________________________________________________________
  if (isTRUE(tolower(params$y.trans) %in% c("sqrt", "log1p")))
    params$yaxis.trans <- tolower(params$y.trans)
  
  
  
  #________________________________________________________
  # Compute taxa abundances and set up taxa/rank facet.
  #________________________________________________________
  with(params, {
    
    # Compute each rank's abundances separately.
    #________________________________________________________
    .ggdata <- taxa_table(
      biom    = biom, 
      rank    = rank, 
      taxa    = taxa, 
      md      = ".all", 
      unc     = unc, 
      other   = other )
    
    
    # Set .taxa as `x` or `facet.by`.
    #________________________________________________________
    if (nlevels(.ggdata[['.taxa']]) > 1) {
      
      if (is.null(x)) { x <- '.taxa'
      } else          { facet.by %<>% c('.taxa', .) }
    }
    
    
    # Facet on multiple ranks
    #________________________________________________________
    if (length(rank) > 1)
      facet.by %<>% c('.rank', .)
    
  })
  
  
  
  #________________________________________________________
  # Initialize the layers environment.
  #________________________________________________________
  init_boxplot_layers(params)
  
  
  
  #________________________________________________________
  # y-axis title
  #________________________________________________________
  if (any(params$.ggdata[[params$.ycol]] > 1)) {
    
    biom <- params$biom
    if (is.null(biom$depth)) { set_layer(params, 'labs', y = "Unrarefied Counts")
    } else                   { set_layer(params, 'labs', y = "Rarefied Counts") }
    
  } else {
    set_layer(params, 'labs', y = "Relative Abundance")
    set_layer(params, 'yaxis', labels = scales::percent)
  }
  
  
  
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


