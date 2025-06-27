#' Display beta diversities in an all vs all grid.
#' 
#' @inherit documentation_heatmap
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' @family visualization
#' 
#' @param grid   Color palette name, or a list with entries for \code{label}, 
#'        \code{colors}, \code{range}, \code{bins}, \code{na.color}, and/or 
#'        \code{guide}. See the Track Definitions section for details.
#'        Default: `"devon"`
#' 
#' @param tracks   A character vector of metadata fields to display as tracks 
#'        at the top of the plot. Or, a list as expected by the `tracks` 
#'        argument of [plot_heatmap()]. Default: `NULL`
#'        
#' @param ...   Additional arguments to pass on to ggplot2::theme().
#'        For example, \code{labs.subtitle = "Plot subtitle"}.
#'         
#'         
#' @section Annotation Tracks:
#' 
#' Metadata can be displayed as colored tracks above the heatmap. Common use 
#' cases are provided below, with more thorough documentation available at 
#' https://cmmr.github.io/rbiom .
#' 
#' \preformatted{## Categorical ----------------------------
#' tracks = "Body Site"
#' tracks = list('Body Site' = "bright")
#' tracks = list('Body Site' = c('Stool' = "blue", 'Saliva' = "green"))
#' 
#' ## Numeric --------------------------------
#' tracks = "Age"
#' tracks = list('Age' = "reds")
#' 
#' ## Multiple Tracks ------------------------
#' tracks = c("Body Site", "Age")
#' tracks = list('Body Site' = "bright", 'Age' = "reds")
#' tracks = list(
#'   'Body Site' = c('Stool' = "blue", 'Saliva' = "green"),
#'   'Age'       = list('colors' = "reds") )
#' }
#' 
#' The following entries in the track definitions are understood: 
#' 
#' \describe{
#'   \item{\code{colors} - }{ A pre-defined palette name or custom set of colors to map to. }
#'   \item{\code{range} - }{ The c(min,max) to use for scale values. }
#'   \item{\code{label} - }{ Label for this track. Defaults to the name of this list element. }
#'   \item{\code{side} - }{ Options are \code{"top"} (default) or \code{"left"}. }
#'   \item{\code{na.color} - }{ The color to use for \code{NA} values. }
#'   \item{\code{bins} - }{ Bin a gradient into this many bins/steps. }
#'   \item{\code{guide} - }{ A list of arguments for guide_colorbar() or guide_legend(). }
#' }
#' 
#' All built-in color palettes are colorblind-friendly.
#' 
#' Categorical palette names: \code{"okabe"}, \code{"carto"}, \code{"r4"}, 
#' \code{"polychrome"}, \code{"tol"}, \code{"bright"}, \code{"light"}, 
#' \code{"muted"}, \code{"vibrant"}, \code{"tableau"}, \code{"classic"}, 
#' \code{"alphabet"}, \code{"tableau20"}, \code{"kelly"}, and \code{"fishy"}.
#' 
#' Numeric palette names: \code{"reds"}, \code{"oranges"}, \code{"greens"}, 
#' \code{"purples"}, \code{"grays"}, \code{"acton"}, \code{"bamako"}, 
#' \code{"batlow"}, \code{"bilbao"}, \code{"buda"}, \code{"davos"}, 
#' \code{"devon"}, \code{"grayC"}, \code{"hawaii"}, \code{"imola"}, 
#' \code{"lajolla"}, \code{"lapaz"}, \code{"nuuk"}, \code{"oslo"}, 
#' \code{"tokyo"}, \code{"turku"}, \code{"bam"}, \code{"berlin"}, 
#' \code{"broc"}, \code{"cork"}, \code{"lisbon"}, \code{"roma"}, 
#' \code{"tofino"}, \code{"vanimo"}, and \code{"vik"}.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     # Keep and rarefy the 10 most deeply sequenced samples.
#'     hmp10 <- rarefy(hmp50, n = 10)
#'     
#'     bdiv_heatmap(hmp10, tracks=c("Body Site", "Age"))
#'     
#'     bdiv_heatmap(hmp10, bdiv="uni", weighted=c(TRUE,FALSE), tracks="sex")

bdiv_heatmap <- function (
    biom, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, tracks = NULL, 
    grid = "devon", label = TRUE, label_size = NULL, rescale = "none", 
    clust = "complete", trees = TRUE, asp = 1, tree_height = 10, 
    track_height = 10, legend = "right", title = TRUE, xlab.angle = "auto", 
    underscores = FALSE, ...) {
  
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE, underscores = underscores)
  remove('underscores')
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- slurp_env(...)
  cache_file <- get_cache_file('bdiv_heatmap', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  params <- list2env(params)
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  with(params, {
    
    validate_bdiv(max = Inf)
    validate_bool("weighted", max = Inf)
    
    # validate_biom_field('order.by', null_ok = TRUE, max = Inf)
    
    if (!is_list(grid)) grid <- list(label = "Dissimilarity", colors = grid)
    # if (length(order.by) > 0) clust <- FALSE
  })
  
  
  
  #________________________________________________________
  # Handle multiple metrics/weightings with recursion.
  #________________________________________________________
  if (length(params$bdiv) > 1 || length(params$weighted) > 1) {
    
    plots <- list()
    
    for (d in params$bdiv)
      for (w in params$weighted)
        plots[[length(plots) + 1]] <- local({
          
          args <- fun_params(bdiv_heatmap, params)
          
          args[['bdiv']]     <- d
          args[['weighted']] <- w
          
          if (isTRUE(args[['title']]))
            args[['title']] <- paste(ifelse(w, "Weighted", "Unweighted"), d)
          
          do.call(bdiv_heatmap, args)
        })
    
    cmd <- current_cmd('bdiv_heatmap')
    p   <- plot_wrap(plots, cmd)
    
    set_cache_value(cache_file, p)
    return (p)
  }
  
  
  #________________________________________________________
  # Sanity Check
  #________________________________________________________
  if (params$biom$n_samples < 2)
    cli_abort("At least two samples are needed for a bdiv heatmap.")
  
  
  #________________________________________________________
  # Replace title=TRUE with a default title string.
  #________________________________________________________
  with(params, {
    
    if (isTRUE(title))
      title <- ifelse(
        test = isTRUE(weighted), 
        yes  = paste("Weighted",   bdiv, "Distance"), 
        no   = paste("Unweighted", bdiv, "Distance") )
    
    if (!(is.null(title) || (is_scalar_character(title) && !is.na(title))))
      cli_abort("title must be TRUE, NULL, or a character string, not {.type {title}}: {title}.")
  })
  
  
  
  #________________________________________________________
  # Convert `tracks` into a named list of lists.
  #________________________________________________________
  biom_tracks(params)
  
  
  
  
  #________________________________________________________
  # Matrix of samples x samples.
  #________________________________________________________
  with(params, {
    
    mtx <- bdiv_distmat(
        biom     = biom, 
        bdiv     = bdiv, 
        weighted = weighted, 
        tree     = tree ) %>%
      as.matrix()
    
    remove("biom", "bdiv", "weighted", "tree")
  })
  
  
  
  #________________________________________________________
  # Actual plotting is handled by plot_heatmap()
  #________________________________________________________
  p <- do.call(plot_heatmap, as.list(params))
  
  attr(p, 'cmd') <- current_cmd('bdiv_heatmap')
  set_cache_value(cache_file, p)
  
  return (p)
}

