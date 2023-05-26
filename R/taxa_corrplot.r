#' Visualize taxa abundance with scatterplots and trendlines.
#' 
#' @name taxa_corrplot
#' 
#' @inherit adiv_corrplot params sections return
#' 
#' @family plotting
#' 
#'        
#' @param rank   What rank of taxa to display. E.g. "Phylum", "Genus", etc. 
#'        Run \code{taxa_ranks()} to see all options for a given BIOM object. 
#'        The default, \code{NULL}, selects the lowest level.
#'        
#' @param taxa   Which taxa to display. An integer value will show the top n
#'        most abundant taxa. A value 0 <= n < 1 will show any taxa with that 
#'        mean abundance or greater (e.g. 0.1). A character vector of
#'        taxon names will show only those taxa. Default: \code{6}.
#'        
#' @param ...   Additional parameters to pass along to ggplot2
#'        functions. Prefix a parameter name with either \code{p.} or \code{s.}
#'        to ensure it gets passed to (and only to) 
#'        \link[ggplot2]{geom_point} or 
#'        \link[ggplot2]{geom_smooth}, respectively. For instance, 
#'        \code{p.size = 2} ensures only the points have their size set to 
#'        \code{2}.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     taxa_corrplot(rarefy(hmp50), "BMI", color.by="Body Site") 
#'     

taxa_corrplot <- function (
    biom, x, rank = NULL, taxa = 6, points = FALSE,
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    model = "linear", ci = 95, ...) {
  
  with_cache("taxa_corrplot", environment(), list(...), local({
    
    
    #________________________________________________________
    # Record the function call in a human-readable format.
    #________________________________________________________
    params  <- as.list(parent.env(environment()))
    history <- attr(biom, 'history')
    history %<>% c(sprintf("taxa_corrplot(%s)", as.args(params, fun = taxa_corrplot)))
    remove(list = setdiff(ls(), c("params", "history")))
    
    
    #________________________________________________________
    # Sanity Checks
    #________________________________________________________
    params %<>% within({
      if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
      rank %<>% validate_arg(biom, 'rank', n = 1, default = tail(c('OTU', taxa_ranks(biom)), 1))
    })
    
    
    #________________________________________________________
    # Subset biom by requested metadata and aes.
    #________________________________________________________
    params %<>% metadata_params(contraints = list(
      x        = list(n = 1, col_type = "num"),
      color.by = list(n = 0:1),
      facet.by = list(col_type = "cat"),
      limit.by = list() ))
    
    
    #________________________________________________________
    # Compute taxa abundance values.
    #________________________________________________________
    
    if (is_rarefied(params[['biom']]))
      params[['biom']] %<>% as_percent()
    
    data <- taxa_table(
      biom = params[['biom']],
      rank = params[['rank']],
      taxa = params[['taxa']],
      md   = unique(c(
        params[['x']], 
        names(params[['color.by']]), 
        params[['facet.by']] )),
      safe = TRUE )
    
    params[['facet.by']] %<>% c(".taxa")
    
    
    #________________________________________________________
    # Initialize the `layers` object.
    #________________________________________________________
    layers <- structure(
      list(),
      'data'     = data,
      'params'   = params,
      'function' = taxa_corrplot,
      'xcol'     = params[['x']],
      'ycol'     = ".value",
      'xmode'    = "numeric" )
    
    initLayer(c('ggplot', 'smooth', 'labs', 'theme_bw'))
    
    if (!is_null(params[['color.by']])) initLayer("color")
    if (isTRUE(params[['points']]))     initLayer("point")
    
    
    #________________________________________________________
    # Add default layer parameters.
    #________________________________________________________
    setLayer("ggplot", mapping = list(x = params[['x']], y = ".value"))
    
    setLayer("labs", x = params[['x']], y = local({
      ylab <- paste(params[['rank']], "Abundance")
      if (!is_rarefied(params[['biom']])) ylab %<>% paste("[UNRAREFIED]")
      return (ylab) }))
    
    setLayer("facet", scales = "free_y")
    
    
    
    
    # Convert layer definitions into a plot.
    #________________________________________________________
    p <- corrplot_build(layers)
    
    attr(p, 'history') <- history
    
    
    return (p)
    
  }))
}






