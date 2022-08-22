#' Combines rare_corrplot and depths_barplot into a single figure.
#' 
#' @name rare_multiplot
#' 
#' @family plotting
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#'        
#' @param rline   Where to draw a vertical line on the plot, intended to show
#'        a particular rarefaction depth. Default: \code{TRUE} (no line).
#' 
#' @param metric   Alpha diversity metric(s) to use. Options are: 
#'        \code{"OTUs"}, \code{"Shannon"}, \code{"Chao1"}, \code{"Simpson"}, 
#'        and/or \code{"InvSimpson"}. Default: \code{"OTUs"}.
#' 
#' @param depths   Rarefaction depths to show in the plot. Passed to
#'        \link{adiv_table}. The default, \code{"multi_even"}, uses a heuristic
#'        to pick 10 evenly spaced depths.
#' 
#' @param points   Overlay a scatter plot. Default: \code{FALSE}.
#'        
#' @param color.by,facet.by   Metadata column to color and/or facet by. If that
#'        column is a \code{factor}, the ordering of levels will be maintained 
#'        in the plot. Default: \code{color.by = NULL, facet.by = NULL}
#'        
#' @param colors,facets   Names of the colors and/or facets values to use in
#'        the plot. Available values for colors are given by \code{colors()}. 
#'        Use a named character vector to map them to specific factor levels in
#'        the metadata. \code{facets} are coerced to unnamed character vectors.
#'        If the length of these vectors is less than the values present in 
#'        their corresponding metadata column, then the data set will be 
#'        subseted accordingly. Default: \code{colors = NULL, facets = NULL}
#' 
#' @param ci   The confidence interval to display around the fitted curve. Set
#'        to \code{FALSE} to hide the confidence interval. Default: \code{95}.
#'        
#' @param caption   Display information about the method used for trendline
#'        fitting beneath the plot. Default: \code{FALSE}.
#'        
#' @param labels   Show sample names under each bar. Default: \code{FALSE}.
#'        
#' @param trans   Y-axis transformation. Options are \code{"log10"} or 
#'        \code{NULL}.  Default: \code{"log10"}.
#'        
#' @param ...   Additional parameters to pass along to ggplot2 
#'        functions. Prefix a parameter name with either \code{p.}, 
#'        \code{r.}/\code{v.}, or \code{s.} to ensure it gets passed to (and
#'        only to)  \link[ggplot2]{geom_point}, \link[ggplot2]{vline}, or 
#'        \link[ggplot2]{geom_smooth}, respectively. For instance, 
#'        \code{p.size = 2} ensures only the points have their size set to 
#'        \code{2}.
#'        
#' @return A \code{ggplot2} plot. The computed data points and statistics will 
#'         be attached as \code{attr(p, 'data')} and \code{attr(p, 'stats')}, 
#'         respectively.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     rare_multiplot(hmp50, color.by="Body Site") 
#'     

rare_multiplot <- function (
    biom, rline = TRUE, metric = "OTUs", depths = NULL, points = FALSE,
    color.by = NULL, facet.by = NULL, colors = NULL, facets = NULL,
    ci = 95, caption = FALSE, labels = FALSE, trans = "log10", ...) {
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  metric <- as.vector(validate_metrics(biom, metric, "adiv", multi=TRUE))
  stopifnot(is.null(depths) || is.numeric(depths))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  params <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  history <- sprintf("rare_multiplot(%s)", as.args(params, fun = rare_multiplot))
  
  
  #________________________________________________________
  # Build the two subplots, then arrange them together.
  #________________________________________________________
  corrplot_params <- params[intersect(formalArgs(rare_corrplot),  names(params))]
  barplot_params  <- params[intersect(formalArgs(depths_barplot), names(params))]
  
  plots <- list(
    'corrplot' = do.call(rare_corrplot,  c(corrplot_params, list(...))),
    'barplot'  = do.call(depths_barplot, c(barplot_params,  list(...))) )
  
  p <- patchwork::wrap_plots(plots, ncol = 1)
  
  attr(p, 'data') <- lapply(plots, attr, which = 'data', exact = TRUE)
  attr(p, 'cmd')  <- paste(collapse = "\n\n", local({
    cmds <- sapply(seq_along(plots), function (i) {
      sub(
        x           = attr(plots[[i]], 'cmd', exact = TRUE), 
        pattern     = "ggplot(data = data", 
        replacement = sprintf("p%i <- ggplot(data = data[[%s]]", i, glue::single_quote(names(plots)[[i]])),
        fixed       = TRUE )
    })
    c(cmds, sprintf("patchwork::wrap_plots(%s, ncol = 1)", paste0(collapse = ", ", "p", seq_along(plots))))
  }))
  
  
  #________________________________________________________
  # Attach history of biom modifications and this call.
  #________________________________________________________
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  return (p)
}






