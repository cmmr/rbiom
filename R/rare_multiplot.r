#' Combines rare_corrplot and rare_barplot into a single figure.
#' 
#' @name rare_multiplot
#' 
#' @inherit rare_corrplot params
#' @inherit adiv_corrplot params sections return
#' 
#' @family visualization
#'        
#' @param labels   Show sample names under each bar. Default: \code{FALSE}.
#'        
#' @param trans   Y-axis transformation. Options are \code{"log10"} or 
#'        \code{NULL}.  Default: \code{"log10"}.
#'        
#' @param ...   Additional parameters to pass along to ggplot2 
#'        functions. Prefix a parameter name with either \code{t.} or 
#'        \code{s.}/\code{pt.} to ensure it gets passed to (and only to)
#'        \link[ggplot2]{geom_smooth} or \link[ggplot2]{geom_point}, 
#'        respectively. For instance, \code{s.size = 2} ensures only the points 
#'        have their size set to \code{2}.
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
    biom, adiv = "OTUs", depths = NULL, layers = "t", rline = TRUE,
    color.by = NULL, facet.by = NULL, limit.by = NULL, 
    model = "log", stats = "emtrends", p.adj = "fdr", ci = 95, 
    caption = FALSE, labels = FALSE, trans = "log10", ...) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(c(as.list(environment()), list(...)), eval)
  cache_file <- get_cache_file("rare_multiplot", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  dots <- list(...)
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  arg_str <- as.args(params, fun = rare_multiplot, indent = 2)
  history <- paste0(collapse = "\n", c(
    attr(biom, 'history', exact = TRUE),
    sprintf("fig  <- rare_multiplot(%s)", arg_str) ))
  remove(list = setdiff(ls(), c("params", "history", "cache_file", "dots")))
  
  
  #________________________________________________________
  # Build the two subplots, then arrange them together.
  #________________________________________________________
  corrplot_params <- params[intersect(formalArgs(rare_corrplot),  names(params))]
  barplot_params  <- params[intersect(formalArgs(rare_barplot), names(params))]
  
  corrplot <- do.call(rare_corrplot,  c(corrplot_params, dots))
  barplot  <- do.call(rare_barplot, c(barplot_params,  dots))
  plots    <- list(corrplot = corrplot, barplot = barplot)
  
  p <- patchwork::wrap_plots(plots, ncol = 1)
  
  attr(p, 'stats') <- attr(corrplot, 'stats', exact = TRUE)
  attr(p, 'data')  <- lapply(plots, attr, which = 'data', exact = TRUE)
  attr(p, 'cmd')   <- paste(collapse = "\n\n", local({
    cmds <- sapply(seq_along(plots), function (i) {
      sub(
        x           = attr(plots[[i]], 'cmd', exact = TRUE), 
        pattern     = "ggplot(data = data", 
        replacement = sprintf("p%i <- ggplot(data = data[[%s]]", i, single_quote(names(plots)[[i]])),
        fixed       = TRUE )
    })
    c(cmds, sprintf("patchwork::wrap_plots(%s, ncol = 1)", paste0(collapse = ", ", "p", seq_along(plots))))
  }))
  
  
  #________________________________________________________
  # Attach history of biom modifications and this call.
  #________________________________________________________
  attr(p, 'history') <- history
  
  
  set_cache_value(cache_file, p)
  return (p)
}






