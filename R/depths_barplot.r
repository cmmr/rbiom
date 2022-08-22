#' Visualize the number of observations per sample.
#' 
#' @name depths_barplot
#' 
#' @family plotting
#' 
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM}
#'        object, as returned from \link{read_biom}. For matrices, the rows and
#'        columns are assumed to be the taxa and samples, respectively.
#'        
#' @param rline   Where to draw a horizontal line on the plot, intended to show
#'        a particular rarefaction depth. Set to \code{TRUE} to show an 
#'        auto-selected rarefaction depth or \code{NULL} to not show a line.
#'        Default: \code{TRUE}.
#'        
#' @param labels   Show sample names under each bar. Default: \code{TRUE}.
#'        
#' @param trans   Y-axis transformation. Options are \code{"log10"} or 
#'        \code{NULL}.  Default: \code{"log10"}.
#'        
#' @param ...   Additional parameters to pass along to ggplot2 
#'        functions. Prefix a parameter name with either \code{p.}, 
#'        \code{r.}/\code{h.}, or \code{s.} to ensure it gets passed to (and
#'        only to)  \link[ggplot2]{geom_point}, \link[ggplot2]{hline}, or 
#'        \link[ggplot2]{geom_smooth}, respectively. For instance, 
#'        \code{p.size = 2} ensures only the points have their size set to 
#'        \code{2}.
#'        
#' @return A \code{ggplot2} plot. The computed data points will be attached as 
#'         \code{attr(p, 'data')}.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     depths_barplot(hmp50, rline=TRUE) 
#'     

depths_barplot <- function (biom, rline = TRUE, labels = TRUE, trans = "log10", ...) {
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  stopifnot(is.logical(rline) || is.null(rline) || is.numeric(rline))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  params <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  history <- sprintf("depths_barplot(%s)", as.args(params, fun = depths_barplot))
  
  
  #________________________________________________________
  # Find the number of observations per sample.
  #________________________________________________________
  ss  <- if (is(biom, "simple_triplet_matrix")) { slam::col_sums(biom)
  } else if (is(biom, "BIOM"))                  { sample_sums(biom)
  } else if (is(biom, "matrix"))                { colSums(biom)
  } else { stop("biom must be a matrix, simple_triplet_matrix, or BIOM object.") }
  ss <- sort(ss)
  
  
  #________________________________________________________
  # Default rarefaction depth.
  #________________________________________________________
  if (isFALSE(rline)) rline <- NULL
    
  if (isTRUE(rline)) {
    rline <- (sum(ss) * .1) / length(ss)
    rline <- min(ss[ss >= rline])
  }
  
  
  #________________________________________________________
  # Split the counts into kept/dropped.
  #________________________________________________________
  data <- if (is.null(rline)) {
    data.frame(
      check.names = FALSE,
      '.x'    = factor(names(ss), levels = names(ss)), 
      '.xmin' = seq_along(ss) - 0.4,
      '.xmax' = seq_along(ss) + 0.4,
      '.ymin' = 1,
      '.ymax' = as.vector(ss) )
    
  } else {
    data.frame(
      check.names = FALSE,
      '.x'     = factor(names(ss), levels = names(ss)),
      '.group' = factor(rep(c("Excluded", "Retained"), each = length(ss))),
      '.xmin'  = seq_along(ss) - 0.4,
      '.xmax'  = seq_along(ss) + 0.4,
      '.ymin'  = as.vector(c(ifelse(ss  < rline, 1, rline), ifelse(ss < rline, 0, 1))),
      '.ymax'  = as.vector(c(ifelse(ss == rline, 0, ss),    ifelse(ss < rline, 0, rline))) )
  }
  data %<>% subset(.ymin > 0 & .ymax > 0)
  
  
  #________________________________________________________
  # Initialize the `layers` object.
  #________________________________________________________
  layers <- structure(
    list(),
    'data'     = data,
    'params'   = params,
    'function' = depths_barplot,
    'xcol'     = ".sample",
    'xmode'    = "factor" )
  
  initLayer(c('ggplot', 'xaxis', 'yaxis', 'theme_bw'))
  
  
  #________________________________________________________
  # Add default layer parameters.
  #________________________________________________________
  setLayer("rect",  mapping = .qw(xmin, xmax, ymin, ymax), color=NA)
  setLayer("labs",  x = "Sample", y = "Sequencing Depth")
  setLayer("yaxis", expand = c(0,0))
  setLayer("theme", panel.grid.major.x = element_blank())
  
  if (identical(trans, "log10")) {
    setLayer("yaxis", c(loglabels(ss), trans="log10"))
  } else {
    setLayer("yaxis", labels = si_units)
  }
  
  if (isTRUE(labels)) {
    setLayer("point", mapping     = list(x = ".x"), y = 1, alpha = 0)
    setLayer("theme", axis.text.x = element_text(angle=-30, vjust=1, hjust=0) )
    setLayer("labs",  x = NULL)
  }
  
  if (!is.null(rline)) {
    setLayer("rect", 'mapping|fill' = ".group")
    setLayer("labs",  fill = "Reads")
    setLayer("hline", yintercept = rline, color = "red", linetype="dashed")
  }
  
  
  
  # Convert layer definitions into a plot.
  #________________________________________________________
  p <- plot_build(layers)
  
  
  # Attach history of biom modifications and this call.
  #________________________________________________________
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  return (p)
  
  
}






