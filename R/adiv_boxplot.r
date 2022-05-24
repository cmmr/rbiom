#' Visualize alpha diversity with boxplots.
#' 
#' @name adiv_boxplot
#' 
#' @param x   A BIOM object, as returned from \link{read.biom}.
#'           
#' @param layers   \code{"box"}, \code{"bar" ("r")}, \code{"violin"}, 
#'        \code{"dot"}, \code{"strip"}, \code{"crossbar"}, \code{"errorbar"}, 
#'        \code{"linerange"}, and \code{"pointrange"}. Single letter 
#'        abbreviations are also accepted. For instance, \code{c("box", "dot")} 
#'        is equivalent to \code{c("b", "d")} and \code{"bd"}.
#'        Default: \code{"rls"}.
#'        
#' @param color.by,pattern.by,shape.by,facet.by   Metadata 
#'        column to color, pattern, shape, and/or facet by. If 
#'        that column is a \code{factor}, the ordering of levels will be 
#'        maintained in the plot.
#'        
#' @param colors,patterns,shapes,facets,xvals   Names of the colors, patterns,
#'        shapes, facets, and/or x values to use in the plot. Available values 
#'        for colors, patterns, and shapes are given by \code{colors()}, 
#'        \code{gridpattern::names_magick}, and \code{0:25},
#'        respectively. Use a named character vector to map them to specific
#'        factor levels in the metadata. \code{facets} and \code{xvals} are
#'        coerced to unnamed character vectors. If the length of these vectors
#'        is less than the values present in their corresponding metadata 
#'        column, then the data set will be subseted accordingly.
#'
#' @param p.adj   Method to use for multiple comparisons adjustment of p-values.
#'        Run \code{p.adjust.methods} for a list of available options.
#'        (Default: \code{fdr})
#'        
#' @param p.label   Minimum adjusted p-value to display on the plot with a bracket.
#'        Set to \code{Inf} to display all p-values, or \code{-Inf} for no brackets.
#'        If a numeric vector with more than one value is provided, they will be
#'        used as breaks for asterisk notation. For ordinations, \code{p.label} 
#'        applies to biplot taxa. (Default: \code{0.05})
#'     
#' @param ci   How to calculate min/max of the \bold{crossbar}, 
#'        \bold{errorbar}, \bold{linerange}, and \bold{pointrange} layers.
#'        Provide a number between 75 and 100 to define a confidence interval's
#'        confidence level, commonly 95 or 97.5. Other options are: 
#'        \bold{range}, 
#'        \bold{sd} (standard deviation), 
#'        \bold{se} (standard error), and 
#'        \bold{mad} (median absolute deviation). 
#'        The center mark of \code{crossbar} and \code{pointrange} represents
#'        the mean, except for \bold{mad} in which case it represents
#'        the median. Trendlines require a confidence interval value. 
#'        Set to \code{NULL} to disable. Default: \code{95}
#'        
#' @param xlab.angle   How to rotate the tick labels on the x-axis. 
#'        \bold{'auto'} (the default), automatically selects a rotation value. 
#'        \bold{0}, \bold{30}, and \bold{90} sets the angle to horizontal, 
#'        angled, and vertical, respectively.
#'        
#' @param safe   If \code{FALSE}, data.frame column names such as 
#'        \code{".metric"} will be auto-converted to \code{"Metric"} to improve
#'        human-readability. Conversion if aborted if a conflict is found with
#'        a metadata column name. (Default: \code{FALSE})
#'        
#' @param ...   Parameters are matched to formal arguments of ggplot2
#'        functions. Prefixing parameter names with a layer name ensures that
#'        a particular parameter is passed to, and only to, that layer. For
#'        instance, \code{dot.size = 2} or \code{d.size = 2} ensures only the 
#'        dotplot layer has its size set to \code{2}.
#'        
#' @return A \code{ggplot2} plot. The computed data points and statistics will 
#'         be attached as \code{attr(p, 'data')} and \code{attr(p, 'stats')}, 
#'         respectively.
#'         
#' Shapes can also be given as their string values, defined in pch_table here:
#' https://github.com/tidyverse/ggplot2/blob/master/R/geom-point.r . Note that
#' some shapes have a colored outline given by `color`, some are filled with 
#' `color`, and some are outlined in `color` and filled with `fill`. See
#' https://blog.albertkuo.me/post/point-shape-options-in-ggplot/ for details.
#' 
#' 
#' @export
#' @seealso \code{\link{stats.table}}
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     adiv_boxplot(biom, x = "Body Site", y = "Shannon")
#'     adiv_boxplot(biom, x = "Sex", y = c("OTUs", "Shannon"), layers="b", color.by="Body Site", scales="free")
#'     adiv_boxplot(biom, x = "Body Site", y = "Simpson", layers="p", color.by="Sex", xlab.angle=30)
#'     

adiv_boxplot <- function (
    biom, x = ".", y = "Shannon", layers = "rls",
    color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, 
    xvals = NULL, colors = NULL, patterns = NULL, shapes = NULL, facets = NULL, 
    p.adj = "fdr", p.label = 0.05, ci = 95, xlab.angle = 'auto', safe = FALSE, ...) {
  
  
  # Sanity checks
  #________________________________________________________
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  if (!all(x %in% c(".", metrics(biom, 'meta'))))
    stop("`x` argument to adiv_boxplot must be metadata column name(s) or '.'")
  if (!all(y %in% metrics(biom, 'adiv')))
    stop(
      "`y` argument to adiv_boxplot must be alpha diversity metric name(s): ", 
      paste(collapse = ", ", metrics(biom, 'adiv')) )
  
  
  
  # Use the generalized boxplot function to make the plot
  #________________________________________________________
  params <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  p <- boxplot_build(params, adiv_boxplot_data, adiv_boxplot_layers)
  
  
  
  # Attach history of biom modifications and this call
  #________________________________________________________
  history <- sprintf("adiv_boxplot(%s)", as.args(params, fun = adiv_boxplot))
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  return (p)
  
}


#______________________________________________________________
# Convert biom object to a data.frame
#______________________________________________________________
adiv_boxplot_data <- function (params) {
  
  ggdata <- alpha.div(
    biom    = params[['biom']], 
    metrics = params[['y']],
    md      = TRUE,
    long    = TRUE, 
    safe    = TRUE )
  
  names(ggdata)[which(names(ggdata) == ".value")] <- ".y"
  
  
  # Allow multiple `x` metadata fields
  #________________________________________________________
  xcol <- params[['xval.by']]
  if (length(xcol) > 1) {
    
    ggdata <- plyr::ldply(xcol, function (x) {
      df <- ggdata
      df[['.xfield']] <- x
      df[['.x']]      <- ggdata[[x]]
      return (df)
    })
    ggdata[['.xfield']] %<>% factor(levels = xcol)
    ggdata[['.x']]      %<>% factor(
      levels = unique(unlist(sapply(xcol, function (x) {
        levels(ggdata[[x]])
      }))))
    
    params[['facet.by']] %<>% c(".xfield")
    params[['xval.by']] <- xcol <- ".x"
  }
  
  
  # Allow multiple `y` alpha div metrics
  #________________________________________________________
  if (length(params[['y']]) > 1) {
    params[['facet.by']] %<>% c(".metric")
    ggdata[['.metric']]  %<>% factor(levels = params[['y']])
  }
  
  
  
  
  attr(ggdata, 'params') <- params
  attr(ggdata, 'xcol')   <- xcol
  attr(ggdata, 'ycol')   <- ".y"
  
  return (ggdata)
}



#______________________________________________________________
# Make adiv-specific layer tweaks
#______________________________________________________________
adiv_boxplot_layers <- function (layers) {
  
  # y-axis title
  #________________________________________________________
  yraw <- attr(layers, 'params', exact = TRUE)[['y']]
  
  if (length(yraw) > 1) {
    ylab <- structure(
      .Data   = "Diversity (\u03B1)", 
      display = '"Diversity (\\u03B1)"' )
    
    attr(layers, 'free_y') <- TRUE
    
  } else {
    ylab <- switch (
      EXPR    = yraw,
      'OTUs'  = "Observed OTUs",
      'Depth' = "Sequencing Depth",
      paste(yraw, "Diversity")
    )
  }
  
  setLayer("labs", y = ylab)
  
  
  return (layers)
}
