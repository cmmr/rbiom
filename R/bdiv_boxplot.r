#' Visualize BIOM data with boxplots.
#' 
#' @name bdiv_boxplot
#' 
#' @param biom   A BIOM object, as returned from \link{read.biom}.
#' 
#' @param x   Metadata column name(s). Must be categorical. Prefix with 
#'        \code{==} or \code{!=} to limit comparisons to within or between
#'        groups, respectively. The default, \code{NULL} groups all distances 
#'        into a single column.
#' 
#' @param y   Beta diversity metric(s) to use. Options are \code{"Manhattan"},
#'        \code{"Euclidean"}, \code{"Bray-Curtis"}, \code{"Jaccard"}, and
#'        \code{"UniFrac"}. UniFrac requires a phylogenetic tree. Default: 
#'        \code{"Bray-Curtis"}.
#'           
#' @param layers   \code{"box"}, \code{"bar" ("r")}, \code{"violin"}, 
#'        \code{"dot"}, \code{"strip"}, \code{"crossbar"}, \code{"errorbar"}, 
#'        \code{"linerange"}, and \code{"pointrange"}. Single letter 
#'        abbreviations are also accepted. For instance, \code{c("box", "dot")} 
#'        is equivalent to \code{c("b", "d")} and \code{"bd"}.
#'        Default: \code{"rls"}.
#'        
#' @param color.by,pattern.by,shape.by,label.by,sort.by,facet.by   Metadata 
#'        column to color, pattern, shape, label, sort, and/or facet by. If 
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
#' @param weighted   When employing a beta diversity metric, use the weighted
#'        version. (Default: \code{TRUE})
#'        
#' @param tree   A phylogenetic tree for use in calculating UniFrac distance.
#'        The default, \code{NULL}, will use the BIOM object's tree.
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
#' To expand the low end of the y axis, you can set \code{y.trans = "sqrt"} or
#' \code{y.trans = "log1p"}. The former applies a square-root transformation, 
#' and the latter plots log(y + 1). Both of these methods work well with data
#' that contains zeroes. 
#' 
#' 
#' @export
#' @seealso \code{\link{stats.table}}
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     bdiv_boxplot(biom, UniFrac ~ `==Body Site`, color.by="Body Site")
#'     
#'
bdiv_boxplot <- function (
  biom, x = NULL, y = "Bray-Curtis", layers = "rls",
  color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, 
  xvals = NULL, colors = NULL, patterns = NULL, shapes = NULL, facets = NULL, 
  p.adj = "fdr", p.label = 0.05, ci = 95, xlab.angle = 'auto', 
  weighted = TRUE, tree = NULL, ...) {
  
  
  # Sanity checks
  #________________________________________________________
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  if (!(is.null(x) || all(sub("^(!=|==)", "", x) %in% metrics(biom, 'meta'))))
    stop("`x` argument to bdiv_boxplot must be metadata column name(s) or NULL")
  if (!all(y %in% metrics(biom, 'bdiv')))
    stop(
      "`y` argument to taxa_boxplot must be beta diversity metric(s): ", 
      paste(collapse = ", ", metrics(biom, 'bdiv')) )
  
  
  
  # Package and print this call
  #________________________________________________________
  params <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  history <- sprintf("bdiv_boxplot(%s)", as.args(params, fun = bdiv_boxplot))
  
  
  
  # Defaults
  #________________________________________________________
  if (is.null(x)) params[['x']] <- "."
  params[['weighted']] <- as.logical(weighted)
  
  
  
  # Use the generalized boxplot function to make the plot
  #________________________________________________________
  p <- boxplot_build(params, bdiv_boxplot_data, bdiv_boxplot_layers)
  
  
  
  # Attach history of biom modifications and this call
  #________________________________________________________
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  return (p)
}


#______________________________________________________________
# Convert biom object to a data.frame
#______________________________________________________________
bdiv_boxplot_data <- function (params) {
  
  md_cols <- NULL
  for (i in grep("\\.by$", names(params), value = TRUE))
    if (!is.null(params[[i]])) {
      md_cols %<>% c(params[[i]])
      params[[i]] %<>% sub("^[!=]=", "", .)
    }
  
  
  # Compute each bdiv metric separately - bray, etc
  #________________________________________________________
  ggdata <- plyr::ddply(
    .data      = expand.grid(
      '.metric' = params[['y']], 
      '.weight' = params[['weighted']] ), 
    .variables = c(".metric", ".weight"), 
    .fun       =  function (x) {
      beta.div(
        biom     = params[['biom']], 
        method   = x[1, '.metric'], 
        weighted = x[1, '.weight'],
        md       = md_cols,
        tree     = params[['tree']],
        long     = TRUE, 
        safe     = TRUE )
  })
  
  ggdata <- within(
    data = ggdata, 
    expr = {
      .y      <- .value
      .metric <- paste(ifelse(.weight, "Weighted", "Unweighted"), .metric)
      rm(".value", ".weight")
  })
  
  
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
  
  
  # Allow multiple `y` beta div metrics
  #________________________________________________________
  if (length(params[['y']]) > 1 || length(params[['weighted']]) > 1) {
    params[['facet.by']] %<>% c(".metric")
    if (length(params[['y']]) == 1) # "Unweighted UniFrac" => "Unweighted"
      ggdata[['.metric']] %<>% strsplit(" ", TRUE) %>% sapply(`[[`, 1)
    ggdata[['.metric']]  %<>% factor(levels = unique(ggdata[['.metric']]))
  }
  
  
  attr(ggdata, 'params') <- params
  attr(ggdata, 'xcol')   <- xcol
  attr(ggdata, 'ycol')   <- ".y"
  
  return (ggdata)
}



#______________________________________________________________
# Make bdiv-specific layer tweaks
#______________________________________________________________
bdiv_boxplot_layers <- function (layers) {
  
  # y-axis title
  #________________________________________________________
  yraw <- attr(layers, 'params', exact = TRUE)[['y']]
  
  if (length(yraw) > 1) {
    ylab <- structure(
      .Data   = "Dissimilarity (\u03B2)", 
      display = '"Dissimilarity (\\u03B2)"' )
    
  } else {
    ylab <- paste(yraw, "Distance")
  }
  
  setLayer("labs", y = ylab)
  
  
  return (layers)
}


