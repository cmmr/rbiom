
# To-do: add metacoder for overlaying on a phylogenetic tree.

#' Visualize BIOM data with boxplots.
#' 
#' @name taxa_boxplot
#' 
#' @family plotting
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#' 
#' @param x   A categorical metadata column name to use for the x-axis. The 
#'        default, \code{".taxa"} puts the taxa along the x-axis.
#'        
#' @param rank   What rank of taxa to display. E.g. "Phylum", "Genus", etc. Run
#'        \code{taxa_ranks()} to see all options for a given BIOM object. The
#'        default, \code{NULL}, selects the lowest level.
#'        
#' @param taxa   Which taxa to display. An integer value will show the top n
#'        most abundant taxa. A value 0 <= n < 1 will show any taxa with that 
#'        mean abundance or greater (e.g. 0.1). A character vector of
#'        taxon names will show only those taxa. Default: \code{5}.
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
#' @param p.top   Only display taxa with the most significant differences in 
#'        abundance. If \code{p.top} is >= 1, then the \code{p.top} most 
#'        significant taxa are displayed. If \code{p.top} is less than one, all 
#'        taxa with an adjusted p-value <= \code{p.top} are displayed. 
#'        Recommended to be used in combination with the \code{taxa} parameter 
#'        to set a lower bound on the mean abundance of considered taxa. 
#'        Default: \code{Inf}
#'
#' @param p.adj   Method to use for multiple comparisons adjustment of p-values.
#'        Run \code{p.adjust.methods} for a list of available options.
#'        Default: \code{fdr}
#'        
#' @param p.label   Minimum adjusted p-value to display on the plot with a 
#'        bracket. Set to \code{Inf} to display all p-values, or \code{-Inf} 
#'        for no brackets. If a numeric vector with more than one value is 
#'        provided, they will be used as breaks for asterisk notation. 
#'        Default: \code{0.05}
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
#' @param other   Should non-selected taxa be displayed as an "Other" group?
#'        Default: \code{FALSE}
#'        
#' @param safe   If \code{FALSE}, data.frame column names such as 
#'        \code{".metric"} will be auto-converted to \code{"Metric"} to improve
#'        human-readability. Conversion if aborted if a conflict is found with
#'        a metadata column name. Default: \code{FALSE}
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
#' @seealso \code{\link{stats_table}}
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     taxa_boxplot(biom, rank = c("Phylum", "Genus"))
#'     
#'
taxa_boxplot <- function (
    biom, x = ".taxa", rank = NULL, taxa = 5, layers = "rls",
    color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, 
    xvals = NULL, colors = NULL, patterns = NULL, shapes = NULL, facets = NULL, 
    p.top = Inf, p.adj = "fdr", p.label = 0.05, ci = 95, xlab.angle = 'auto', ...) {
  
  
  # Sanity checks
  #________________________________________________________
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  if (!all(x %in% c(".taxa", metrics(biom, 'meta'))))
    stop("`x` argument to taxa_boxplot must be metadata column name(s) or '.taxa'")
  if (is.null(rank)) rank <- tail(taxa_ranks(biom), 1)
  if (!all(rank %in% metrics(biom, 'rank')))
    stop(
      "`rank` argument to taxa_boxplot must be taxa rank(s): ", 
      paste(collapse = ", ", metrics(biom, 'rank')) )
  
  
  
  # Use the generalized boxplot function to make the plot
  #________________________________________________________
  params <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  p <- boxplot_build(params, taxa_boxplot, taxa_boxplot_data, taxa_boxplot_layers)
  
  
  
  # Attach history of biom modifications and this call
  #________________________________________________________
  history <- sprintf("taxa_boxplot(%s)", as.args(params, fun = taxa_boxplot))
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  return (p)
  
}


#______________________________________________________________
# Convert biom object to a data.frame
#______________________________________________________________
taxa_boxplot_data <- function (params) {
  
  biom  <- params[['biom']]
  ranks <- params[['rank']]
  taxa  <- params[['taxa']]
  
  # Convert abundance spec to taxa names
  if (is.numeric(taxa))
    taxa <- as.vector(sapply(ranks, function (rank) {
      means <- taxa_means(as_percent(biom), rank)
      if (taxa < 1) { return (names(which(means >= taxa)))
      } else        { return (head(names(means), taxa)) }
    }))
  
  if (is_rarefied(biom))
    biom <- as_percent(biom)
  
  ggdata <- plyr::ldply(ranks, function (rank) {
    
    df <- taxa_table(
      biom = biom,
      rank = rank,
      md   = TRUE,
      safe = TRUE )
    
    df <- df[df[['.taxa']] %in% taxa,,drop=FALSE]
    if (nrow(df) == 0) return (NULL)
    
    df[['.rank']] <- rank
    df[['.y']]    <- df[['.value']]
    return (df)
  })
  
  
  ggdata[['.taxa']] %<>% factor(levels = taxa)
  ggdata[['.rank']] %<>% factor(levels = ranks)
  
  if (params[['x']] != ".taxa") { params[['facet.by']] %<>% c(".taxa", .)
  } else if (length(ranks) > 1) { params[['facet.by']] %<>% c(".rank", .) }
  
  
  attr(ggdata, 'params') <- params
  attr(ggdata, 'xcol')   <- params[['x']]
  attr(ggdata, 'ycol')   <- ".y"
  
  return (ggdata)
}


#______________________________________________________________
# Make taxa-specific layer tweaks
#______________________________________________________________
taxa_boxplot_layers <- function (layers) {
  
  setLayer("labs", y = "Relative Abundance")
  setLayer("yaxis", labels = scales::percent)
  
  return (layers)
  
}
