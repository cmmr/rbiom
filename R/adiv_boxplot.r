#' Visualize alpha diversity with boxplots.
#' 
#' @name adiv_boxplot
#' 
#' @family plotting
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#' 
#' @param x   A categorical metadata column name to use for the x-axis. The 
#'        default, \code{NULL}, groups all samples into a single category. 
#' 
#' @param metric   Alpha diversity metric(s) to use. Options are: \code{"OTUs"}, 
#'        \code{"Shannon"}, \code{"Chao1"}, \code{"Simpson"}, and/or 
#'        \code{"InvSimpson"}. Default: \code{"Shannon"}.
#'           
#' @param layers   \code{"box" ("x")}, \code{"bar"}, \code{"violin"}, 
#'        \code{"dot"}, \code{"strip"}, \code{"crossbar"}, \code{"errorbar"}, 
#'        \code{"linerange"}, and \code{"pointrange"}. Single letter 
#'        abbreviations are also accepted. For instance, \code{c("box", "dot")} 
#'        is equivalent to \code{c("x", "d")} and \code{"xd"}.
#'        See \code{vignette("boxplots")} for examples of each.
#'        Default: \code{"lsb"}.
#'        
#' @param color.by,pattern.by,shape.by,facet.by,limit.by   Metadata columns to 
#'        use for aesthetics and partitioning. See below for details.
#'        Default: \code{NULL}
#'
#' @param flip   Transpose the axes, so that taxa are present as rows instead
#'        of columns. Default: \code{FALSE}
#'
#' @param stripe   Shade every other x position. Default: \emph{same as flip}
#'
#' @param p.adj   Method to use for multiple comparisons adjustment of p-values.
#'        Run \code{p.adjust.methods} for a list of available options.
#'        (Default: \code{fdr})
#'        
#' @param p.label   Minimum adjusted p-value to display on the plot with a 
#'        bracket.
#'        \itemize{
#'          \item{\code{p.label = 0.05} - }{ Use this specific value. }
#'          \item{\code{p.label = TRUE} - }{ equivalent to \code{p.label = 0.05} }
#'          \item{\code{p.label = FALSE} - }{ do not show any stats on the plot }
#'          \item{\code{p.label = Inf} - }{ display all p-values }
#'          \item{\code{p.label = NULL} - }{ do not calculate stats }
#'        }
#'        If a numeric vector with more than one value is 
#'        provided, they will be used as breaks for asterisk notation.
#'        Default: \code{TRUE}
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
#' @param outliers   Show boxplot outliers? \code{TRUE} to always show. 
#'        \code{FALSE} to always hide. \code{NULL} to only hide them when
#'        overlaying a dot or strip chart.  Default: \code{NULL}
#'        
#' @param xlab.angle   How to rotate the tick labels on the x-axis. 
#'        \bold{'auto'} (the default), automatically selects a rotation value. 
#'        \bold{0}, \bold{30}, and \bold{90} sets the angle to horizontal, 
#'        angled, and vertical, respectively.
#'        
#' @param ...   Parameters are matched to formal arguments of ggplot2
#'        functions. Prefixing parameter names with a layer name ensures that
#'        a particular parameter is passed to, and only to, that layer. For
#'        instance, \code{dot.size = 2} or \code{d.size = 2} ensures only the 
#'        dotplot layer has its size set to \code{2}. The special prefix
#'        \code{pt.} will control both the dot and strip layers.
#'        
#' @return A \code{ggplot2} plot. The computed data points and statistics will 
#'         be attached as \code{attr(p, 'data')} and \code{attr(p, 'stats')}, 
#'         respectively.
#' 
#' 
#' @section Aesthetics and Partitions:
#' 
#' Metadata can be used to flexibly subset, partition, and apply aesthetics 
#' when creating a plot. Common use cases are provided below. More thorough 
#' documentation is available at \url{https://cmmr.github.io/rbiom}.
#' 
#' 
#' \preformatted{  ## Colors ----------------------------
#'   color.by = "Body Site"
#'   color.by = list('Body Site' = "bright")
#'   color.by = list('Body Site' = c("Stool", "Saliva"))
#'   color.by = list('Body Site' = list('values' = c("Stool", "Saliva"), 'colors' = "bright"))
#'   color.by = list('Body Site' = c('Stool' = "blue", 'Saliva' = "green"))
#'   
#'   ## Patterns --------------------------
#'   pattern.by = "Body Site"
#'   pattern.by = list('Body Site' = c("Stool", "Saliva"))
#'   pattern.by = list('Body Site' = c('Stool' = "left45", 'Saliva' = "hs_cross"))
#'   
#'   ## Shapes ----------------------------
#'   shape.by = "Body Site"
#'   shape.by = list('Body Site' = c("Stool", "Saliva"))
#'   shape.by = list('Body Site' = c('Stool' = 7, 'Saliva' = 8))
#'   
#'   ## Facets ----------------------------
#'   facet.by = "Body Site"
#'   facet.by = c("Body Site", "Sex")
#'   facet.by = list('Body Site' = c("Stool", "Saliva"), "Sex")
#'   
#'   ## Limits ----------------------------
#'   limit.by = list('Sex' = "Male", 'Age' = c(20,40))
#'   limit.by = list('Body Site' = c("Saliva", "Anterior nares"), 'Age' = c(NA,35))
#' }
#' 
#' \itemize{
#'   \item{\code{color.by} - }{A categorical metadata column. (Max 1)}
#'   \item{\code{pattern.by} - }{A categorical metadata column. (Max 1)}
#'   \item{\code{shape.by} - }{A categorical metadata column. (Max 1)}
#'   \item{\code{facet.by} - }{Categorical metadata column(s) only.}
#'   \item{\code{limit.by} - }{Any metadata column(s).}
#' }
#' 
#' All built-in color palettes are colorblind-friendly. The available 
#' categorical palette names are: \code{"okabe"}, \code{"carto"}, \code{"r4"}, 
#' \code{"polychrome"}, \code{"tol"}, \code{"bright"}, \code{"light"}, 
#' \code{"muted"}, \code{"vibrant"}, \code{"tableau"}, \code{"classic"}, 
#' \code{"alphabet"}, \code{"tableau20"}, \code{"kelly"}, and \code{"fishy"}.
#' 
#' Patterns are sourced from the magick R package. Pattern names are: 
#' \code{"bricks"}, \code{"hexagons"}, \code{"horizontalsaw"}, 
#' \code{"hs_fdiagonal"}, \code{"fishscales"}, \code{"verticalsaw"}, 
#' \code{"checkerboard"}, \code{"octagons"}, \code{"right45"}, 
#' \code{"hs_cross"}, \code{"hs_bdiagonal"}, \code{"hs_diagcross"}, 
#' \code{"hs_horizontal"}, \code{"hs_vertical"}, \code{"left45"}, 
#' \code{"leftshingle"}, \code{"rightshingle"}, \code{"verticalbricks"}, 
#' \code{"verticalleftshingle"}, and \code{"verticalrightshingle"}.
#' 
#' Shapes can be given as per base R - numbers 0 through 17 for various shapes,
#' or the decimal value of an ascii character, e.g. a-z = 65:90; A-Z = 97:122 to use 
#' letters instead of shapes on the plot. Character strings may used as well.
#' 
#' @export
#' @seealso \code{\link{stats_table}}
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     adiv_boxplot(biom, x = "Body Site", metric = "Shannon")
#'     adiv_boxplot(biom, x = "Sex", metric = c("OTUs", "Shannon"), layers="b", color.by="Body Site", scales="free")
#'     adiv_boxplot(biom, x = "Body Site", metric = "Simpson", layers="p", color.by="Sex", xlab.angle=30)
#'     

adiv_boxplot <- function (
    biom, x = NULL, metric = "Shannon", layers = "lsb",
    color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, limit.by = NULL, 
    flip = FALSE, stripe = flip, p.adj = "fdr", p.label = 0.05, ci = 95, outliers = NULL,
    xlab.angle = 'auto', ...) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(c(as.list(environment()), list(...)), eval)
  cache_file <- get_cache_file("adiv_boxplot", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  arg_str <- as.args(params, fun = adiv_boxplot, indent = 2)
  history <- paste0(collapse = "\n", c(
    attr(biom, 'history', exact = TRUE),
    sprintf("fig  <- adiv_boxplot(%s)", arg_str) ))
  remove(list = setdiff(ls(), c("params", "history", "cache_file")))
  
  
  #________________________________________________________
  # Sanity checks. x and *.by are checked by boxplot_build.
  #________________________________________________________
  params %<>% within({
    if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
    metric %<>% validate_arg(biom, 'metric', 'adiv', n = c(1,Inf))
  })
  
  
  #________________________________________________________
  # Use the generalized boxplot function to make the plot
  #________________________________________________________
  p <- boxplot_build(params, adiv_boxplot, adiv_boxplot_data, adiv_boxplot_layers)
  
  
  #________________________________________________________
  # Attach history of biom modifications and this call
  #________________________________________________________
  attr(p, 'history') <- history
  
  
  
  set_cache_value(cache_file, p)
  return (p)
}


#______________________________________________________________
# Convert biom object to a data.frame
#______________________________________________________________
adiv_boxplot_data <- function (params) {
  
  ggdata <- adiv_table(
    biom    = params[['biom']], 
    metrics = params[['metric']],
    md      = TRUE,
    long    = TRUE )
  
  
  
  # Always put .metric last in facet.by list.
  #________________________________________________________
  metric <- params[['metric']]
  if (length(metric) > 1) {
    facet.by <- setdiff(params[['facet.by']], ".metric")
    params[['facet.by']] <- c(facet.by, ".metric")
    ggdata[['.metric']] %<>% factor(levels = metric)
    remove("facet.by")
  }
  remove("metric")
  
  
  
  
  attr(ggdata, 'params') <- params
  attr(ggdata, 'xcol')   <- params[['x']]
  attr(ggdata, 'ycol')   <- attr(ggdata, 'response', exact = TRUE)
  
  return (ggdata)
}



#______________________________________________________________
# Make adiv-specific layer tweaks
#______________________________________________________________
adiv_boxplot_layers <- function (layers) {
  
  # y-axis title
  #________________________________________________________
  yraw <- attr(layers, 'params', exact = TRUE)[['metric']]
  
  if (length(yraw) > 1) {
    ylab <- structure(
      .Data   = "Diversity (\u03B1)", 
      display = '"Diversity (\\u03B1)"' )
    
    # setLayer("labs", scales = "free_y")
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
