#' Visualize diversity or abundance as a boxplot, dotplot, etc.
#' 
#' @name plot
#' 
#' @param x   A BIOM object, as returned from \link{read.biom}.
#' 
#' @param formula   Definition of what to plot on the x- and y- axes, given as 
#'        \code{y ~ x}. For example, \bold{Shannon ~ `Body Site`}. 
#'        X-axis options are:
#'        \describe{
#'            \item{Metadata}{
#'              The name of a column from \code{metadata(biom)}.
#'            }
#'            \item{Ordination Method}{
#'              \bold{PCoA}, \bold{tSNE}, or \bold{NMDS}.
#'            }
#'        }
#'        Y-axis options are:
#'        \describe{
#'            \item{Alpha Diversity Metrics (one or more)}{
#'              \bold{OTUs}, \bold{Shannon}, \bold{Chao1}, \bold{Simpson}, 
#'              and/or \bold{InvSimpson}.
#'            }
#'            \item{Beta Diversity Metrics (one only)}{
#'              \bold{Manhattan}, \bold{Euclidean}, \bold{Bray-Curtis}, 
#'              \bold{Jaccard}, or \bold{UniFrac}. \bold{Distance} will use 
#'              \bold{UniFrac} if a phylogenetic tree is present, or 
#'              \bold{Bray-Curtis} otherwise. Use in combination with the 
#'              \code{weighted} parameter. Metadata column names can be 
#'              prefixed with \bold{==} or \bold{!=} to limit distance 
#'              calculations to \emph{within} or \emph{between}, respectively, 
#'              those categories. See examples below.
#'            }
#'            \item{Taxa Abundances (one only)}{
#'              \bold{Kingdom}, \bold{Phylum}, \bold{Class}, \bold{Order}, 
#'              \bold{Family}, \bold{Genus}, \bold{Species}, \bold{Strain}, or 
#'              \bold{OTU}. Supported ranks will vary by biom. Run 
#'              \code{taxa.ranks(biom)} to see the available options. 
#'              Specifying \bold{Abundance} will default to the most precise 
#'              rank possible.
#'            }
#'           }
#'           
#' @param layers   What kind of plot to create. Options are \bold{box}, 
#'        \bold{violin}, \bold{dot}, \bold{strip}, \bold{crossbar}, 
#'        \bold{errorbar}, \bold{linerange}, and/or \bold{pointrange}. 
#'        Single letter abbreviations are also accepted. For instance,
#'        \code{c("box", "dot")} is equivalent to \code{c("b", "d")} and 
#'        \code{"bd"}. Default: "rls".
#'                 
#' @param color.by,pattern.by,shape.by,facet.by   Metadata column to color, 
#'        pattern, shape, and/or facet by. If that column is a \code{factor}, 
#'        the ordering of levels will be maintained in the plot.
#'        
#' @param colors,patterns,shapes   Names of the colors, patterns, and/or shapes
#'        to use in the plot. Available names can be found by running 
#'        \code{colors()}, \code{ggpattern::magick_pattern_names}, and 
#'        \code{0:25}, respectively. Use a named character vector to map them 
#'        to specific factor levels in the metadata.
#'        
#' @param p.min   Minimum adjusted p-value to display on the plot with a bracket.
#'        Set to \code{Inf} to display all p-values, or \code{-Inf} for no brackets.
#'        If a numeric vector with more than one value is provided, they will be
#'        used as breaks for asterisk notation. (Default: \code{0.05})
#'
#' @param p.adj   Method to use for multiple comparisons adjustment of p-values.
#'        Run \code{p.adjust.methods} for a list of available options.
#'        (Default: \code{fdr})
#'     
#' @param vline   How to calculate min/max of the \bold{crossbar}, 
#'        \bold{errorbar}, \bold{linerange}, and \bold{pointrange} layers. 
#'        Options are \bold{range}, \bold{ci} (confidence interval), \bold{sd}
#'        (standard deviation), \bold{se} (standard error), and \bold{mad}
#'        (median absolute deviation). You may optionally append a number to 
#'        \bold{ci} to specify the confidence level, for instance 
#'        \code{vline = "ci95"} will calculate the 95% confidence interval, 
#'        whereas \code{vline = "ci99"} will give the 99% confidence interval.
#'        The center mark of \code{crossbar} and \code{pointrange} represents
#'        the mean, except when \code{vline="mad"} in which case it represents
#'        the median.
#'        Default: \code{ci95}
#'        
#' @param xlab.angle   How to rotate the tick labels on the x-axis. 
#'        \bold{'auto'} (the default), automatically selects a rotation value. 
#'        \bold{0}, \bold{30}, and \bold{90} sets the angle to horizontal, 
#'        angled, and vertical, respectively.
#'        
#' @param rline   For plots with \code{Reads} as the \code{x} or \code{y} axis,
#'        highlight this rarefaction depth with a dashed line.
#'        
#' @param ...   Parameters passed on to ggplot2 functions. Prefixing a 
#'        parameter name with \bold{b.}, \bold{v.}, \bold{d.}, or \bold{s.} 
#'        forces that parameter to be passed to, and only to, 
#'        \bold{geom_boxplot}, \bold{geom_violin}, \bold{geom_dotplot}, or 
#'        \bold{geom_jitter}, respectively. Otherwise, parameters are passed 
#'        along by matching against formal arguments.
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
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile) %>% rarefy()
#'     
#'     plot(biom, Shannon ~ `Body Site`)
#'     plot(biom, Shannon ~ Sex, layers="vb", color.by="Body Site")
#'     
#'     plot(biom, Simpson ~ `Body Site`, layers="p", color.by="Sex", xlab.angle=30)
#'     
#'     # Ordination
#'     plot(biom, bray ~ nmds)
#'     
#'     # Dissimilarity boxplots
#'     plot(biom, UniFrac ~ `==Body Site`)
#'     
#'     # Dissimilarity Heatmap
#'     plot(biom, UniFrac ~ heatmap)
#'     
#'     # Taxa abundance boxplots
#'     plot(biom, Phylum ~ .)
#'     
#'     # Taxa stacked abundance
#'     plot(biom, Phylum ~ stacked)
#'     
#'     # Taxa abundance heatmap
#'     plot(biom, Phylum ~ heatmap)
#'     
#'
plot.BIOM <- function (x, formula, layers = "rls", color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, colors = NULL, patterns = NULL, shapes = NULL, p.min = 0.05, p.adj = "fdr", vline = "ci95", xlab.angle = 'auto', rline = NULL, ...) {
  
  biom <- x
  dots <- list(...)
  
  
  #--------------------------------------------------------------
  # pattern.by requires non-cran ggpattern package
  #--------------------------------------------------------------
  if (!is.null(pattern.by) && !nzchar(system.file(package = "ggpattern")))
      stop(simpleError(paste0(
        "\n",
        "Error: rbiom requires the R package 'ggpattern' to be\n",
        "installed in order to create plots with patterned fills.\n\n",
        "Please run the following commands to install 'ggpattern':\n",
        "   install.packages('remotes')\n",
        "   remotes::install_github('coolbutuseless/ggpattern')\n\n" )))
  
  
  #--------------------------------------------------------------
  # Sanity checks
  #--------------------------------------------------------------
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  
  
  #--------------------------------------------------------------
  # Extract the x and y variables from the formula/vector/list
  #--------------------------------------------------------------
  if (is(formula, "formula") && length(all.vars(formula, unique = FALSE)) == 2) {
    x <- all.vars(formula[[3]])
    y <- all.vars(formula[[2]])
    
  } else if ((is.character(formula) || is.list(formula)) && length(formula) == 2) {
    x <- formula[[1]]
    y <- formula[[2]]
    
  } else {
    stop("Please provide a valid formula.")
  }
  
  
  
  #--------------------------------------------------------------
  # Replace shortcut keywords
  #--------------------------------------------------------------
  if (identical(tolower(x), "abundance")) x <- tail(c('OTU', taxa.ranks(biom)), 1)
  if (identical(tolower(y), "abundance")) y <- tail(c('OTU', taxa.ranks(biom)), 1)
  if (identical(tolower(x), "distance"))  x <- ifelse(has.phylogeny(biom), 'unifrac', 'bray-curtis')
  if (identical(tolower(y), "distance"))  y <- ifelse(has.phylogeny(biom), 'unifrac', 'bray-curtis')
  if (identical(x, "."))                { x <- ".all"; biom$metadata[['.all']] <- ".all" }
  
  
  #-----------------------------------------------
  # Assign a 'mode' attribute to x and y.
  # Possible modes: taxa, adiv, bdiv, meta, ord,
  #                 Reads, Samples, .
  #-----------------------------------------------
  x    <- validate_metrics(biom, x)
  y    <- validate_metrics(biom, y)
  mode <- paste(attr(y, 'mode', exact = TRUE), "~", attr(x, 'mode', exact = TRUE))
  
  
  
  #==================================================================
  # Combination handling
  #------------------------------------------------------------------
  # Rarefied ~ Reads         = Reads Retained
  # Rarefied ~ Samples       = Samples Retained
  # Rarefied ~ Shannon       = rarefaction curve
  # Shannon  ~ `Body Site`   = boxplot*
  # UniFrac  ~ PCoA          = scatter
  # UniFrac  ~ `==Body Site` = between/within boxplots
  # UniFrac  ~ Samples       = heatmap (samples vs samples)
  # Reads    ~ Phylum        = boxplot
  # Reads    ~ Samples       = stacked barplot (color.by='Phylum')
  # Phylum   ~ Samples       = heatmap (taxa vs samples)
  # Phylum   ~ `Body Site`   = heatmap (taxa vs groups of samples)
  #
  # *should work without any metadata column
  #==================================================================
  
  fn <-  if (mode == "Rarefied ~ Reads")   { plot_rarefied   #o
  } else if (mode == "Rarefied ~ Samples") { plot_rarefied   #o
  } else if (mode == "Rarefied ~ adiv")    { plot_rarefied   #o
  } else if (mode == "adiv ~ .")           { plot_factor     #-
  } else if (mode == "adiv ~ factor")      { plot_factor     #-
  } else if (mode == "adiv ~ numeric")     { plot_numeric
  } else if (mode == "bdiv ~ ord")         { plot_ordination
  } else if (mode == "bdiv ~ factor")      { plot_factor     #-
  } else if (mode == "bdiv ~ clust")       { plot_heatmap    #>
  } else if (mode == "rank ~ .")           { plot_factor     #-
  } else if (mode == "rank ~ factor")      { plot_factor     #-
  } else if (mode == "rank ~ stacked")     { plot_stacked
  } else if (mode == "rank ~ clust")       { plot_heatmap    #>
  } else if (mode == "taxon ~ factor")     { plot_factor     #-
  } else if (mode == "taxon ~ numeric")    { plot_numeric
  } else { stop("Invalid formula of form '", mode, "'") }
  
  
  args <- as.list(match.call())
  args <- args[names(args)]
  args[['biom']] <- biom
  args[['x']]    <- x
  args[['y']]    <- y
  
  p <- do.call(fn, args)
  return (p)
  
}



#--------------------------------------------------------------
# Assign colors to categorical values.
#--------------------------------------------------------------
assign_colors <- function (vals, keys) {
  keys <- levels(keys)
  n    <- length(keys)
  
  if (is.null(vals)) {
    vals <- if (n <= 2) { # Colorblind friendly palette
      c('#00B9EB', '#ED5F16')
    } else if (n <= 8) {  # Colorblind friendly palette of 8 (jfly.iam.u-tokyo.ac.jp/color/)
      c("#0072B2", "#D55E00", "#CC79A7", "#009E73", "#F0E442", "#56B4E9", "#E69F00", "#999999")
    } else if (n <= 11) {   # Andrea's set of 11
      c('#66CDAA', '#FF8C69', '#8DB6CD', '#008B8B', '#FF6EB4', '#A2CD5A', '#FF6347', '#FFC125', 
        '#EEC591', '#BEBEBE', '#CD96CD')
    } else if (n <= 20) {   # Andrea's set of 20
      c('#8B4500', '#CD853F', '#FF1493', '#FFB5C5', '#CDC5BF', '#6C7B8B', '#B22222', '#FF0000', 
        '#FF7F24', '#FFD700', '#00FF00', '#2E8B57', '#00FFFF', '#63B8FF', '#0000FF', '#191970', 
        '#8B008B', '#A020F0', '#DA70D6', '#F08080')
    } else {                # Dan's set of 24
      c('#1C86EE', '#E31A1C', '#008B00', '#6A3D9A', '#A52A2A', '#FF7F00', '#FFD700', '#7EC0EE', 
        '#FB9A99', '#90EE90', '#CAB2D6', '#FDBF6F', '#B3B3B3', '#EEE685', '#B03060', '#FF83FA', 
        '#FF1493', '#0000FF', '#36648B', '#00CED1', '#00FF00', '#8B8B00', '#CDCD00', '#8B4500')
    }
  }
  
  assign_cleanup("colors", vals, keys)
}


#--------------------------------------------------------------
# Assign patterns to categorical values.
#--------------------------------------------------------------
assign_patterns <- function (vals, keys) {
  keys <- levels(keys)
  
  if (is.null(vals)) {
    vals <- unique(c(
      'bricks', 'fishscales', 'right45', 'horizonal_saw', 
      'hs_cross', 'crosshatch45', 'crosshatch',
      ggpattern::magick_pattern_names ))
  }
  
  assign_cleanup("patterns", vals, keys)
}


#--------------------------------------------------------------
# Assign shapes to categorical values.
#--------------------------------------------------------------
assign_shapes <- function (vals, keys) {
  keys <- levels(keys)
  n    <- length(keys)
  
  if (is.null(vals))
    vals <- if (n <= 4) 15:18 else 0:14
  
  assign_cleanup("shapes", vals, keys)
}


#--------------------------------------------------------------
# Ensure correct length/keys for colors, patterns, & shapes
#--------------------------------------------------------------
assign_cleanup <- function (mode, vals, keys) {
  n <- length(keys)
  
  if (is.null(names(vals))) {
    if (length(vals) < n) vals <- rep_len(vals, n)
    if (length(vals) > n) vals <- vals[seq_len(n)]
    vals <- setNames(vals, keys)
    
  } else {
    missing <- setdiff(keys, names(vals))
    if (length(missing) > 3) missing <- c(head(missing, 3), "...")
    if (length(missing) > 0)
      stop("Missing ", mode, " assignments for: ", paste(sep=", ", missing))
    vals <- vals[keys]
  }
  
  return (vals)
}

