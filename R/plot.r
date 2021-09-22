#' Visualize the BIOM data.
#' 
#' Provide two terms to \code{plot()} and it will automatically produce the
#' appropriate chart type. The terms that \code{plot} understands and valid 
#' combinations are listed in the Details section. Not all parameters used by
#' all chart types.
#' 
#' @md
#' @section Terms:
#' Terms come in the nine categories given below. The values are 
#' case-insensitive and can be unambiguously abbreviated.
#' \itemize{
#'   \item{\bold{Ordination Method (ORD): }}{
#'     \code{PCoA}, \code{tSNE}, or \code{NMDS}.
#'   }
#'   \item{\bold{Alpha Diversity Metric (ADIV): }}{
#'     \code{OTUs}, \code{Shannon}, \code{Chao1}, \code{Simpson}, 
#'     or \code{InvSimpson}.
#'   }
#'   \item{\bold{Beta Diversity Metric (BDIV): }}{
#'     \code{Manhattan}, \code{Euclidean}, \code{Bray-Curtis}, \code{Jaccard}, 
#'     or \code{UniFrac}.
#'   }
#'   \item{\bold{Distance Formula (DIST): }}{
#'     \code{correlation}, \code{euclidean}, \code{maximum}, \code{manhattan}, 
#'     \code{canberra}, \code{binary}, or \code{minkowski}.
#'   }
#'   \item{\bold{Taxonomic Rank (RANK): }}{
#'     \code{Kingdom}, \code{Phylum}, \code{Class}, \code{Order}, 
#'     \code{Family}, \code{Genus}, \code{Species}, \code{Strain}, or 
#'     \code{OTU}. Supported ranks will vary by BIOM. Run 
#'     \code{taxa.ranks(biom)} to see the available options.
#'   }
#'   \item{\bold{Individual Taxon (TAXON): }}{
#'     The name of a taxon from the BIOM's \code{taxonomy()}.
#'     For instance, \code{Firmicutes} or \code{Prevotella}.
#'   }
#'   \item{\bold{Metadata (FACTOR or NUMERIC): }}{
#'     The name of a column from the BIOM's \code{metadata()}.
#'     For instance, \code{`Body Site`} or \code{Age}.
#'   }
#'   \item{\bold{Clustering Method (CLUST): }}{
#'     \code{average}, \code{ward}, \code{mcquitty}, \code{single}, 
#'     \code{median}, \code{complete}, or \code{centroid}. The following 
#'     aliases are also understood: \code{heatmap} = \code{complete}, 
#'     \code{UPGMA} = \code{average}, \code{WPGMA} = \code{mcquitty}, 
#'     \code{WPGMC} = \code{median}, and \code{UPGMC} = \code{centroid}.
#'   }
#'   \item{\bold{Special: }}{
#'     \code{Rarefied}, \code{Reads}, \code{Samples}, \code{.}, \code{stacked}
#'   }
#' }
#' 
#' @section Term Combinations:
#' 
#' | Combination        | Chart Type  | Example                                 | 
#' | ------------------ | ----------- | --------------------------------------- | 
#' | ADIV ~ .           | Super box   | `` plot(biom, Shannon ~ .)           `` |
#' | ADIV ~ FACTOR      | Super box   | `` plot(biom, Shannon ~ `Body Site`) `` |
#' | BDIV ~ FACTOR      | Super box   | `` plot(biom, Bray ~ Sex)            `` |
#' | RANK ~ .           | Super box   | `` plot(biom, Phylum ~ .)            `` |
#' | TAXON ~ FACTOR     | Super box   | `` plot(biom, Firmicutes ~ Sex)      `` |
#' | ADIV ~ NUMERIC     | Line chart  | `` plot(biom, Shannon ~ Age)         `` |
#' | TAXON ~ NUMERIC    | Line chart  | `` plot(biom, Prevoltella ~ BMI)     `` |
#' | BDIV ~ CLUST       | Heatmap     | `` plot(biom, UniFrac ~ ward)        `` |
#' | RANK ~ CLUST       | Heatmap     | `` plot(biom, Genus ~ heatmap)       `` |
#' | BDIV ~ ORD         | Ordination  | `` plot(biom, Bray ~ NMDS)           `` |
#' | RANK ~ stacked     | Stacked bar | `` plot(biom, Family ~ stacked)      `` |
#' | Rarefied ~ Reads   | Rarefaction | `` plot(biom, Rarefied ~ Reads)      `` |
#' | Rarefied ~ Samples | Rarefaction | `` plot(biom, Rarefied ~ Samples)    `` |
#' | Rarefied ~ ADIV    | Rarefaction | `` plot(biom, Rarefied ~ Shannon)    `` |
#' 
#' 
#' 
#' @name plot
#' 
#' @param x   A BIOM object, as returned from \link{read.biom}.
#' 
#' @param formula   Combination of terms to plot, either in \code{y ~ x} or
#'        \code{c("x", "y")} form.
#'           
#' @param layers   See "Layers" section for details. Options for super box
#'        plots are \bold{box}, \bold{bar} (r), \bold{violin}, \bold{dot}, 
#'        \bold{strip}, \bold{crossbar}, \bold{errorbar}, \bold{linerange}, and
#'        \bold{pointrange}. Options for ordination plots are: \bold{point}, 
#'        \bold{centroid}, \bold{ellipse}, and \bold{name} for samples, 
#'        and \bold{mean}, \bold{taxon}, and \bold{arrow} for taxa biplots.
#'        Single letter abbreviations are also accepted. For instance,
#'        \code{c("box", "dot")} is equivalent to \code{c("b", "d")} and 
#'        \code{"bd"}. Default: "rls" for super box plots and "pce"/"p" for
#'        ordinations with/without a \code{color.by} argument.
#'                 
#' @param color.by,pattern.by,shape.by,label.by,sort.by,facet.by   Metadata 
#'        column to color, pattern, shape, label, sort, and/or facet by. If 
#'        that column is a \code{factor}, the ordering of levels will be 
#'        maintained in the plot.
#'        
#' @param colors,patterns,shapes,facets,xvals   Names of the colors, patterns,
#'        shapes, facets, and/or x values to use in the plot. Available values 
#'        for colors, patterns, and shapes are given by \code{colors()}, 
#'        \code{ggpattern::magick_pattern_names}, and \code{0:25},
#'        respectively. Use a named character vector to map them  to specific
#'        factor levels in the metadata. \code{facets} and \code{xvals} are
#'        coerced to unnamed character vectors. If the length of these vectors
#'        is less than the values present in their corresponding metadata 
#'        column, then the data set will be subseted accordingly.
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
#' @param se   How to calculate min/max of the \bold{crossbar}, 
#'        \bold{errorbar}, \bold{linerange}, and \bold{pointrange} layers. 
#'        Options are \bold{range}, \bold{ci} (confidence interval), \bold{sd}
#'        (standard deviation), \bold{se} (standard error), and \bold{mad}
#'        (median absolute deviation). You may optionally append a number to 
#'        \bold{ci} to specify the confidence level, for instance 
#'        \code{se = "ci95"} will calculate the 95% confidence interval, 
#'        whereas \code{se = "ci99"} will give the 99% confidence interval.
#'        The center mark of \code{crossbar} and \code{pointrange} represents
#'        the mean, except when \code{se="mad"} in which case it represents
#'        the median. In the case of trendlines, \code{se} must be in "ciXX"
#'        form. Set to \code{NULL} to disable. Default: \code{ci95}
#'        
#' @param rline   On rarefaction plots, highlight this rarefaction depth with a
#'        dashed line. (Default: NULL)
#'        
#' @param xlab.angle   How to rotate the tick labels on the x-axis. 
#'        \bold{'auto'} (the default), automatically selects a rotation value. 
#'        \bold{0}, \bold{30}, and \bold{90} sets the angle to horizontal, 
#'        angled, and vertical, respectively.
#'        
#' @param weighted   When employing a beta diversity metric, use the weighted
#'        version. (Default: \code{TRUE})
#'        
#' @param rank   What rank of taxa to display. E.g. "Phylum", "Genus", etc. Run
#'        \code{taxa.ranks()} to see all options for a given BIOM object. The
#'        default, \code{"auto"}, selects the lowest level.
#'        
#' @param taxa   Which taxa to display. An integer value will show the top n
#'        most abundant taxa. A value 0 <= n < 1 will show any taxa with at
#'        least that level of statistical significance. A character vector of
#'        taxon names will show only those taxa. (Default: \code{5} for 
#'        ordination biplots, \code{10} for heatmaps.)
#'        
#' @param abbr   When selecting taxa by name, allow abbreviated 'taxa' values, 
#'        e.g. \code{c('staph', 'lact')}. (Default: TRUE).
#'        
#' @param other   Should non-selected taxa be displayed as an "Other" group?
#'        (Default: \code{FALSE})
#'        
#' @param anno   Annotations to include on the plot. Options are: \bold{title},
#'        \bold{method}, \bold{p.value}, \bold{r.squared}, \bold{statistic},
#'        \bold{aic}, \bold{bic}, and \bold{axes}. See the 'stats' attribute of 
#'        the returned plot for the complete list, which changes based on the 
#'        test function used. Single letter abbreviations are also accepted. 
#'        (Default: \code{"tmpr"})
#'        
#' @param gradient   For heatmaps, the color gradient to use for the cells. 
#'        (Default: \code{heat.colors(20)}) 
#'        
#' @param normalize.rows   For heatmaps, should each row (taxon) have its 
#'        values rescaled from min-max to 0-1. (Default: \code{TRUE})
#'        
#' @param dist   For \code{RANK ~ CLUST}, the distance metric to use.
#'        (Default: \code{"euclidean"})
#'        
#' @param model   For regressions, the formula to use in the smoothing function.
#'        For example: \code{y ~ x}, \code{y ~ log(x)}, or \code{y ~ poly(x,2)}.
#'        (Default: \code{y ~ x})
#'        
#' @param regr   For regressions, the smoothing function to use.
#'        For example: \code{lm}, \code{glm}, \code{gam}, \code{loess},
#'        or a function with similar input parameters.
#'        (Default: \code{lm})
#'        
#' @param ...   Parameters for underlying functions.
#'        See "Additional Parameters" section for details. 
#'        For heatmaps, \code{...} is handled by \code{pheatmap::pheatmap()}.
#'        Otherwise, parameters are matched to formal arguments of ggplot2
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
#'     
#'     plot(biom, Shannon ~ `Body Site`)
#'     plot(biom, Shannon ~ Sex, layers="vb", color.by="Body Site")
#'     
#'     plot(biom, Simpson ~ `Body Site`, layers="p", color.by="Sex", xlab.angle=30)
#'     
#'     # Ordination
#'     plot(biom, bray ~ nmds, color.by="Body Site")
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
plot.BIOM <- function (
  x, formula, layers = "rls", 
  color.by = NULL, pattern.by = NULL, shape.by = NULL, label.by = NULL, 
  sort.by = NULL, facet.by = NULL, xvals = nuLL, colors = NULL, 
  patterns = NULL, shapes = NULL, facets = NULL, p.min = 0.05, p.adj = "fdr", 
  se = "ci95", rline = NULL, xlab.angle = 'auto', weighted = TRUE, 
  rank = "auto", taxa = NULL, abbr = TRUE, other = FALSE, anno = "tmprf", 
  gradient = heat.colors(20), normalize.rows = TRUE, dist = "euclidean", 
  model = y ~ x, regr = lm, ...) {
  
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
  if (identical(tolower(x), "abundance")) x <- tail(unique(c('OTU', taxa.ranks(biom))), 1)
  if (identical(tolower(y), "abundance")) y <- tail(unique(c('OTU', taxa.ranks(biom))), 1)
  if (identical(tolower(x), "distance"))  x <- ifelse(has.phylogeny(biom), 'unifrac', 'bray-curtis')
  if (identical(tolower(y), "distance"))  y <- ifelse(has.phylogeny(biom), 'unifrac', 'bray-curtis')
  if (identical(x, "."))                { x <- ".all"; biom$metadata[['.all']] <- ".all" }
  
  
  #-----------------------------------------------
  # Assign a 'mode' attribute to x and y.
  # Possible modes: taxa, adiv, bdiv, meta, ord,
  #                 Reads, Samples, .
  #-----------------------------------------------
  x     <- validate_metrics(biom, x)
  y     <- validate_metrics(biom, y)
  xmode <- attr(x, 'mode', exact = TRUE)
  ymode <- attr(y, 'mode', exact = TRUE)
  mode  <- paste(ymode, "~", ifelse(x == ".all", ".", xmode))
  remove("xmode", "ymode")
  
  
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
  } else if (mode == "rank ~ stacked")     { plot_stacked
  } else if (mode == "rank ~ clust")       { plot_heatmap    #>
  } else if (mode == "taxon ~ factor")     { plot_factor     #-
  } else if (mode == "taxon ~ numeric")    { plot_numeric
  } else { stop("Invalid formula of form '", mode, "'") }
  
  
  args <- as.list(match.call())
  args <- args[setdiff(names(args), 'formula')]
  args[['biom']] <- biom
  args[['x']]    <- x
  args[['y']]    <- y
  
  p <- do.call(fn, args)
  
  
  if (interactive() && !is.null(attr(p, 'err')))
    warning(attr(p, 'err'), call. = FALSE)
  
  
  return (p)
  
}
