#' Visualize BIOM data.
#' 
#' Provide two terms to \code{plot()} and it will automatically produce the
#' appropriate chart type. The terms that \code{plot} understands and valid 
#' combinations are listed in the Details section. Not all parameters are used 
#' by all chart types.
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
#'     and \code{InvSimpson}. Multiple allowed.
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
#'     \code{Family}, \code{Genus}, \code{Species}, \code{Strain}, and 
#'     \code{OTU}. Multiple allowed. Supported ranks will vary by BIOM. Run 
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
#' | Combination        | Chart Type   | Example                                  | 
#' | ------------------ | ------------ | ---------------------------------------- | 
#' | ADIV ~ .           | Box Plot     | `` plot(hmp50, Shannon ~ .)           `` |
#' | ADIV ~ FACTOR      | Box Plot     | `` plot(hmp50, Shannon ~ `Body Site`) `` |
#' | BDIV ~ FACTOR      | Box Plot     | `` plot(hmp50, Bray ~ Sex)            `` |
#' | RANK ~ .           | Box Plot     | `` plot(hmp50, Phylum ~ .)            `` |
#' | TAXON ~ FACTOR     | Box Plot     | `` plot(hmp50, Firmicutes ~ Sex)      `` |
#' | ADIV ~ NUMERIC     | Scatter Plot | `` plot(hmp50, Shannon ~ Age)         `` |
#' | ADIV ~ ADIV        | Scatter Plot | `` plot(hmp50, Depth ~ OTUs)          `` |
#' | TAXON ~ NUMERIC    | Scatter Plot | `` plot(hmp50, Prevoltella ~ BMI)     `` |
#' | BDIV ~ CLUST       | Heatmap      | `` plot(hmp50, UniFrac ~ ward)        `` |
#' | RANK ~ CLUST       | Heatmap      | `` plot(hmp50, Genus ~ heatmap)       `` |
#' | BDIV ~ ORD         | Ordination   | `` plot(hmp50, Bray ~ NMDS)           `` |
#' | RANK ~ stacked     | Stacked bar  | `` plot(hmp50, Family ~ stacked)      `` |
#' | Rarefied ~ Reads   | Rarefaction  | `` plot(hmp50, Rarefied ~ Reads)      `` |
#' | Rarefied ~ Samples | Rarefaction  | `` plot(hmp50, Rarefied ~ Samples)    `` |
#' | Rarefied ~ ADIV    | Rarefaction  | `` plot(hmp50, Rarefied ~ Shannon)    `` |
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
#'        \bold{spider}, \bold{ellipse}, and \bold{name} for samples, 
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
#'        \code{gridpattern::names_magick}, and \code{0:25},
#'        respectively. Use a named character vector to map them  to specific
#'        factor levels in the metadata. \code{facets} and \code{xvals} are
#'        coerced to unnamed character vectors. If the length of these vectors
#'        is less than the values present in their corresponding metadata 
#'        column, then the data set will be subseted accordingly.
#'        
#' @param rank   What rank of taxa to display. E.g. "Phylum", "Genus", etc. Run
#'        \code{taxa.ranks()} to see all options for a given BIOM object. The
#'        default, \code{NULL}, selects the lowest level.
#'        
#' @param taxa   Which taxa to display. An integer value will show the top n
#'        most abundant taxa. A value 0 <= n < 1 will show any taxa with that 
#'        mean abundance or greater (e.g. 0.1). A character vector of
#'        taxon names will show only those taxa. (Default: \code{10} for 
#'        heatmaps, \code{5} otherwise.)
#'
#' @param p.top   For \code{RANK ~ .} (taxa abundance) plots or 
#'        \code{BDIV ~ ORD} ordination biplots, only display taxa with the most 
#'        significant differences in abundance. If \code{p.top} is >= 1, then 
#'        the \code{p.top} most significant taxa are displayed. If \code{p.top} 
#'        is less than one, all taxa with an adjusted p-value <= \code{p.top} 
#'        are displayed. Recommended to be used in combination with the 
#'        \code{taxa} parameter to set a lower bound on the mean abundance of 
#'        considered taxa. (Default: \code{Inf})
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
#' @param perms   Number of random permutations to use for estimating statistical
#'        significance in \code{BDIV ~ ORD} ordinations. (Default: \code{1000})
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
#' @param abbr   When selecting taxa by name, allow abbreviated 'taxa' values, 
#'        e.g. \code{c('staph', 'lact')}. (Default: TRUE).
#'        
#' @param other   Should non-selected taxa be displayed as an "Other" group?
#'        (Default: \code{FALSE})
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
#' @param regr   For \code{BDIV ~ ORD}, the number of random permutations to
#'        use when calculating adonis statistics. (Default: \code{999})
#'        
#' @param safe   If \code{FALSE}, data.frame column names such as 
#'        \code{".metric"} will be auto-converted to \code{"Metric"} to improve
#'        human-readability. Conversion if aborted if a conflict is found with
#'        a metadata column name. (Default: \code{FALSE})
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
#'     plot(biom, c(OTUs, Shannon) ~ Sex, layers="vb", color.by="Body Site")
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
#'     plot(biom, c(Phylum, Genus) ~ .)
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
  sort.by = NULL, facet.by = NULL, xvals = NULL, colors = NULL,
  patterns = NULL, shapes = NULL, facets = NULL, 
  rank = NULL, taxa = NULL, 
  p.top = Inf, p.adj = "fdr", p.label = 0.05, perms = 1000,
  ci = 95, rline = NULL, xlab.angle = 'auto', weighted = TRUE,
  abbr = TRUE, other = FALSE,
  gradient = heat.colors(20), normalize.rows = TRUE, dist = "euclidean",
  model = y ~ x, regr = lm, safe = FALSE, ...) {
  
  
  #--------------------------------------------------------------
  # Reconstruct the command for provenance tracking
  #--------------------------------------------------------------
  all_args <- intersect(names(formals(plot.BIOM)), names(match.call()))
  all_args <- sapply(all_args, function (i) get(i))
  all_args <- c(all_args, list(...))
  history <- sprintf("plot(%s)", as.args(all_args, fun = plot.BIOM))
  remove("all_args")
  
  
  #--------------------------------------------------------------
  # Sanity checks
  #--------------------------------------------------------------
  biom <- x
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  
  
  params <- list(
    xval.by = as.vector(x), color.by = color.by, pattern.by = pattern.by, 
    shape.by = shape.by, label.by = label.by, sort.by = sort.by, facet.by = facet.by, 
    xvals = xvals, colors = colors, patterns = patterns, shapes = shapes, facets = facets, 
    rank = rank, taxa = taxa, 
    p.top = p.top, p.adj = p.adj, p.label = p.label, perms = perms,
    ci = ci, rline = rline, xlab.angle = xlab.angle, weighted = weighted,
    abbr = abbr, other = other, 
    gradient = gradient, normalize.rows = normalize.rows, dist = dist,
    model = model, regr = regr, safe = safe, ... )
  
  
  #--------------------------------------------------------------
  # Process the x/y values
  #--------------------------------------------------------------
  parsed <- parse_formula(biom, formula, .x = params[['.x']], .y = params[['.y']])
  x      <- parsed[['x']]
  y      <- parsed[['y']]
  mode   <- parsed[['mode']]
  if (identical(x, ".")) { x <- ".all"; biom$metadata[['.all']] <- ".all" }
  
  
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
  
  fun <- if (mode == "Rarefied ~ Reads")   { plot_rarefied   #o
  } else if (mode == "Rarefied ~ Samples") { plot_rarefied   #o
  } else if (mode == "Rarefied ~ adiv")    { plot_rarefied   #o
  } else if (mode == "adiv ~ .")           { boxplot     #-
  } else if (mode == "adiv ~ factor")      { boxplot     #-
  } else if (mode == "adiv ~ numeric")     { plot_numeric
  } else if (mode == "adiv ~ adiv")        { plot_numeric
  } else if (mode == "bdiv ~ ord")         { ordination_plot
  } else if (mode == "bdiv ~ factor")      { boxplot     #-
  } else if (mode == "bdiv ~ clust")       { plot_heatmap    #>
  } else if (mode == "rank ~ .")           { boxplot     #-
  } else if (mode == "rank ~ stacked")     { plot_stacked
  } else if (mode == "rank ~ clust")       { plot_heatmap    #>
  } else if (mode == "taxon ~ factor")     { boxplot     #-
  } else if (mode == "taxon ~ numeric")    { plot_numeric
  } else { stop("Invalid formula of form '", mode, "'") }
  
  p <- fun(biom = biom, x = x, y = y, layers = layers, mode = mode, params = params)
  
  
  if (interactive() && !is.null(attr(p, 'err', exact = TRUE)))
    warning(attr(p, 'err', exact = TRUE), call. = FALSE)
  
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  return (p)
  
}
