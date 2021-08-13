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
#' | RANK ~ FACTOR      | Super box   | `` plot(biom, Phylum ~ `Body Site`)  `` |
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
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile) %>% rarefy()
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
  color.by = NULL, pattern.by = NULL, shape.by = NULL, label.by = NULL, sort.by = NULL, facet.by = NULL, 
  colors = NULL, patterns = NULL, shapes = NULL, p.min = 0.05, p.adj = "fdr", se = "ci95", rline = NULL, 
  xlab.angle = 'auto', weighted = TRUE, rank = "auto", taxa = NULL, abbr = TRUE, other = FALSE, anno = "tmprf",
  gradient = heat.colors(20), normalize.rows = TRUE, dist = "euclidean", model = y ~ x, method = lm, ...) {
  
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
  args <- args[setdiff(names(args), 'formula')]
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
  
  vals <- c(16, 17, 15, 3, 7, 8)
  if (n > 6)  vals <- c(0:14)
  if (n > 15) vals <- c(65:90, 97:122)
  if (n > 52) vals <- rep(vals, ceiling(n / 52))
  
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


#--------------------------------------------------------------
# Expand layer string / char vector to possible values
# x = "dr" or c("d", "bar") or more examples below
# choices = c(d = "dot", r = "bar")
# default = "dr"
#--------------------------------------------------------------
layer_match <- function (x, choices, default) {
  
  if (is.null(default))        default        <- choices[[1]]
  if (is.null(x))              x              <- default 
  # if (is.null(names(choices))) names(choices) <- unlist(substr(choices, 1, 1))
  
  x <- x[nchar(x) > 0]
  if (length(x) == 0) return (default)
  
  x    <- tolower(x)
  vals <- unname(choices)
  
  if (length(x) > 1) {
    x <- vals[pmatch(x, tolower(vals))] # c("d", "bar") => c("dot", "bar")
    
  } else if (nchar(x) == 1) {
    if (is.null(names(choices))) {
      x <- vals[pmatch(x, tolower(vals))]
    } else {
      x <- unname(choices[x])  # "r" => "bar"
    }
    
  } else if (is.na(pmatch(x, tolower(vals)))) {
    x <- strsplit(x, '')[[1]]
    
    if (is.null(names(choices))) {
      x <- vals[pmatch(x, tolower(vals))]
    } else {
      x <- unname(choices[x]) # "dr" => c("dot", "bar")
    }
    
  } else {
    x <- vals[pmatch(x, tolower(vals))] # "do" => "dot"
  }
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) return (default)
  
  return (x)
}


#--------------------------------------------------------------
# Create the plot and add each layer with its arguments.
# Also attach a human-readable version of the plot command.
#--------------------------------------------------------------
layers_toPlot <- function (layers, dots) {
  
  stopifnot('ggplot' %in% names(layers))
  
  #--------------------------------------------------------------
  # 'layers' can draw on the layers listed here
  #--------------------------------------------------------------
  layerspecs <- list(
    'ggplot'     = list('.fn' = "ggplot",                    '.regex' = "^g(|gplot)\\."),
    'point'      = list('.fn' = "geom_point",                '.regex' = "^p(|oint)\\."),
    'bar'        = list('.fn' = "geom_bar",                  '.regex' = "^(|ba)r\\."),
    'violin'     = list('.fn' = "geom_violin",               '.regex' = "^v(|iolin)\\."),
    'strip'      = list('.fn' = "geom_quasirandom",          '.regex' = "^s(|trip)\\."),
    'box'        = list('.fn' = "geom_boxplot",              '.regex' = "^b(|ox)\\."),
    'dot'        = list('.fn' = "geom_beeswarm",             '.regex' = "^d(|ot)\\."),
    'errorbar'   = list('.fn' = "geom_errorbar",             '.regex' = "^e(|rrorbar)\\."),
    'crossbar'   = list('.fn' = "geom_crossbar",             '.regex' = "^c(|rossbar)\\."),
    'linerange'  = list('.fn' = "geom_linerange",            '.regex' = "^l(|inerange)\\."),
    'pointrange' = list('.fn' = "geom_pointrange",           '.regex' = "^p(|ointrange)\\."),
    'regression' = list('.fn' = "stat_smooth",               '.regex' = "^r(|egression)\\."),
    'stats_text' = list('.fn' = "geom_text",                 '.regex' = "^a(|nno)\\."),
    'brackets'   = list('.fn' = "geom_segment",              '.regex' = "^h(|line)\\."),
    'fill'       = list('.fn' = "scale_fill_manual",         '.regex' = "^fill\\."),
    'color'      = list('.fn' = "scale_color_manual",        '.regex' = "^color\\."),
    'shape'      = list('.fn' = "scale_shape_manual",        '.regex' = "^shape\\."),
    'pattern'    = list('.fn' = "scale_pattern_type_manual", '.regex' = "^pattern\\."),
    'facet'      = list('.fn' = "facet_wrap",                '.regex' = "^f(|acet)\\."),
    'labs'       = list('.fn' = "labs",                      '.regex' = "^labs\\."),
    'x_disc'     = list('.fn' = "scale_x_discrete",          '.regex' = "^x(|axis)\\."),
    'x_cont'     = list('.fn' = "scale_x_continuous",        '.regex' = "^x(|axis)\\."),
    'y_cont'     = list('.fn' = "scale_y_continuous",        '.regex' = "^y(|axis)\\."),
    'theme_bw'   = list('.fn' = "theme_bw",                  '.regex' = "^(theme_|)bw\\."),
    'theme'      = list('.fn' = "theme",                     '.regex' = "^t(|heme)\\.") )
  
  
  #--------------------------------------------------------------
  # Other common settings
  #--------------------------------------------------------------
  layerspecs[['stats_text']][['show.legend']] <- FALSE
  layerspecs[['errorbar']][['width']] <- 0.5
  layerspecs[['crossbar']][['width']] <- 0.5
  
  layerspecs[['stats_text']] %<>% c(list(
    'color'   = "black",
    'mapping' = list(x=".x", y=".y", label=".label") ))
  layerspecs[['brackets']] %<>% c(list(
    'color'   = "black",
    'mapping' = list(x=".x", y=".y", xend=".xend", yend=".yend") ))
  
  for (i in c('errorbar', 'crossbar', 'linerange', 'pointrange')) {
    layerspecs[[i]][['mapping']] <- list(y=".y", ymin=".ymin", ymax=".ymax")
  }
  
  
  #--------------------------------------------------------------
  # Only subset the data when multiple data layers are present
  #--------------------------------------------------------------
  if (isTRUE(length(unique(layers[['ggplot']][['data']][['.src']])) > 1)) {
    layerspecs[['stats_text']][['data']] <- ~ subset(., .src == "stats")
    layerspecs[['brackets']][['data']]   <- ~ subset(., .src == "brackets")
    
    for (i in c('point', 'bar', 'violin', 'strip', 'box', 'dot', 'regression')) {
      layerspecs[[i]][['data']] <- ~ subset(., .src == "points")
    }
    
    for (i in c('errorbar', 'crossbar', 'linerange', 'pointrange')) {
      layerspecs[[i]][['data']] <- ~ subset(., .src == "vline")
    }
  }
  
  
  #--------------------------------------------------------------
  # Special cases for ggbeeswarm and ggpattern functions
  #--------------------------------------------------------------
  patterned <- isTRUE("pattern" %in% names(layers))
  for (layer in names(layers)) {
    
    if (!layer %in% names(layerspecs)) {
      stopifnot(isTRUE(nzchar(layers[[layer]][['.fn']])))
      stopifnot(is.function(layers[[layer]][['.fun']]))
    }
    
    if (patterned && layer %in% c("bar", "box", "crossbar", "violin")) {
      layerspecs[[layer]][['.fn']] %<>% paste0(., "_pattern")
    }
    
    if (patterned && layer %in% c("fill", "color", "shape")) {
      layerspecs[[layer]][['.fn']] <- paste0("scale_pattern_", layer, "_manual")
    }
    
    pkg <- "ggplot2"
    if (grepl("_beeswarm",    layerspecs[[layer]][['.fn']], fixed = TRUE)) pkg <- "ggbeeswarm"
    if (grepl("_quasirandom", layerspecs[[layer]][['.fn']], fixed = TRUE)) pkg <- "ggbeeswarm"
    if (grepl("_pattern",     layerspecs[[layer]][['.fn']], fixed = TRUE)) pkg <- "ggpattern"
    
    layerspecs[[layer]][['.fun']] <- do.call(`::`, list(pkg, layerspecs[[layer]][['.fn']]))
    if (pkg != "ggplot2") layerspecs[[layer]][['.fn']] %<>% paste0(pkg, "::", .)
    
    # Merge custom and default top-level parameters
    for (i in setdiff(names(layerspecs[[layer]]), names(layers[[layer]])))
      layers[[layer]][[i]] <- layerspecs[[layer]][[i]]
    
    # Merge custom and default mapping parameters
    for (i in setdiff(names(layerspecs[[layer]][['mapping']]), names(layers[[layer]][['mapping']])))
      layers[[layer]][['mapping']][[i]] <- layerspecs[[layer]][['mapping']][[i]]
  }
  
  
  
  #--------------------------------------------------------------
  # Standardize the list order: ggplot() first, theme() last, etc
  #--------------------------------------------------------------
  layers <- layers[c(
    intersect(names(layerspecs), names(layers)),
    setdiff(names(layers), names(layerspecs))
  )]
  
  
  p    <- NULL
  cmds <- c()
  
  for (layer in names(layers)) {
    
    args <- layers[[layer]]
    fn   <- args[['.fn']]
    fun  <- args[['.fun']]
    regx <- args[['.regex']]
    args <- args[grep("^\\.", names(args), invert = TRUE)]
    
    
    # Unprefixed dot arguments, e.g. 'scales'="free_x"
    #--------------------------------------------------------------
    for (i in intersect(names(dots), formalArgs(fun)))
      args[[i]] <- dots[[i]]
    
    
    # Prefixed dot arguments, e.g. 'facet.scales'="free_x"
    #--------------------------------------------------------------
    for (i in grep(regx, names(dots), value = TRUE))
      args[[sub(regx, "", i, perl = TRUE)]] <- dots[[i]]
    
    
    # Don't specify an argument if it's already the default
    #--------------------------------------------------------------
    defaults <- formals(fun)
    for (i in intersect(names(defaults), names(args)))
      if (identical(defaults[[i]], args[[i]]))
        args[[i]] <- NULL
    
    
    # Create the aes object for mapping=
    #--------------------------------------------------------------
    if ('mapping' %in% names(args))
      args[['mapping']] <- do.call(
        what  = aes_string, 
        args  = args[['mapping']], 
        quote = TRUE )
    
    
    # Rewrite layer with the updated args
    #--------------------------------------------------------------
    layers[[layer]] <- args
    
    
    # Skip theme() unless it has arguments
    #--------------------------------------------------------------
    if (fn == "theme" && length(args) == 0) next
    
    
    # Show ggplot() layer as "ggplot(data)", rest more verbosely
    #--------------------------------------------------------------
    if (fn == "ggplot") {
      p <- do.call(fun, args) 
      cmds <- sprintf("%s(%s)", fn, as.args(args, fun=fun))
    } else {
      p    <- p + do.call(fun, args)
      cmds <- c(cmds, sprintf("%s(%s)", fn, as.args(args, indent = 4, fun=fun)))
    }
    
  }
  
  #--------------------------------------------------------------
  # Attach the number of facet rows and cols as plot attributes
  #--------------------------------------------------------------
  if ('nfacets' %in% names(attributes(dots))) {
    rc <- ggplot2::wrap_dims(
      n    = attr(dots, 'nfacets'), 
      nrow = layers[['facet']][['nrow']], 
      ncol = layers[['facet']][['ncol']] )
    attr(dots, 'facet.nrow') <- min(rc[[1]], attr(dots, 'nfacets'))
    attr(dots, 'facet.ncol') <- min(rc[[2]], attr(dots, 'nfacets'))
  }
  attr(p, 'facet.nrow') <- attr(dots, 'facet.nrow') %||% 1
  attr(p, 'facet.ncol') <- attr(dots, 'facet.ncol') %||% 1
  
  
  attr(p, 'cmd')  <- paste(collapse=" +\n  ", cmds)
  attr(p, 'data') <- layers[['ggplot']][['data']]
  
  return (p)
}


#------------------------------------------------------------------
# Convert a list of arguments to character strings.
# When indent > 0, produces a multi-line string.
# When fun is a function, puts the arguments in the expected order.
#------------------------------------------------------------------
as.args <- function (args = list(), indent = 0, fun = NULL) {
  
  stopifnot(is.list(args))
  stopifnot(is.numeric(indent))
  
  # Right-pad parameter names.
  fmt <- "%s = %s"
  if (isTRUE(indent > 0 && length(args) > 1))
    fmt <- paste(collapse="", rep(" ", indent)) %>%
      paste0("\n", ., "%-", max(nchar(names(args))), "s = %s")
  
  # Re-arrange parameter order; omit names where implicitly known.
  if (is.function(fun)) {
    f_args <- formalArgs(fun)
    args   <- args[c(intersect(f_args, names(args)), sort(setdiff(names(args), f_args)))]
    if (isTRUE(indent == 0))
      for (i in seq_along(args))
        if (names(args)[[i]] == f_args[[i]]) names(args)[i] <- "" else break
  }
  
  # Convert arguments to `eval`-able string representations
  strs <- c()
  for (i in seq_along(args)) {
    
    key <- names(args)[[i]]
    val <- args[[i]]
    
    val <- if (is.null(val))                    { "NULL"
    } else if (!is.null(attr(val, 'display')))  { attr(val, 'display') 
    } else if (is.character(val))               { glue::double_quote(val) 
    } else if (is.logical(val))                 { as.character(val)
    } else if (is.numeric(val))                 { as.character(val)
    } else if (is(val, 'BIOM'))                 { "biom"
    } else if (is.data.frame(val))              { "data"
    } else if (is(val, 'formula'))              { capture.output(val)[[1]]
    } else if (is.function(val))                { fun_toString(val)
    } else if (is(val, 'uneval'))               { aes_toString(val)
    } else if (is.factor(val))                  { as.character(val)
    } else                                      { capture.output(val) }
    
    if (length(val) > 1)
      val <- paste0("c(", paste(collapse=", ", val), ")")
    
    if (nzchar(key)) {
      key <- capture.output(as.name(key))
      val <- sprintf(fmt = fmt, key, val)
    }
    
    strs <- c(strs, val)
  }
  
  strs <- paste(strs, collapse = ", ")
  
  if (isTRUE(indent > 0 && length(args) > 1))
    strs <- paste0(strs, " ")
  
  return (strs)
}


#------------------------------------------------------------------
# Find a function's name
#------------------------------------------------------------------
fun_toString <- function (x) {
  
  pkg <- environment(x)[['.packageName']]
  if (!is.null(pkg))
    for (fn in getNamespaceExports(pkg))
      if (identical(x, do.call(`::`, list(pkg, fn)))) {
        if (!pkg %in% getOption("defaultPackages"))
          fn <- sprintf("%s::%s", pkg, fn)
        return (fn)
      }
  
  chr <- capture.output(x)
  if (nchar(chr) < 50) return (chr)
  
  return ("_Custom_Function")
}


#------------------------------------------------------------------
# Convert an aes object to a string
#------------------------------------------------------------------
aes_toString <- function (x) {
  
  # Consistently order the aes parameters
  keys <- names(x)
  sort <- c(
    "x", "xend", "xmin", "xmax", 
    "y", "yend", "ymin", "ymax", 
    "colour", "fill", "shape", "group", 
    "pattern_type", "pattern_fill", "label" )
  keys <- c(intersect(sort, keys), setdiff(keys, sort))
  
  results <- c()
  for (key in keys) {
    
    val <- x[[key]]
    val <- if (is(val, 'formula')) { capture.output(as.name(all.vars(val)))
    } else                         { glue::double_quote(val) }
    
    key %<>% sub(pattern = "colour", replacement = "color")
    key <- capture.output(as.name(key))
    
    results %<>% c(sprintf("%s = %s", key, val))
  }
  
  return (sprintf("aes(%s)", paste(collapse = ", ", results)))
}


#------------------------------------------------------------------
# rbind(), but add/rearrange columns as needed
#------------------------------------------------------------------
append_df <- function (x, y) {
  xy <- unique(c(names(x), names(y)))
  for (i in setdiff(xy, names(x))) x[[i]] <- NA
  for (i in setdiff(xy, names(y))) y[[i]] <- NA
  rbind(x[,xy,drop=F], y[,xy,drop=F])
}


#------------------------------------------------------------------
# Explicitly define the code to be displayed in cmd
#------------------------------------------------------------------
as.cmd <- function (expr, env=NULL) {
  if (is.null(env)) {
    cmd <- capture.output(substitute(expr))
  } else {
    cmd <- do.call(substitute, list(expr=substitute(expr), env=env))
    cmd <- capture.output(cmd)
  }
  
  if (length(cmd) > 1)
    cmd <- paste(trimws(cmd), collapse = " ")
  
  structure(expr, 'display' = cmd)
}


#------------------------------------------------------------------
# Identify and remove rows of df with bad values
#------------------------------------------------------------------
finite_check <- function (df, col=".y", metric=NULL) {
  
  if (all(is.finite(df[[col]])))
    return (NULL)
    
  bad <- which(!is.finite(df[[col]]))
  n   <- length(bad)
  
  if (is.null(metric))
    if (".metric" %in% names(df))
      if (length(unique(df[['.metric']])) == 1)
        metric <- df[1,'.metric']
  
  metric <- ifelse(is.null(metric), "", paste0(metric, " "))
  msg    <- ifelse(n == 1, 
    paste0("One sample had a non-finite ", metric, "value and was excluded from this plot:\n"),
    paste0(n, " samples had non-finite ", metric, "values and were excluded from this plot:\n") )
  
  
  if ('.sample' %in% names(df)) {
    msg %<>% paste0(
      glue::glue_collapse(
        x = paste(df[bad,'.sample'], "=", df[bad,col]), 
        width = 100, sep = ", ", last = ", and " ))
  } else {
    msg %<>% paste0(
      glue::glue_collapse(
        x = as.character(unique(df[bad,col])), 
        width = 100, sep = ", ", last = ", and " ))
  }
  
  list('bad' = bad, 'msg' = msg)
}


