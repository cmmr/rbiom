
# default ====
#' documentation_default
#' 
#' @name documentation_default
#' @keywords internal
#' 
#' @param biom   An \code{rbiom}-class object, or data coercible with 
#'        [as_rbiom()].
#'        
#' @param tree  A \code{phylo} object representing the phylogenetic
#'        relationships of the taxa in \code{biom}. Only required when 
#'        computing UniFrac distances. Default: \code{otu_tree(biom)}
#'     
#' @param md  A character vector naming the metadata fields to include in the 
#'        output data frame, or \code{'.all'} to include all metadata fields.
#'        Default: \code{'.all'}
#' 
#' @param adiv   Alpha diversity metric(s) to use. Options are: \code{"OTUs"}, 
#'        \code{"Shannon"}, \code{"Chao1"}, \code{"Simpson"}, and/or 
#'        \code{"InvSimpson"}. Set \code{adiv=".all"} to use all metrics.
#'        Default: \code{"Shannon"} \cr\cr
#'        Multiple values allowed. Non-ambiguous abbreviations are allowed.
#' 
#' @param bdiv  Beta diversity distance algorithm(s) to use. Options are:
#'        \code{"Bray-Curtis"}, \code{"Manhattan"}, \code{"Euclidean"}, 
#'        \code{"Jaccard"}, and \code{"UniFrac"}. For \code{"UniFrac"}, a 
#'        phylogenetic tree must be present in \code{biom} or explicitly 
#'        provided via \code{tree=}. Default: \code{"Bray-Curtis"} \cr\cr
#'        Multiple/abbreviated values allowed.
#' 
#' @param taxa   Which taxa to display. An integer value will show the top n
#'        most abundant taxa. A value 0 <= n < 1 will show any taxa with that 
#'        mean abundance or greater (e.g. \code{0.1} implies >= 10%). A 
#'        character vector of taxa names will show only those named taxa. 
#'        Default: \code{6}.
#'        
#' @param ord    Method for reducing dimensionality. Options are:
#'        \itemize{
#'            \item{\code{"UMAP"} - }{ Uniform manifold approximation and projection; [uwot::umap()]. }
#'            \item{\code{"PCoA"} - }{ Principal coordinate analysis; [ape::pcoa()]. }
#'            \item{\code{"NMDS"} - }{ Nonmetric multidimensional scaling; [vegan::metaMDS()]. }
#'            \item{\code{"tSNE"} - }{ t-distributed stochastic neighbor embedding; [tsne::tsne()]. }
#'        }
#'        Default: \code{"UMAP"} \cr\cr
#'        Multiple values allowed. Non-ambiguous abbreviations are allowed.
#'     
#' @param weighted  Take relative abundances into account. When 
#'        \code{weighted=FALSE}, only presence/absence is considered.
#'        Default: \code{TRUE} \cr\cr
#'        Multiple values allowed.
#'        
#' @param rank   What rank(s) of taxa to display. E.g. \code{"Phylum"}, 
#'        \code{"Genus"}, \code{".otu"}, etc. An integer vector can also be 
#'        given, where \code{1} is the highest rank, \code{2} is the second 
#'        highest, \code{-1} is the lowest rank, \code{-2} is the second 
#'        lowest, and \code{0} is the OTU "rank". Run \code{taxa_ranks()} to 
#'        see all options for a given \code{rbiom} object. Default: \code{-1}.
#'        
#' @param lineage  Include all ranks in the name of the taxa. For instance,
#'        setting to \code{TRUE} will produce 
#'        \code{Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales}. 
#'        Otherwise the taxa name will simply be \code{Coriobacteriales}. You 
#'        want to set this to TRUE when \code{unc = "asis"} and you have taxa 
#'        names (such as \emph{Incertae_Sedis}) that map to multiple higher 
#'        level ranks. Default: \code{FALSE}
#'        
#' @param unc  How to handle unclassified, uncultured, and similarly ambiguous
#'        taxa names. Options are: 
#'        \itemize{
#'          \item{\code{"singly"} - }{ Replaces them with the OTU name. }
#'          \item{\code{"grouped"} - }{ Replaces them with a higher rank's name. }
#'          \item{\code{"drop"} - }{ Excludes them from the result. }
#'          \item{\code{"asis"} - }{ To not check/modify any taxa names. }
#'        }
#'        Default: \code{"singly"} \cr\cr
#'        Abbreviations are allowed.
#'        
#' @param other  Sum all non-itemized taxa into an "Other" taxa. When 
#'        \code{FALSE}, only returns taxa matched by the \code{taxa} 
#'        argument. Specifying \code{TRUE} adds "Other" to the returned set.
#'        A string can also be given to imply \code{TRUE}, but with that
#'        value as the name to use instead of "Other".
#'        Default: \code{FALSE}
#'        
#' @param sparse  If true, returns a sparse matrix as described by 
#'        [slam::simple_triplet_matrix()], otherwise returns a normal R
#'        matrix object. Default: \code{FALSE}
#'        
#' @param p.top   Only display taxa with the most significant differences in 
#'        abundance. If \code{p.top} is >= 1, then the \code{p.top} most 
#'        significant taxa are displayed. If \code{p.top} is less than one, all 
#'        taxa with an adjusted p-value <= \code{p.top} are displayed. 
#'        Recommended to be used in combination with the \code{taxa} parameter 
#'        to set a lower bound on the mean abundance of considered taxa. 
#'        Default: \code{Inf}
#'        
#' @param y.trans   The transformation to apply to the y-axis. Visualizing 
#'        differences of both high- and low-abundance taxa is best done with
#'        a non-linear axis. Options are: 
#'        \itemize{
#'          \item{\code{"sqrt"} - }{ square-root transformation }
#'          \item{\code{"log1p"} - }{ log(y + 1) transformation }
#'          \item{\code{NULL} - }{ no transformation }
#'        }
#'        These methods allow visualization of both high- and low-abundance
#'        taxa simultaneously, without complaint about 'zero' count
#'        observations. Default: \code{"sqrt"}
#'        
#' @param flip   Transpose the axes, so that taxa are present as rows instead
#'        of columns. Default: \code{FALSE}
#'
#' @param stripe   Shade every other x position. Default: \emph{same as flip}
#'
#' @param p.label   Minimum adjusted p-value to display on the plot with a 
#'        bracket.
#'        \itemize{
#'          \item{\code{p.label = 0.05} - }{ Show p-values that are <= 0.05. }
#'          \item{\code{p.label = 0} - }{ Don't show any p-values on the plot. }
#'          \item{\code{p.label = 1} - }{ Show all p-values on the plot. }
#'        }
#'        If a numeric vector with more than one value is 
#'        provided, they will be used as breaks for asterisk notation.
#'        Default: \code{0.05}
#'     
#' @param level   The confidence level for calculating a confidence interval. 
#'        Default: \code{0.95}
#'        
#' @param caption   Add methodology caption beneath the plot.
#'        Default: \code{TRUE}.
#'        
#' @param outliers   Show boxplot outliers? \code{TRUE} to always show. 
#'        \code{FALSE} to always hide. \code{NULL} to only hide them when
#'        overlaying a dot or strip chart.  Default: \code{NULL}
#' 
#' @param xlab.angle   Angle of the labels at the bottom of the plot. 
#'        Options are \code{"auto"}, \code{'0'}, \code{'30'}, and \code{'90'}. 
#'        Default: \code{"auto"}.
#'        
#' @param k   Number of ordination dimensions to return. Either \code{2L} or 
#'        \code{3L}. Default: \code{2L}
#'        
#' @param split.by   Name(s) of metadata columns that the data should be split
#'        by prior to any calculations. Default: \code{NULL}
#' 
#' @param dm   A \code{dist}-class distance matrix, as returned from 
#'        [bdiv_distmat()] or [stats::dist()]. Required.
#'        
#' @param groups  A named vector of grouping values. The names should 
#'        correspond to \code{attr(dm, 'Labels')}. Values can be either 
#'        categorical or numeric. Required.
#' 
#' @param stat.by   The categorical metadata field defining the statistical 
#'        groups. Required.
#' 
#' @param regr   To run a regression analysis, set \code{regr} to the numeric 
#'        metadata field with the "x-axis" values. Leaving \code{regr=NULL} 
#'        will generate boxplot-like statistics; when non-NULL, corrplot-like 
#'        statistics will be returned. Default: \code{NULL}
#'        
#' @param seed  Random seed for permutations. Default: \code{0}
#'        
#' @param permutations  Number of random permutations to use. 
#'        Default: \code{999}
#' 
#' @param p.adj   Method to use for multiple comparisons adjustment of 
#'        p-values. Run \code{p.adjust.methods} for a list of available 
#'        options. Default: \code{"fdr"}.
#' 
#' @param depths   Rarefaction depths to show in the plot, or \code{NULL} to 
#'        auto-select. Default: \code{NULL}.
#'        
#' @param rline   Where to draw a horizontal line on the plot, intended to show
#'        a particular rarefaction depth. Set to \code{TRUE} to show an 
#'        auto-selected rarefaction depth or \code{FALSE} to not show a line.
#'        Default: \code{NULL}.
#' 
#' @param labels   Show sample names under each bar. Default: \code{FALSE}.
#'   
NULL



# taxa - 4 ====
#' documentation_taxa.4
#' 
#' @name documentation_taxa.4
#' @keywords internal
#' 
#' @param taxa   Which taxa to display. An integer value will show the top n
#'        most abundant taxa. A value 0 <= n < 1 will show any taxa with that 
#'        mean abundance or greater (e.g. \code{0.1} implies >= 10%). A 
#'        character vector of taxa names will show only those named taxa. 
#'        Default: \code{4}.
#' 
NULL


# cmp ====
#' documentation_cmp
#' 
#' @name documentation_cmp
#' @keywords internal
#' 
#' @param x   A categorical metadata column name. Prefix the column name with 
#'        \code{==} or \code{!=} to limit comparisons to within or between
#'        groups, respectively. The default, \code{NULL} groups all distances 
#'        into a single column.
#'        
#' @param color.by,pattern.by,shape.by,facet.by,limit.by   Metadata columns to 
#'        use for data partitioning. Prefix the column name with 
#'        \code{==} or \code{!=} to limit comparisons to within or between
#'        groups, respectively. Default: \code{NULL}
#'        
#' @param within,between   Metadata field(s) for intra- or inter- sample 
#'        comparisons. Default: \code{NULL}
#'          
NULL




# boxplot ====
#' documentation_boxplot
#' 
#' @name documentation_boxplot
#' @keywords internal
#' 
#' @inherit documentation_plot_stats_return return
#' 
#' 
#' @param x   A categorical metadata column name to use for the x-axis. The 
#'        default, \code{NULL}, groups all samples into a single category. 
#'        
#' @param layers   \code{"bar"}, \code{"box" ("x")}, \code{"violin"}, 
#'        \code{"dot"}, \code{"strip"}, \code{"crossbar"}, \code{"errorbar"}, 
#'        \code{"linerange"}, and \code{"pointrange"}. Single letter 
#'        abbreviations are also accepted. For instance, \code{c("box", "dot")} 
#'        is equivalent to \code{c("x", "d")} and \code{"xd"}.
#'        See \code{vignette("boxplots")} for examples of each.
#'        Default: \code{"bld"}.
#'        
#' @param color.by,pattern.by,shape.by,facet.by,limit.by   Metadata columns to 
#'        use for aesthetics and partitioning. Default: \code{NULL}
#'     
#' @param ci   How to calculate min/max of the \bold{crossbar}, 
#'        \bold{errorbar}, \bold{linerange}, and \bold{pointrange} layers.
#'        Options are: \code{"ci"} (confidence interval), \code{"range"}, 
#'        \code{"sd"} (standard deviation), \code{"se"} (standard error), and 
#'        \code{"mad"} (median absolute deviation). 
#'        The center mark of \bold{crossbar} and \bold{pointrange} represents
#'        the mean, except for code{"mad"} in which case it represents
#'        the median. Default: \code{"ci"}
#'        
#' @param ...   Parameters are matched to formal arguments of ggplot2
#'        functions. Prefixing parameter names with a layer name ensures that
#'        a particular parameter is passed to, and only to, that layer. For
#'        instance, \code{dot.size = 2} or \code{d.size = 2} ensures only the 
#'        dotplot layer has its size set to \code{2}. The special prefix
#'        \code{pt.} will control both the dot and strip layers.
#' 
NULL




# boxplot aes - section ====
#' documentation_boxplot_aes_section
#' 
#' @name documentation_boxplot_aes_section
#' @keywords internal
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
#'   \item{\code{facet.by} - }{Categorical metadata column(s).}
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
#' 
NULL



# corrplot ====
#' documentation_corrplot
#' 
#' @name documentation_corrplot
#' @keywords internal
#' 
#' @inherit documentation_plot_stats_return return
#' 
#' 
#' @param x   A numeric metadata field to use for the x-axis. Required.
#'           
#' @param layers   \code{"trend"}, \code{"scatter"}. Single letter 
#'        abbreviations are also accepted. For instance, 
#'        \code{c("trend", "scatter")} is equivalent to \code{c("t", "s")} and 
#'        \code{"ts"}. See \code{vignette("corrplots")} for examples of each. 
#'        Default: \code{"t"}.
#'        
#' @param color.by,facet.by,limit.by   Metadata columns to use for aesthetics 
#'        and partitioning. See below for details. Default: \code{NULL}
#'        
#' @param test   Which test statistic to display on the plot. Options are: 
#'        \itemize{
#'          \item{\code{"fit"} - }{ How well does the model fit the data? }
#'          \item{\code{"terms"} - }{ How strongly does 'x' influence 'y'? }
#'          \item{\code{"means"} - }{ Is the average 'y' value non-zero? }
#'          \item{\code{"trends"} - }{ Does any trendline have a non-zero slope? }
#'          \item{\code{"pw_means"} - }{ Are the means of any trendlines different? }
#'          \item{\code{"pw_trends"} - }{ Are the slopes of any trendlines different? }
#'          \item{\code{"none"} - }{ Do not compute or show statistics. }
#'        }
#'        Default: \code{"trends"} \cr\cr
#'        Note: \code{"pw_means"} and \code{"pw_trends"} can only be calculated
#'        when using a \code{color.by} metadata column with more than one level. \cr\cr
#'        Statistical tests are run separately on each facet. P-values are 
#'        adjusted for multiple comparisons by considering all facets together.
#'        
#' @param ...   Additional parameters to pass along to ggplot2
#'        functions. Prefix a parameter name with either \code{t.} or 
#'        \code{s.}/\code{pt.} to ensure it gets passed to (and only to) 
#'        \link[ggplot2]{geom_smooth} or \link[ggplot2]{geom_point}, 
#'        respectively. For instance, \code{s.size = 2} ensures only the 
#'        scatterplot points have their size set to \code{2}.
#' 
#' 
#' @section Aesthetics and Partitions:
#' 
#' Metadata can be used to flexibly subset, partition, and apply aesthetics 
#' when creating a plot. Common use cases are provided below. More thorough 
#' documentation is available at \url{https://cmmr.github.io/rbiom}.
#' 
#' \preformatted{  ## Colors ----------------------------
#'   color.by = "Body Site"
#'   color.by = list('Body Site' = "bright")
#'   color.by = list('Body Site' = c("Stool", "Saliva"))
#'   color.by = list('Body Site' = list('values' = c("Stool", "Saliva"), 'colors' = "bright"))
#'   color.by = list('Body Site' = c('Stool' = "blue", 'Saliva' = "green"))
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
#'   \item{\code{color.by} - }{Any metadata column. (Max 1)}
#'   \item{\code{facet.by} - }{Only categorical metadata column(s).}
#'   \item{\code{limit.by} - }{Any metadata column(s).}
#' }
#' 
#' All built-in color palettes are colorblind-friendly.
#' 
#' The available categorical palette names are: \code{"okabe"}, \code{"carto"}, 
#' \code{"r4"}, \code{"polychrome"}, \code{"tol"}, \code{"bright"}, 
#' \code{"light"}, \code{"muted"}, \code{"vibrant"}, \code{"tableau"}, 
#' \code{"classic"}, \code{"alphabet"}, \code{"tableau20"}, \code{"kelly"}, 
#' and \code{"fishy"}.
#' 
#' The available numeric palette names are: \code{"reds"}, \code{"oranges"}, 
#' \code{"greens"}, \code{"purples"}, \code{"grays"}, \code{"acton"}, 
#' \code{"bamako"}, \code{"batlow"}, \code{"bilbao"}, \code{"buda"}, 
#' \code{"davos"}, \code{"devon"}, \code{"grayC"}, \code{"hawaii"}, 
#' \code{"imola"}, \code{"lajolla"}, \code{"lapaz"}, \code{"nuuk"}, 
#' \code{"oslo"}, \code{"tokyo"}, \code{"turku"}, \code{"bam"}, 
#' \code{"berlin"}, \code{"broc"}, \code{"cork"}, \code{"lisbon"}, 
#' \code{"roma"}, \code{"tofino"}, \code{"vanimo"}, \code{"vik"}
#' 
NULL



# heatmap ====
#' documentation_heatmap
#' 
#' @name documentation_heatmap
#' @keywords internal
#' 
#' @inherit documentation_plot_return return
#' 
#' @param grid   Color palette name, or a list with entries for \code{label}, 
#'        \code{colors}, \code{range}, \code{bins}, \code{na.color}, and/or 
#'        \code{guide}. See the Track Definitions section for details.
#'        Default: \code{list(label = "Grid Value", colors = "imola")}.
#'        
#' @param label   Label the matrix rows and columns. You can supply a list
#'        or logical vector of length two to control row labels and column 
#'        labels separately, for example 
#'        \code{label = c(rows = TRUE, cols = FALSE)}, or simply 
#'        \code{label = c(T, F)}. Other valid options are \code{"rows"},
#'        \code{"cols"}, \code{"both"}, \code{"bottom"}, \code{"right"},
#'        and \code{"none"}.
#'        Default: \code{TRUE}.
#'        
#' @param label_size   The font size to use for the row and column labels. You 
#'        can supply a numeric vector of length two to control row label sizes 
#'        and column label sizes separately, for example 
#'        \code{c(rows = 20, cols = 8)}, or simply \code{c(20, 8)}.
#'        Default: \code{NULL}, which computes: 
#'        \code{pmax(8, pmin(20, 100 / dim(mtx)))}.
#'        
#' @param rescale   Rescale rows or columns to all have a common min/max.
#'        Options: \code{"none"}, \code{"rows"}, or \code{"cols"}.
#'        Default: \code{"none"}.
#'        
#' @param trees   Draw a dendrogram for rows (left) and columns (top). You can 
#'        supply a list or logical vector of length two to control the row tree 
#'        and column tree separately, for example 
#'        \code{trees = c(rows = T, cols = F)}, or simply \code{trees = c(T, F)}. 
#'        Other valid options are \code{"rows"}, \code{"cols"}, \code{"both"}, 
#'        \code{"left"}, \code{"top"}, and \code{"none"}.
#'        Default: \code{TRUE}.
#'        
#' @param clust   Clustering algorithm for reordering the rows and columns by 
#'        similarity. You can supply a list or character vector of length two to 
#'        control the row and column clustering separately, for example 
#'        \code{clust = c(rows = "complete", cols = NA)}, or simply 
#'        \code{clust = c("complete", NA)}. Options are:
#'        \itemize{
#'          \item{\code{FALSE} or \code{NA} - }{ Disable reordering. }
#'          \item{An \code{hclust} class object}{ E.g. from [stats::hclust()]. }
#'          \item{A method name - }{ \code{"ward.D"}, 
#'            \code{"ward.D2"}, \code{"single"}, \code{"complete"}, 
#'            \code{"average"}, \code{"mcquitty"}, \code{"median"}, or 
#'            \code{"centroid"}. }
#'        }
#'        Default: \code{"complete"}.
#'        
#' @param dist   Distance algorithm to use when reordering the rows and columns 
#'        by similarity. You can supply a list or character vector of length
#'        two to control the row and column clustering separately, for example 
#'        \code{dist = c(rows = "euclidean", cols = "maximum")}, or simply 
#'        \code{dist = c("euclidean", "maximum")}. Options are:
#'        \itemize{
#'          \item{A \code{dist} class object}{ E.g. from [stats::dist()] or [bdiv_distmat()]. }
#'          \item{A method name - }{ \code{"euclidean"}, 
#'            \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, 
#'            \code{"binary"}, or \code{"minkowski"}. }
#'        }
#'        Default: \code{"euclidean"}.
#'        
#' @param tree_height,track_height   The height of the dendrogram or annotation
#'        tracks in multiples (or fractions) of the smaller dimension of the
#'        grid cell size. Use a numeric vector of length two to assign
#'        \code{c(top, left)} independently. 
#'        Default: \code{NULL}, which computes:
#'        \code{tree_height = sqrt(min(dim(mtx))), track_height = tree_height / 4}.
#'        
#' @param ratio   Height/width ratio for entire grid. 
#'        Default: \code{1} (square).
#'        
#' @param legend   Where to place the legend. Options are: \code{"right"} or
#'        \code{"bottom"}. Default: \code{"right"}.
#'        
#' @param ...   Additional arguments to pass on to ggplot2::theme().
#' 
NULL



# dist test ====
#' documentation_dist_test
#' 
#' @name documentation_dist_test
#' @keywords internal
#' 
#' @param stat.by   The categorical or numeric metadata field over which 
#'        statistics should be calculated. Required.
#'        
#' @param test    Permutational test for accessing significance. Options are:
#'        \itemize{
#'            \item{\code{"adonis2"} - }{ Permutational MANOVA; [vegan::adonis2()]. }
#'            \item{\code{"mrpp"} - }{ Multiple response permutation procedure; [vegan::mrpp()]. }
#'        }
#'        Default: \code{"adonis2"} \cr\cr
#'        Abbreviations are allowed.
#'        
NULL




# test - ifelse ====
#' documentation_test.ifelse
#' 
#' @name documentation_test.ifelse
#' @keywords internal
#' 
#' 
#' @param test   The statistic of interest. An overview of options can be 
#'        found in the "Statistical Tests" section below. One of
#'        \code{"predict"}, \code{"terms"}, \code{"fit"}, \code{"means"}, 
#'        \code{"trends"}, \code{"es_means"}, \code{"es_trends"}, 
#'        \code{"pw_means"}, or \code{"pw_trends"}. 
#'        Default: \code{ifelse(is.null(regr), "means", "trends")}
#' 
#' 
#' @section Statistical Tests:
#' 
#' \bold{When \code{regr} is NULL:}
#' \itemize{
#'   \item{\code{"means"} - }{ 
#'     Considers all groups at once using the Kruskal-Wallis
#'     non-parametric test: [stats::kruskal.test()]. }
#'   \item{\code{"pw_means"} - }{
#'     Pairwise comparison of \code{stat.by} groups using the Wilcox rank sum 
#'     (aka Mann-Whitney) non-parametric test: [stats::wilcox.test()]. }
#' }
#' 
#' \bold{When \code{regr} is not NULL:}
#' \itemize{
#'   \item{\code{"predict"} - }{
#'     Augments original data with fitted information. See 
#'     [broom::augment.lm()]. }
#'   \item{\code{"terms"} - }{ Summary information about the model's 
#'     terms, including p-value, r-squared, AIC, BIC, and confidence 
#'     interval at \code{level}. See [broom::tidy.lm()]. }
#'   \item{\code{"fit"} - }{ Goodness of fit measures, p-values, and 
#'     more for the overall model. See [broom::glance.lm()]. }
#'   \item{\code{"means"} - }{ The estimated marginal mean (EMM) for 
#'     each \code{stat.by} group, along with confidence intervals (CI), 
#'     standard errors (SE), t-ratios, and p-values testing for mean = zero. 
#'     See [emmeans::emmeans()] and [emmeans::summary.emmGrid()]. }
#'   \item{\code{"trends"} - }{ The trendline slope for each 
#'     \code{stat.by} group, along with CI and SE and p-value testing for
#'     slope = zero. See [emmeans::emtrends()] and 
#'     [emmeans::summary.emmGrid()]. }
#'   \item{\code{"pw_means"} - }{ Pairwise means. All \code{stat.by} groups are 
#'     compared to each other and the difference in means is estimated 
#'     along with SE, t-ratios, and p-values testing if the two means 
#'     are the same. See [emmeans::pairs.emmGrid()]. }
#'   \item{\code{"pw_trends"} - }{ Pairwise trends. As above, but
#'     comparing trendline slopes instead of means. }
#'   \item{\code{"es_means"} - }{ Effect sizes for pairwise differences 
#'     of means, including SE and CI. See [emmeans::eff_size()]. }
#'   \item{\code{"es_trends"} - }{ Effect sizes for pairwise differences 
#'     of slopes, including SE and CI. See [emmeans::eff_size()]. }
#' }
#'        
NULL



# test - pw_means ====
#' documentation_test.pw_means
#' 
#' @name documentation_test.pw_means
#' @keywords internal
#' 
#' 
#' @param test   The statistic of interest. An overview of options can be 
#'        found in the "Statistical Tests" section below. One of
#'        \code{"predict"}, \code{"terms"}, \code{"fit"}, \code{"means"}, 
#'        \code{"trends"}, \code{"es_means"}, \code{"es_trends"}, 
#'        \code{"pw_means"}, or \code{"pw_trends"}. 
#'        Default: \code{"pw_means"}
#' 
#' 
#' @section Statistical Tests:
#' 
#' \itemize{
#'   \item{\code{"predict"} - }{
#'     Augments original data with fitted information. See 
#'     [broom::augment.lm()]. }
#'   \item{\code{"terms"} - }{ Summary information about the model's 
#'     terms, including p-value, r-squared, AIC, BIC, and confidence 
#'     interval at \code{level}. See [broom::tidy.lm()]. }
#'   \item{\code{"fit"} - }{ Goodness of fit measures, p-values, and 
#'     more for the overall model. See [broom::glance.lm()]. }
#'   \item{\code{"means"} - }{ The estimated marginal mean (EMM) for 
#'     each \code{stat.by} group, along with confidence intervals (CI), 
#'     standard errors (SE), t-ratios, and p-values testing for mean = zero. 
#'     See [emmeans::emmeans()] and [emmeans::summary.emmGrid()]. }
#'   \item{\code{"trends"} - }{ The trendline slope for each 
#'     \code{stat.by} group, along with CI and SE and p-value testing for
#'     slope = zero. See [emmeans::emtrends()] and 
#'     [emmeans::summary.emmGrid()]. }
#'   \item{\code{"pw_means"} - }{ Pairwise means. All \code{stat.by} groups are 
#'     compared to each other and the difference in means is estimated 
#'     along with SE, t-ratios, and p-values testing if the two means 
#'     are the same. See [emmeans::pairs.emmGrid()]. }
#'   \item{\code{"pw_trends"} - }{ Pairwise trends. As above, but
#'     comparing trendline slopes instead of means. }
#'   \item{\code{"es_means"} - }{ Effect sizes for pairwise differences 
#'     of means, including SE and CI. See [emmeans::eff_size()]. }
#'   \item{\code{"es_trends"} - }{ Effect sizes for pairwise differences 
#'     of slopes, including SE and CI. See [emmeans::eff_size()]. }
#' }
#' 
NULL



# test - trends ====
#' documentation_test.trends
#' 
#' @name documentation_test.trends
#' @keywords internal
#' 
#' 
#' @param test   The statistic of interest. An overview of options can be 
#'        found in the "Statistical Tests" section below. One of
#'        \code{"predict"}, \code{"terms"}, \code{"fit"}, \code{"means"}, 
#'        \code{"trends"}, \code{"es_means"}, \code{"es_trends"}, 
#'        \code{"pw_means"}, or \code{"pw_trends"}. 
#'        Default: \code{"trends"}
#' 
#' 
#' @section Statistical Tests:
#' 
#' \itemize{
#'   \item{\code{"predict"} - }{
#'     Augments original data with fitted information. See 
#'     [broom::augment.lm()]. }
#'   \item{\code{"terms"} - }{ Summary information about the model's 
#'     terms, including p-value, r-squared, AIC, BIC, and confidence 
#'     interval at \code{level}. See [broom::tidy.lm()]. }
#'   \item{\code{"fit"} - }{ Goodness of fit measures, p-values, and 
#'     more for the overall model. See [broom::glance.lm()]. }
#'   \item{\code{"means"} - }{ The estimated marginal mean (EMM) for 
#'     each \code{stat.by} group, along with confidence intervals (CI), 
#'     standard errors (SE), t-ratios, and p-values testing for mean = zero. 
#'     See [emmeans::emmeans()] and [emmeans::summary.emmGrid()]. }
#'   \item{\code{"trends"} - }{ The trendline slope for each 
#'     \code{stat.by} group, along with CI and SE and p-value testing for
#'     slope = zero. See [emmeans::emtrends()] and 
#'     [emmeans::summary.emmGrid()]. }
#'   \item{\code{"pw_means"} - }{ Pairwise means. All \code{stat.by} groups are 
#'     compared to each other and the difference in means is estimated 
#'     along with SE, t-ratios, and p-values testing if the two means 
#'     are the same. See [emmeans::pairs.emmGrid()]. }
#'   \item{\code{"pw_trends"} - }{ Pairwise trends. As above, but
#'     comparing trendline slopes instead of means. }
#'   \item{\code{"es_means"} - }{ Effect sizes for pairwise differences 
#'     of means, including SE and CI. See [emmeans::eff_size()]. }
#'   \item{\code{"es_trends"} - }{ Effect sizes for pairwise differences 
#'     of slopes, including SE and CI. See [emmeans::eff_size()]. }
#' }
#' 
NULL



# model - lm ====
#' documentation_model.lm
#' 
#' @name documentation_model.lm
#' @keywords internal
#' 
#' @inherit documentation_model_section sections
#' 
#' @param model   What type of trend model to fit to the data. Options are: 
#'        \code{"lm"} (linear), \code{"log"} (logarithmic), or \code{"gam"} 
#'        (generalized additive). See the "Model Options" section below for
#'        additional details. Default: \code{"lm"}
#' 
NULL



# model - log ====
#' documentation_model.log
#' 
#' @name documentation_model.log
#' @keywords internal
#' 
#' @inherit documentation_model_section sections
#' 
#' @param model   What type of trend model to fit to the data. Options are: 
#'        \code{"lm"} (linear), \code{"log"} (logarithmic), or \code{"gam"} 
#'        (generalized additive). See the "Model Options" section below for
#'        additional details. Default: \code{"log"}
#' 
NULL



# model - section ====
#' documentation_model_section
#' 
#' @name documentation_model_section
#' @keywords internal
#' 
#' 
#' @section Model Options:
#' 
#' The predefined options are: 
#' \itemize{
#'   \item{\code{"lm"} - }{  Linear model: \code{stats::lm(formula = y ~ x)}.) }
#'   \item{\code{"log"} - }{ Logarithmic model: \code{stats::lm(formula = y ~ log(x))}. }
#'   \item{\code{"gam"} - }{ Generalized additive model: \code{mgcv::gam(formula = y ~ s(x, bs = "cs"), method = "REML")}. }
#' }
#' 
#' You can alternatively provide a list of length two where the first 
#' element is a character vector of length 1 naming a function, and the 
#' second element is a list of arguments to pass to that function. One 
#' of the function's arguments must be named 'formula'. 
#' For example, \code{model = list("stats::lm", list(formula = y ~ x))}.
#' 
NULL




# plot - return ====
#' documentation_plot_return
#' 
#' @name documentation_plot_return
#' @keywords internal
#' 
#' @return A \code{ggplot2} plot. \cr The computed data points, ggplot command, 
#'         and object history are available as \code{$data}, \code{$code}, and 
#'         \code{$history}, respectively.
#' 
NULL


# plot w/ stats - return ====
#' documentation_plot_stats_return
#' 
#' @name documentation_plot_stats_return
#' @keywords internal
#' 
#' @return A \code{ggplot2} plot. \cr The computed data points, statistics, 
#'         ggplot command, and object history are available as \code{$data}, 
#'         \code{$stats}, \code{$code}, and \code{$history}, respectively.
#' 
NULL


# stats - return ====
#' documentation_stats_return
#' 
#' @name documentation_stats_return
#' @keywords internal
#' 
#' @return A tibble data frame with summary statistics. \cr
#'         The R code or generating these statistics is in \code{$code}, and 
#'         the object history is in \code{$history}.
#' 
NULL

