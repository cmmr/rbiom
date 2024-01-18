
# default ====
#' documentation_default
#' 
#' @name documentation_default
#' @keywords internal
#' 
#' @param biom   An rbiom-class object, or data coercible with 
#'        [as_rbiom()].
#' 
#' @param mtx   A matrix-like object.
#'        
#' @param tree  A `phylo` object representing the phylogenetic
#'        relationships of the taxa in `biom`. Only required when 
#'        computing UniFrac distances. Default: `biom$tree`
#'     
#' @param md  A character vector naming the metadata fields to include in the 
#'        output data frame, or `'.all'` to include all metadata fields.
#'        Default: `'.all'`
#' 
#' @param adiv   Alpha diversity metric(s) to use. Options are: `"OTUs"`, 
#'        `"Shannon"`, `"Chao1"`, `"Simpson"`, and/or 
#'        `"InvSimpson"`. Set `adiv=".all"` to use all metrics.
#'        Default: `"Shannon"` \cr\cr
#'        Multiple values allowed. Non-ambiguous abbreviations are allowed.
#' 
#' @param bdiv  Beta diversity distance algorithm(s) to use. Options are:
#'        `"Bray-Curtis"`, `"Manhattan"`, `"Euclidean"`, 
#'        `"Jaccard"`, and `"UniFrac"`. For `"UniFrac"`, a 
#'        phylogenetic tree must be present in `biom` or explicitly 
#'        provided via `tree=`. Default: `"Bray-Curtis"` \cr\cr
#'        Multiple/abbreviated values allowed.
#' 
#' @param taxa   Which taxa to display. An integer value will show the top n
#'        most abundant taxa. A value 0 <= n < 1 will show any taxa with that 
#'        mean abundance or greater (e.g. `0.1` implies >= 10%). A 
#'        character vector of taxa names will show only those named taxa. 
#'        Default: `6`.
#'        
#' @param ord    Method for reducing dimensionality. Options are:
#'        \itemize{
#'            \item{`"UMAP"` - }{ Uniform manifold approximation and projection; [uwot::umap()]. }
#'            \item{`"PCoA"` - }{ Principal coordinate analysis; [ape::pcoa()]. }
#'            \item{`"NMDS"` - }{ Nonmetric multidimensional scaling; [vegan::metaMDS()]. }
#'            \item{`"tSNE"` - }{ t-distributed stochastic neighbor embedding; [tsne::tsne()]. }
#'        }
#'        Default: `"UMAP"` \cr\cr
#'        Multiple values allowed. Non-ambiguous abbreviations are allowed.
#'     
#' @param weighted  Take relative abundances into account. When 
#'        `weighted=FALSE`, only presence/absence is considered.
#'        Default: `TRUE` \cr\cr
#'        Multiple values allowed.
#'        
#' @param rank   What rank(s) of taxa to display. E.g. `"Phylum"`, 
#'        `"Genus"`, `".otu"`, etc. An integer vector can also be 
#'        given, where `1` is the highest rank, `2` is the second 
#'        highest, `-1` is the lowest rank, `-2` is the second 
#'        lowest, and `0` is the OTU "rank". Run `biom$ranks` to 
#'        see all options for a given rbiom object. Default: `-1`.
#'        
#' @param lineage  Include all ranks in the name of the taxa. For instance,
#'        setting to `TRUE` will produce 
#'        `Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales`. 
#'        Otherwise the taxa name will simply be `Coriobacteriales`. You 
#'        want to set this to TRUE when `unc = "asis"` and you have taxa 
#'        names (such as \emph{Incertae_Sedis}) that map to multiple higher 
#'        level ranks. Default: `FALSE`
#'        
#' @param unc  How to handle unclassified, uncultured, and similarly ambiguous
#'        taxa names. Options are: 
#'        \itemize{
#'          \item{`"singly"` - }{ Replaces them with the OTU name. }
#'          \item{`"grouped"` - }{ Replaces them with a higher rank's name. }
#'          \item{`"drop"` - }{ Excludes them from the result. }
#'          \item{`"asis"` - }{ To not check/modify any taxa names. }
#'        }
#'        Default: `"singly"` \cr\cr
#'        Abbreviations are allowed.
#'        
#' @param other  Sum all non-itemized taxa into an "Other" taxa. When 
#'        `FALSE`, only returns taxa matched by the `taxa` 
#'        argument. Specifying `TRUE` adds "Other" to the returned set.
#'        A string can also be given to imply `TRUE`, but with that
#'        value as the name to use instead of "Other".
#'        Default: `FALSE`
#'        
#' @param sparse  If true, returns a sparse matrix as described by 
#'        [slam::simple_triplet_matrix()], otherwise returns a normal R
#'        matrix object. Default: `FALSE`
#'        
#' @param p.top   Only display taxa with the most significant differences in 
#'        abundance. If `p.top` is >= 1, then the `p.top` most 
#'        significant taxa are displayed. If `p.top` is less than one, all 
#'        taxa with an adjusted p-value <= `p.top` are displayed. 
#'        Recommended to be used in combination with the `taxa` parameter 
#'        to set a lower bound on the mean abundance of considered taxa. 
#'        Default: `Inf`
#'        
#' @param y.trans   The transformation to apply to the y-axis. Visualizing 
#'        differences of both high- and low-abundance taxa is best done with
#'        a non-linear axis. Options are: 
#'        \itemize{
#'          \item{`"sqrt"` - }{ square-root transformation }
#'          \item{`"log1p"` - }{ log(y + 1) transformation }
#'          \item{`NULL` - }{ no transformation }
#'        }
#'        These methods allow visualization of both high- and low-abundance
#'        taxa simultaneously, without complaint about 'zero' count
#'        observations. Default: `"sqrt"`
#'        
#' @param flip   Transpose the axes, so that taxa are present as rows instead
#'        of columns. Default: `FALSE`
#'
#' @param stripe   Shade every other x position. Default: \emph{same as flip}
#'
#' @param p.label   Minimum adjusted p-value to display on the plot with a 
#'        bracket.
#'        \itemize{
#'          \item{`p.label = 0.05` - }{ Show p-values that are <= 0.05. }
#'          \item{`p.label = 0` - }{ Don't show any p-values on the plot. }
#'          \item{`p.label = 1` - }{ Show all p-values on the plot. }
#'        }
#'        If a numeric vector with more than one value is 
#'        provided, they will be used as breaks for asterisk notation.
#'        Default: `0.05`
#'     
#' @param level   The confidence level for calculating a confidence interval. 
#'        Default: `0.95`
#'        
#' @param caption   Add methodology caption beneath the plot.
#'        Default: `TRUE`
#'        
#' @param outliers   Show boxplot outliers? `TRUE` to always show. 
#'        `FALSE` to always hide. `NULL` to only hide them when
#'        overlaying a dot or strip chart.  Default: `NULL`
#' 
#' @param xlab.angle   Angle of the labels at the bottom of the plot. 
#'        Options are `"auto"`, `'0'`, `'30'`, and `'90'`. 
#'        Default: `"auto"`.
#'        
#' @param k   Number of ordination dimensions to return. Either `2L` or 
#'        `3L`. Default: `2L`
#'        
#' @param split.by   Name(s) of metadata columns that the data should be split
#'        by prior to any calculations. Default: `NULL`
#' 
#' @param dm   A `dist`-class distance matrix, as returned from 
#'        [bdiv_distmat()] or [stats::dist()]. Required.
#'        
#' @param groups  A named vector of grouping values. The names should 
#'        correspond to `attr(dm, 'Labels')`. Values can be either 
#'        categorical or numeric. Required.
#' 
#' @param stat.by   The categorical metadata field defining the statistical 
#'        groups. Required.
#' 
#' @param regr   To run a regression analysis, set `regr` to the numeric 
#'        metadata field with the "x-axis" values. Leaving `regr=NULL` 
#'        will generate boxplot-like statistics; when non-NULL, corrplot-like 
#'        statistics will be returned. Default: `NULL`
#'        
#' @param seed  Random seed for permutations. Default: `0`
#'        
#' @param permutations  Number of random permutations to use. 
#'        Default: `999`
#' 
#' @param p.adj   Method to use for multiple comparisons adjustment of 
#'        p-values. Run `p.adjust.methods` for a list of available 
#'        options. Default: `"fdr"`
#' 
#' @param depths   Rarefaction depths to show in the plot, or `NULL` to 
#'        auto-select. Default: `NULL`
#'        
#' @param rline   Where to draw a horizontal line on the plot, intended to show
#'        a particular rarefaction depth. Set to `TRUE` to show an 
#'        auto-selected rarefaction depth or `FALSE` to not show a line.
#'        Default: `NULL`
#'        
#' @param clone   Create a copy of `biom` before modifying. If `FALSE`, `biom` 
#'        is modified in place as a side-effect. See [speed ups][speed] for 
#'        use cases. Default: `TRUE`
#' 
#' @param labels   Show sample names under each bar. Default: `FALSE`
#' 
#' @param trans   Transformation to apply. Options are: 
#'        `c("none", "rank", "log", "log1p", "sqrt")`. `"rank"` is useful for 
#'        correcting for non-normally distributions before applying regression 
#'        statistics. Default: `"none"`
#'   
NULL



# biom - rbiom object ====
#' documentation_biom.rbiom
#' 
#' @name documentation_biom.rbiom
#' @keywords internal
#' 
#' @param biom    An [rbiom object][rbiom_objects], such as from [as_rbiom()].
#' 
NULL



# return - biom ====
#' documentation_return.biom
#' 
#' @name documentation_return.biom
#' @keywords internal
#' 
#' @return An [rbiom object][rbiom_objects].
#' 
NULL



# rank - NULL ====
#' documentation_rank.NULL
#' 
#' @name documentation_rank.NULL
#' @keywords internal
#'        
#' @param rank   What rank(s) of taxa to compute biplot coordinates and 
#'        statistics for, or `NULL` to disable. E.g. `"Phylum"`, 
#'        `"Genus"`, `".otu"`, etc. An integer vector can also be 
#'        given, where `1` is the highest rank, `2` is the second 
#'        highest, `-1` is the lowest rank, `-2` is the second 
#'        lowest, and `0` is the OTU "rank". Run `biom$ranks` to 
#'        see all options for a given rbiom object. Default: `NULL`.
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
#'        mean abundance or greater (e.g. `0.1` implies >= 10%). A 
#'        character vector of taxa names will show only those named taxa. 
#'        Default: `4`.
#' 
NULL


# cmp ====
#' documentation_cmp
#' 
#' @name documentation_cmp
#' @keywords internal
#' 
#' @param x   A categorical metadata column name. Prefix the column name with 
#'        `==` or `!=` to limit comparisons to within or between
#'        groups, respectively. The default, `NULL` groups all distances 
#'        into a single column.
#'        
#' @param color.by,pattern.by,shape.by,facet.by,limit.by   Metadata columns to 
#'        use for data partitioning. Prefix the column name with 
#'        `==` or `!=` to limit comparisons to within or between
#'        groups, respectively. Default: `NULL`
#'        
#' @param within,between   Metadata field(s) for intra- or inter- sample 
#'        comparisons. Default: `NULL`
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
#'        default, `NULL`, groups all samples into a single category. 
#'        
#' @param layers   One or more of 
#'        `c("bar", "box" ("x"), "violin", "dot", "strip", "crossbar", "errorbar", "linerange", "pointrange")`. 
#'        Single letter abbreviations are also accepted. For instance, 
#'        `c("box", "dot")` is equivalent to `c("x", "d")` and `"xd"`.
#'        See [plot types][plots] for examples of each. Default: `"bld"`
#'        
#' @param color.by,pattern.by,shape.by,facet.by,limit.by   Metadata columns to 
#'        use for aesthetics and partitioning. Default: `NULL`
#'     
#' @param ci   How to calculate min/max of the \bold{crossbar}, 
#'        \bold{errorbar}, \bold{linerange}, and \bold{pointrange} layers.
#'        Options are: `"ci"` (confidence interval), `"range"`, 
#'        `"sd"` (standard deviation), `"se"` (standard error), and 
#'        `"mad"` (median absolute deviation). 
#'        The center mark of \bold{crossbar} and \bold{pointrange} represents
#'        the mean, except for code{"mad"} in which case it represents
#'        the median. Default: `"ci"`
#'        
#' @param ...   Parameters are matched to formal arguments of ggplot2
#'        functions. Prefixing parameter names with a layer name ensures that
#'        a particular parameter is passed to, and only to, that layer. For
#'        instance, `dot.size = 2` or `d.size = 2` ensures only the 
#'        dotplot layer has its size set to `2`. The special prefix
#'        `pt.` will control both the dot and strip layers.
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
#'   \item{`color.by` - }{A categorical metadata column. (Max 1)}
#'   \item{`pattern.by` - }{A categorical metadata column. (Max 1)}
#'   \item{`shape.by` - }{A categorical metadata column. (Max 1)}
#'   \item{`facet.by` - }{Categorical metadata column(s).}
#'   \item{`limit.by` - }{Any metadata column(s).}
#' }
#' 
#' All built-in color palettes are colorblind-friendly. The available 
#' categorical palette names are: `"okabe"`, `"carto"`, `"r4"`, 
#' `"polychrome"`, `"tol"`, `"bright"`, `"light"`, 
#' `"muted"`, `"vibrant"`, `"tableau"`, `"classic"`, 
#' `"alphabet"`, `"tableau20"`, `"kelly"`, and `"fishy"`.
#' 
#' Patterns are sourced from the magick R package. Pattern names are: 
#' `"bricks"`, `"hexagons"`, `"horizontalsaw"`, 
#' `"hs_fdiagonal"`, `"fishscales"`, `"verticalsaw"`, 
#' `"checkerboard"`, `"octagons"`, `"right45"`, 
#' `"hs_cross"`, `"hs_bdiagonal"`, `"hs_diagcross"`, 
#' `"hs_horizontal"`, `"hs_vertical"`, `"left45"`, 
#' `"leftshingle"`, `"rightshingle"`, `"verticalbricks"`, 
#' `"verticalleftshingle"`, and `"verticalrightshingle"`.
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
#' @param layers   One or more of 
#'        `c("trend", "confidence", "scatter", "residual", "name")`. Single 
#'        letter abbreviations are also accepted. For instance, 
#'        `c("trend", "scatter")` is equivalent to `c("t", "s")` and `"ts"`. 
#'        See [plot types][plots] for examples of each. Default: `"t"`
#'        
#' @param color.by,facet.by,limit.by   Metadata columns to use for aesthetics 
#'        and partitioning. See below for details. Default: `NULL`
#'        
#' @param test   Which test statistic to display on the plot. Options are: 
#'        \itemize{
#'          \item{`"fit"` - }{ How well does the model fit the data? }
#'          \item{`"terms"` - }{ How strongly does 'x' influence 'y'? }
#'          \item{`"means"` - }{ Is the average 'y' value non-zero? }
#'          \item{`"trends"` - }{ Does any trendline have a non-zero slope? }
#'          \item{`"pw_means"` - }{ Are the means of any trendlines different? }
#'          \item{`"pw_trends"` - }{ Are the slopes of any trendlines different? }
#'          \item{`"none"` - }{ Do not compute or show statistics. }
#'        }
#'        Default: `"trends"` \cr\cr
#'        Note: `"pw_means"` and `"pw_trends"` can only be calculated
#'        when using a `color.by` metadata column with more than one level. \cr\cr
#'        Statistical tests are run separately on each facet. P-values are 
#'        adjusted for multiple comparisons by considering all facets together.
#'        
#' @param ...   Additional parameters to pass along to ggplot2
#'        functions. Prefix a parameter name with either `t.` or 
#'        `s.`/`pt.` to ensure it gets passed to (and only to) 
#'        \link[ggplot2]{geom_smooth} or \link[ggplot2]{geom_point}, 
#'        respectively. For instance, `s.size = 2` ensures only the 
#'        scatterplot points have their size set to `2`.
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
#'   \item{`color.by` - }{Any metadata column. (Max 1)}
#'   \item{`facet.by` - }{Only categorical metadata column(s).}
#'   \item{`limit.by` - }{Any metadata column(s).}
#' }
#' 
#' All built-in color palettes are colorblind-friendly.
#' 
#' The available categorical palette names are: `"okabe"`, `"carto"`, 
#' `"r4"`, `"polychrome"`, `"tol"`, `"bright"`, 
#' `"light"`, `"muted"`, `"vibrant"`, `"tableau"`, 
#' `"classic"`, `"alphabet"`, `"tableau20"`, `"kelly"`, 
#' and `"fishy"`.
#' 
#' The available numeric palette names are: `"reds"`, `"oranges"`, 
#' `"greens"`, `"purples"`, `"grays"`, `"acton"`, 
#' `"bamako"`, `"batlow"`, `"bilbao"`, `"buda"`, 
#' `"davos"`, `"devon"`, `"grayC"`, `"hawaii"`, 
#' `"imola"`, `"lajolla"`, `"lapaz"`, `"nuuk"`, 
#' `"oslo"`, `"tokyo"`, `"turku"`, `"bam"`, 
#' `"berlin"`, `"broc"`, `"cork"`, `"lisbon"`, 
#' `"roma"`, `"tofino"`, `"vanimo"`, `"vik"`
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
#' @param grid   Color palette name, or a list with entries for `label`, 
#'        `colors`, `range`, `bins`, `na.color`, and/or 
#'        `guide`. See the Track Definitions section for details.
#'        Default: `list(label = "Grid Value", colors = "imola")`
#'        
#' @param label   Label the matrix rows and columns. You can supply a list
#'        or logical vector of length two to control row labels and column 
#'        labels separately, for example 
#'        `label = c(rows = TRUE, cols = FALSE)`, or simply 
#'        `label = c(T, F)`. Other valid options are `"rows"`,
#'        `"cols"`, `"both"`, `"bottom"`, `"right"`,
#'        and `"none"`.
#'        Default: `TRUE`
#'        
#' @param label_size   The font size to use for the row and column labels. You 
#'        can supply a numeric vector of length two to control row label sizes 
#'        and column label sizes separately, for example 
#'        `c(rows = 20, cols = 8)`, or simply `c(20, 8)`.
#'        Default: `NULL`, which computes: 
#'        `pmax(8, pmin(20, 100 / dim(mtx)))`
#'        
#' @param rescale   Rescale rows or columns to all have a common min/max.
#'        Options: `"none"`, `"rows"`, or `"cols"`.
#'        Default: `"none"`
#'        
#' @param trees   Draw a dendrogram for rows (left) and columns (top). You can 
#'        supply a list or logical vector of length two to control the row tree 
#'        and column tree separately, for example 
#'        `trees = c(rows = T, cols = F)`, or simply `trees = c(T, F)`. 
#'        Other valid options are `"rows"`, `"cols"`, `"both"`, 
#'        `"left"`, `"top"`, and `"none"`.
#'        Default: `TRUE`
#'        
#' @param clust   Clustering algorithm for reordering the rows and columns by 
#'        similarity. You can supply a list or character vector of length two to 
#'        control the row and column clustering separately, for example 
#'        `clust = c(rows = "complete", cols = NA)`, or simply 
#'        `clust = c("complete", NA)`. Options are:
#'        \itemize{
#'          \item{`FALSE` or `NA` - }{ Disable reordering. }
#'          \item{An `hclust` class object}{ E.g. from [stats::hclust()]. }
#'          \item{A method name - }{ `"ward.D"`, 
#'            `"ward.D2"`, `"single"`, `"complete"`, 
#'            `"average"`, `"mcquitty"`, `"median"`, or 
#'            `"centroid"`. }
#'        }
#'        Default: `"complete"`
#'        
#' @param dist   Distance algorithm to use when reordering the rows and columns 
#'        by similarity. You can supply a list or character vector of length
#'        two to control the row and column clustering separately, for example 
#'        `dist = c(rows = "euclidean", cols = "maximum")`, or simply 
#'        `dist = c("euclidean", "maximum")`. Options are:
#'        \itemize{
#'          \item{A `dist` class object}{ E.g. from [stats::dist()] or [bdiv_distmat()]. }
#'          \item{A method name - }{ `"euclidean"`, 
#'            `"maximum"`, `"manhattan"`, `"canberra"`, 
#'            `"binary"`, or `"minkowski"`. }
#'        }
#'        Default: `"euclidean"`
#'        
#' @param tree_height,track_height   The height of the dendrogram or annotation
#'        tracks as a percentage of the overall grid size. Use a numeric vector 
#'        of length two to assign `c(top, left)` independently. 
#'        Default: `10` (10% of the grid's height)
#'        
#' @param asp   Aspect ratio (height/width) for entire grid.
#'        Default: `1` (square)
#'        
#' @param legend   Where to place the legend. Options are: `"right"` or
#'        `"bottom"`. Default: `"right"`
#'        
#' @param title   Plot title. Set to `TRUE` for a default title, `NULL` for 
#'        no title, or any character string. Default: `TRUE`
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
#'            \item{`"adonis2"` - }{ Permutational MANOVA; [vegan::adonis2()]. }
#'            \item{`"mrpp"` - }{ Multiple response permutation procedure; [vegan::mrpp()]. }
#'            \item{`"none"` - }{ Don't run any statistics. }
#'        }
#'        Default: `"adonis2"` \cr\cr
#'        Abbreviations are allowed.
#'        
NULL




# transform.y - ifelse ====
#' documentation_tranform.y.ifelse
#' 
#' @name documentation_tranform.y.ifelse
#' @keywords internal
#' 
#' @param transform.y   How to transform y-axis values. Useful for correcting 
#'        for non-normally distributions before applying regression statistics. 
#'        Options are `"none"` or `"rank"`. 
#'        Default: `ifelse(is.null(regr), "rank", "none")`
NULL



# test - ifelse ====
#' documentation_test.ifelse
#' 
#' @name documentation_test.ifelse
#' @keywords internal
#' 
#' 
#' @param test   The statistic of interest. An overview of options can be 
#'        found in the "Statistical Tests" section below. One of `"predict"`, 
#'        `"terms"`, `"fit"`, `"means"`, `"trends"`, `"es_means"`, 
#'        `"es_trends"`, `"pw_means"`, or `"pw_trends"`. 
#'        Default: `ifelse(is.null(regr), "means", "trends")`
#' 
#' 
#' @section Statistical Tests:
#' 
#' \bold{When `regr` is NULL:}
#' \itemize{
#'   \item{`"means"` - }{ 
#'     Considers all groups at once using the Kruskal-Wallis
#'     non-parametric test: [stats::kruskal.test()]. }
#'   \item{`"pw_means"` - }{
#'     Pairwise comparison of `stat.by` groups using the Wilcox rank sum 
#'     (aka Mann-Whitney) non-parametric test: [stats::wilcox.test()]. }
#' }
#' 
#' \bold{When `regr` is not NULL:}
#' \itemize{
#'   \item{`"predict"` - }{
#'     Augments original data with fitted information. See 
#'     [broom::augment.lm()]. }
#'   \item{`"terms"` - }{ Summary information about the model's 
#'     terms, including p-value, r-squared, AIC, BIC, and confidence 
#'     interval at `level`. See [broom::tidy.lm()]. }
#'   \item{`"fit"` - }{ Goodness of fit measures, p-values, and 
#'     more for the overall model. See [broom::glance.lm()]. }
#'   \item{`"means"` - }{ The estimated marginal mean (EMM) for 
#'     each `stat.by` group, along with confidence intervals (CI), 
#'     standard errors (SE), t-ratios, and p-values testing for mean = zero. 
#'     See [emmeans::emmeans()] and [emmeans::summary.emmGrid()]. }
#'   \item{`"trends"` - }{ The trendline slope for each 
#'     `stat.by` group, along with CI and SE and p-value testing for
#'     slope = zero. See [emmeans::emtrends()] and 
#'     [emmeans::summary.emmGrid()]. }
#'   \item{`"pw_means"` - }{ Pairwise means. All `stat.by` groups are 
#'     compared to each other and the difference in means is estimated 
#'     along with SE, t-ratios, and p-values testing if the two means 
#'     are the same. See [emmeans::pairs.emmGrid()]. }
#'   \item{`"pw_trends"` - }{ Pairwise trends. As above, but
#'     comparing trendline slopes instead of means. }
#'   \item{`"es_means"` - }{ Effect sizes for pairwise differences 
#'     of means, including SE and CI. See [emmeans::eff_size()]. }
#'   \item{`"es_trends"` - }{ Effect sizes for pairwise differences 
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
#'        `"predict"`, `"terms"`, `"fit"`, `"means"`, 
#'        `"trends"`, `"es_means"`, `"es_trends"`, 
#'        `"pw_means"`, or `"pw_trends"`. 
#'        Default: `"pw_means"`
#' 
#' 
#' @section Statistical Tests:
#' 
#' \itemize{
#'   \item{`"predict"` - }{
#'     Augments original data with fitted information. See 
#'     [broom::augment.lm()]. }
#'   \item{`"terms"` - }{ Summary information about the model's 
#'     terms, including p-value, r-squared, AIC, BIC, and confidence 
#'     interval at `level`. See [broom::tidy.lm()]. }
#'   \item{`"fit"` - }{ Goodness of fit measures, p-values, and 
#'     more for the overall model. See [broom::glance.lm()]. }
#'   \item{`"means"` - }{ The estimated marginal mean (EMM) for 
#'     each `stat.by` group, along with confidence intervals (CI), 
#'     standard errors (SE), t-ratios, and p-values testing for mean = zero. 
#'     See [emmeans::emmeans()] and [emmeans::summary.emmGrid()]. }
#'   \item{`"trends"` - }{ The trendline slope for each 
#'     `stat.by` group, along with CI and SE and p-value testing for
#'     slope = zero. See [emmeans::emtrends()] and 
#'     [emmeans::summary.emmGrid()]. }
#'   \item{`"pw_means"` - }{ Pairwise means. All `stat.by` groups are 
#'     compared to each other and the difference in means is estimated 
#'     along with SE, t-ratios, and p-values testing if the two means 
#'     are the same. See [emmeans::pairs.emmGrid()]. }
#'   \item{`"pw_trends"` - }{ Pairwise trends. As above, but
#'     comparing trendline slopes instead of means. }
#'   \item{`"es_means"` - }{ Effect sizes for pairwise differences 
#'     of means, including SE and CI. See [emmeans::eff_size()]. }
#'   \item{`"es_trends"` - }{ Effect sizes for pairwise differences 
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
#'        `"predict"`, `"terms"`, `"fit"`, `"means"`, 
#'        `"trends"`, `"es_means"`, `"es_trends"`, 
#'        `"pw_means"`, or `"pw_trends"`. 
#'        Default: `"trends"`
#' 
#' 
#' @section Statistical Tests:
#' 
#' \itemize{
#'   \item{`"predict"` - }{
#'     Augments original data with fitted information. See 
#'     [broom::augment.lm()]. }
#'   \item{`"terms"` - }{ Summary information about the model's 
#'     terms, including p-value, r-squared, AIC, BIC, and confidence 
#'     interval at `level`. See [broom::tidy.lm()]. }
#'   \item{`"fit"` - }{ Goodness of fit measures, p-values, and 
#'     more for the overall model. See [broom::glance.lm()]. }
#'   \item{`"means"` - }{ The estimated marginal mean (EMM) for 
#'     each `stat.by` group, along with confidence intervals (CI), 
#'     standard errors (SE), t-ratios, and p-values testing for mean = zero. 
#'     See [emmeans::emmeans()] and [emmeans::summary.emmGrid()]. }
#'   \item{`"trends"` - }{ The trendline slope for each 
#'     `stat.by` group, along with CI and SE and p-value testing for
#'     slope = zero. See [emmeans::emtrends()] and 
#'     [emmeans::summary.emmGrid()]. }
#'   \item{`"pw_means"` - }{ Pairwise means. All `stat.by` groups are 
#'     compared to each other and the difference in means is estimated 
#'     along with SE, t-ratios, and p-values testing if the two means 
#'     are the same. See [emmeans::pairs.emmGrid()]. }
#'   \item{`"pw_trends"` - }{ Pairwise trends. As above, but
#'     comparing trendline slopes instead of means. }
#'   \item{`"es_means"` - }{ Effect sizes for pairwise differences 
#'     of means, including SE and CI. See [emmeans::eff_size()]. }
#'   \item{`"es_trends"` - }{ Effect sizes for pairwise differences 
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
#'        `"lm"` (linear), `"log"` (logarithmic), or `"gam"` 
#'        (generalized additive). See the "Model Options" section below for
#'        additional details. Default: `"lm"`
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
#'        `"lm"` (linear), `"log"` (logarithmic), or `"gam"` 
#'        (generalized additive). See the "Model Options" section below for
#'        additional details. Default: `"log"`
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
#'   \item{`"lm"` - }{  Linear model: `stats::lm(formula = y ~ x)`.) }
#'   \item{`"log"` - }{ Logarithmic model: `stats::lm(formula = y ~ log(x))`. }
#'   \item{`"gam"` - }{ Generalized additive model: `mgcv::gam(formula = y ~ s(x, bs = "cs"), method = "REML")`. }
#' }
#' 
#' You can alternatively provide a list of length two where the first 
#' element is a character vector of length 1 naming a function, and the 
#' second element is a list of arguments to pass to that function. One 
#' of the function's arguments must be named 'formula'. 
#' For example, `model = list("stats::lm", list(formula = y ~ x))`.
#' 
#' By default, the same arguments are used for both plotting and calculating 
#' statistics. Arguments just for statistics can be set as a third list 
#' element: 
#' `model = list("stats::lm", list(formula = y ~ x), list(formula = rank(y) ~ x))`.
#' 
NULL




# plot - return ====
#' documentation_plot_return
#' 
#' @name documentation_plot_return
#' @keywords internal
#' 
#' @return A `ggplot2` plot. \cr The computed data points and ggplot 
#'         command are available as `$data` and `$code`, 
#'         respectively.
#' 
NULL


# plot w/ stats - return ====
#' documentation_plot_stats_return
#' 
#' @name documentation_plot_stats_return
#' @keywords internal
#' 
#' @return A `ggplot2` plot. \cr The computed data points, statistics, 
#'         and ggplot command are available as `$data`, `$stats`, and 
#'         `$code`, respectively.
#' 
NULL


# stats - return ====
#' documentation_stats_return
#' 
#' @name documentation_stats_return
#' @keywords internal
#' 
#' @return A tibble data frame with summary statistics. \cr
#'         The R code or generating these statistics is in `$code`.
#' 
NULL

