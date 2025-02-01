
# default ====
#' documentation_default
#' 
#' @name documentation_default
#' @keywords internal
#' 
#' @param biom   An [rbiom object][rbiom_objects], such as from [as_rbiom()]. 
#'        Any value accepted by [as_rbiom()] can also be given here.
#' 
#' @param mtx   A matrix-like object.
#'        
#' @param tree  A `phylo` object representing the phylogenetic
#'        relationships of the taxa in `biom`. Only required when 
#'        computing UniFrac distances. Default: `biom$tree`
#'
#' @param underscores   When parsing the tree, should underscores be kept as 
#'        is? By default they will be converted to spaces (unless the entire ID
#'        is quoted). Default `FALSE`
#'     
#' @param md  Dataset field(s) to include in the output data frame, or `'.all'` 
#'        to include all metadata fields. Default: `'.all'`
#' 
#' @param adiv   Alpha diversity metric(s) to use. Options are: `"OTUs"`, 
#'        `"Shannon"`, `"Chao1"`, `"Simpson"`, and/or 
#'        `"InvSimpson"`. Set `adiv=".all"` to use all metrics.
#'        Multiple/abbreviated values allowed.
#'        Default: `"Shannon"`
#' 
#' @param bdiv  Beta diversity distance algorithm(s) to use. Options are:
#'        `"Bray-Curtis"`, `"Manhattan"`, `"Euclidean"`, 
#'        `"Jaccard"`, and `"UniFrac"`. For `"UniFrac"`, a 
#'        phylogenetic tree must be present in `biom` or explicitly 
#'        provided via `tree=`. Multiple/abbreviated values allowed. 
#'        Default: `"Bray-Curtis"`
#'       
#' 
#' @param taxa   Which taxa to display. An integer value will show the top n
#'        most abundant taxa. A value 0 <= n < 1 will show any taxa with that 
#'        mean abundance or greater (e.g. `0.1` implies >= 10%). A 
#'        character vector of taxa names will show only those named taxa. 
#'        Default: `6`.
#'        
#' @param ord    Method for reducing dimensionality. Options are:
#'        \describe{
#'            \item{`"PCoA"` - }{ Principal coordinate analysis; [ape::pcoa()]. }
#'            \item{`"UMAP"` - }{ Uniform manifold approximation and projection; [uwot::umap()]. }
#'            \item{`"NMDS"` - }{ Nonmetric multidimensional scaling; [vegan::metaMDS()]. }
#'            \item{`"tSNE"` - }{ t-distributed stochastic neighbor embedding; [tsne::tsne()]. }
#'        }
#'        Multiple/abbreviated values allowed. Default: `"PCoA"`
#'        
#'     
#' @param weighted  Take relative abundances into account. When 
#'        `weighted=FALSE`, only presence/absence is considered.
#'        Multiple values allowed. Default: `TRUE`
#'     
#' @param normalized  Only changes the "Weighted UniFrac" calculation.
#'        Divides result by the total branch weights. Default: `TRUE`
#'       
#' 
#' @param delta  For numeric metadata, report the absolute difference in values 
#'        for the two samples, for instance `2` instead of `"10 vs 12"`. 
#'        Default: `TRUE`
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
#'        \describe{
#'          \item{`"singly"` - }{ Replaces them with the OTU name. }
#'          \item{`"grouped"` - }{ Replaces them with a higher rank's name. }
#'          \item{`"drop"` - }{ Excludes them from the result. }
#'          \item{`"asis"` - }{ To not check/modify any taxa names. }
#'        }
#'        Abbreviations are allowed. Default: `"singly"`
#'       
#'        
#' @param other  Sum all non-itemized taxa into an "Other" taxa. When 
#'        `FALSE`, only returns taxa matched by the `taxa` 
#'        argument. Specifying `TRUE` adds "Other" to the returned set.
#'        A string can also be given to imply `TRUE`, but with that
#'        value as the name to use instead of "Other".
#'        Default: `FALSE`
#'        
#' @param sparse  If `TRUE`, returns a 
#'        [slam::simple_triplet_matrix()], otherwise returns a 
#'        normal R matrix object. Default: `FALSE`
#'        
#' @param p.top   Only display taxa with the most significant differences in 
#'        abundance. If `p.top` is >= 1, then the `p.top` most 
#'        significant taxa are displayed. If `p.top` is less than one, all 
#'        taxa with an adjusted p-value <= `p.top` are displayed. 
#'        Recommended to be used in combination with the `taxa` parameter 
#'        to set a lower bound on the mean abundance of considered taxa. 
#'        Default: `Inf`
#'        
#' @param y.transform   The transformation to apply to the y-axis. Visualizing 
#'        differences of both high- and low-abundance taxa is best done with
#'        a non-linear axis. Options are: 
#'        \describe{
#'          \item{`"sqrt"` - }{ square-root transformation }
#'          \item{`"log1p"` - }{ log(y + 1) transformation }
#'          \item{`"none"` - }{ no transformation }
#'        }
#'        These methods allow visualization of both high- and low-abundance
#'        taxa simultaneously, without complaint about 'zero' count
#'        observations. Default: `"sqrt"`
#'        Use `xaxis.transform` or `yaxis.transform` to pass custom values 
#'        directly to ggplot2's `scale_*` functions.
#'        
#' @param flip   Transpose the axes, so that taxa are present as rows instead
#'        of columns. Default: `FALSE`
#'
#' @param stripe   Shade every other x position. Default: \emph{same as flip}
#'     
#' @param ci   How to calculate min/max of the \bold{crossbar}, 
#'        \bold{errorbar}, \bold{linerange}, and \bold{pointrange} layers.
#'        Options are: `"ci"` (confidence interval), `"range"`, 
#'        `"sd"` (standard deviation), `"se"` (standard error), and 
#'        `"mad"` (median absolute deviation). 
#'        The center mark of \bold{crossbar} and \bold{pointrange} represents
#'        the mean, except for `"mad"` in which case it represents the median. 
#'        Default: `"ci"`
#'
#' @param p.label   Minimum adjusted p-value to display on the plot with a 
#'        bracket.
#'        \describe{
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
#' @param split.by   Dataset field(s) that the data should be split by prior to 
#'        any calculations. Must be categorical. Default: `NULL`
#' 
#' @param dm   A `dist`-class distance matrix, as returned from 
#'        [bdiv_distmat()] or [stats::dist()]. Required.
#'        
#' @param groups  A named vector of grouping values. The names should 
#'        correspond to `attr(dm, 'Labels')`. Values can be either 
#'        categorical or numeric. Required.
#' 
#' @param df   The dataset (data.frame or tibble object). "Dataset fields" 
#'        mentioned below should match column names in `df`. Required.
#'        
#' @param regr   Dataset field with the x-axis (independent; predictive) 
#'        values. Must be numeric. Default: `NULL`
#' 
#' @param resp   Dataset field with the y-axis (dependent; response) values, 
#'        such as taxa abundance or alpha diversity. 
#'        Default: `attr(df, 'response')`
#'        
#' @param stat.by   Dataset field with the statistical groups. Must be 
#'        categorical. Default: `NULL`
#'        
#' @param color.by   Dataset field with the group to color by. Must be 
#'        categorical. Default: `stat.by`
#'        
#' @param shape.by   Dataset field with the group for shapes. Must be 
#'        categorical. Default: `stat.by`
#'        
#' @param facet.by   Dataset field(s) to use for faceting. Must be categorical. 
#'        Default: `NULL`
#'        
#' @param colors   How to color the groups. Options are:
#'        \describe{
#'            \item{`TRUE` - }{ Automatically select colorblind-friendly colors. }
#'            \item{`FALSE` or `NULL` - }{ Don't use colors. }
#'            \item{a palette name - }{ Auto-select colors from this set. E.g. `"okabe"` }
#'            \item{character vector - }{ Custom colors to use. E.g. `c("red", "#00FF00")` }
#'            \item{named character vector - }{ Explicit mapping. E.g. `c(Male = "blue", Female = "red")` }
#'        }
#'        See "Aesthetics" section below for additional information.
#'        Default: `TRUE`
#'        
#' @param shapes   Shapes for each group. 
#'        Options are similar to `colors`'s: `TRUE`, `FALSE`, `NULL`, shape 
#'        names (typically integers 0 - 17), or a named vector mapping 
#'        groups to specific shape names.
#'        See "Aesthetics" section below for additional information.
#'        Default: `TRUE`
#'        
#' @param patterns   Patterns for each group. 
#'        Options are similar to `colors`'s: `TRUE`, `FALSE`, `NULL`, pattern 
#'        names (`"brick"`, `"chevron"`, `"fish"`, `"grid"`, etc), or a named 
#'        vector mapping groups to specific pattern names.
#'        See "Aesthetics" section below for additional information.
#'        Default: `FALSE`
#' 
#' @param test   Method for computing p-values: `'wilcox'`, `'kruskal'`, 
#'        `'emmeans'`, or `'emtrends'`. Default: `'emmeans'`
#' 
#' @param fit   How to fit the trendline. `'lm'`, `'log'`, or `'gam'`. 
#'        Default: `'gam'`
#'        
#' @param at   Position(s) along the x-axis where the means or slopes should be 
#'        evaluated. Default: `NULL`, which samples 100 evenly spaced positions 
#'        and selects the position where the p-value is most significant.
#'        
#' @param alt   Alternative hypothesis direction. Options are `'!='` 
#'        (two-sided; not equal to `mu`), `'<'` (less than `mu`), or `'>'` 
#'        (greater than `mu`). Default: `'!='`
#'        
#' @param mu   Reference value to test against. Default: `0`
#' 
#' @param check   Generate additional plots to aid in assessing data normality. 
#'        Default: `FALSE`
#'        
#' @param within,between   Dataset field(s) for intra- or inter- sample 
#'        comparisons. Alternatively, dataset field names given elsewhere can 
#'        be prefixed with `'=='` or `'!='` to assign them to `within` or 
#'        `between`, respectively. Default: `NULL`
#'        
#' @param seed  Random seed for permutations. Must be a non-negative integer. 
#'              Default: `0`
#'        
#' @param cpus  The number of CPUs to use. Set to `NULL` to use all available, 
#'        or to `1` to disable parallel processing. Default: `NULL`
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
#' @param transform   Transformation to apply. Options are: 
#'        `c("none", "rank", "log", "log1p", "sqrt", "percent")`. `"rank"` is 
#'        useful for correcting for non-normally distributions before applying 
#'        regression statistics. Default: `"none"`
#' 
#' @param ties   When `transform="rank"`, how to rank identical values.
#'        Options are: `c("average", "first", "last", "random", "max", "min")`. 
#'        See `rank()` for details. Default: `"random"`
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
#' @param .data    An [rbiom object][rbiom_objects], such as from [as_rbiom()].
#' 
#' @param x    An [rbiom object][rbiom_objects], such as from [as_rbiom()].
#' 
#' @param object    An [rbiom object][rbiom_objects], such as from [as_rbiom()].
#' 
#' @param data    An [rbiom object][rbiom_objects], such as from [as_rbiom()].
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



# rank - 2 ====
#' documentation_rank.2
#' 
#' @name documentation_rank.2
#' @keywords internal
#'        
#' @param rank   What rank(s) of taxa to compute biplot coordinates and 
#'        statistics for, or `NULL` to disable. E.g. `"Phylum"`, 
#'        `"Genus"`, `".otu"`, etc. An integer vector can also be 
#'        given, where `1` is the highest rank, `2` is the second 
#'        highest, `-1` is the lowest rank, `-2` is the second 
#'        lowest, and `0` is the OTU "rank". Run `biom$ranks` to 
#'        see all options for a given rbiom object. Default: `2`.
#' 
NULL



# clusters ====
#' documentation_clusters
#' 
#' @name documentation_clusters
#' @keywords internal
#'        
#' @param k   Number of clusters. Default: `5L`
#'        
#' @param rank   Which taxa rank to use. E.g. `"Phylum"`, 
#'        `"Genus"`, `".otu"`, etc. An integer can also be 
#'        given, where `1` is the highest rank, `2` is the second 
#'        highest, `-1` is the lowest rank, `-2` is the second 
#'        lowest, and `0` is the OTU "rank". Run `biom$ranks` to 
#'        see all options for a given rbiom object. Default: `.otu`.
#' 
#' @return A numeric factor assigning samples to clusters.
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
#' @section Metadata Comparisons:
#' 
#' Prefix metadata fields with `==` or `!=` to limit comparisons to within or 
#' between groups, respectively. For example, `stat.by = '==Sex'` will 
#' run calculations only for intra-group comparisons, returning "Male" and
#' "Female", but NOT "Female vs Male". Similarly, setting 
#' `stat.by = '!=Body Site'` will only show the inter-group comparisons, such 
#' as "Saliva vs Stool", "Anterior nares vs Buccal mucosa", and so on.
#' 
#' The same effect can be achieved by using the `within` and `between` 
#' parameters. `stat.by = '==Sex'` is equivalent to 
#' `stat.by = 'Sex', within = 'Sex'`.
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
#'        `label = c(TRUE, FALSE)`. Other valid options are `"rows"`,
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
#'        `trees = c(rows = TRUE, cols = FALSE)`, 
#'        or simply `trees = c(TRUE, FALSE)`. 
#'        Other valid options are `"rows"`, `"cols"`, `"both"`, 
#'        `"left"`, `"top"`, and `"none"`.
#'        Default: `TRUE`
#'        
#' @param clust   Clustering algorithm for reordering the rows and columns by 
#'        similarity. You can supply a list or character vector of length two to 
#'        control the row and column clustering separately, for example 
#'        `clust = c(rows = "complete", cols = NA)`, or simply 
#'        `clust = c("complete", NA)`. Options are:
#'        \describe{
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
#'        \describe{
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
#'        \describe{
#'            \item{`"adonis2"` - }{ Permutational MANOVA; [vegan::adonis2()]. }
#'            \item{`"mrpp"` - }{ Multiple response permutation procedure; [vegan::mrpp()]. }
#'            \item{`"none"` - }{ Don't run any statistics. }
#'        }
#'        Abbreviations are allowed. Default: `"adonis2"`
#'       
#'        
NULL




# plot - return ====
#' documentation_plot_return
#' 
#' @name documentation_plot_return
#' @keywords internal
#' 
#' @return A `ggplot2` plot. The computed data points and ggplot 
#'         command are available as `$data` and `$code`, 
#'         respectively.
#' 
NULL

