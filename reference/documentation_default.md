# documentation_default

documentation_default

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- mtx:

  A matrix-like object.

- tree:

  A `phylo` object representing the phylogenetic relationships of the
  taxa in `biom`. Only required when computing UniFrac distances.
  Default: `biom$tree`

- underscores:

  When parsing the tree, should underscores be kept as is? By default
  they will be converted to spaces (unless the entire ID is quoted).
  Default `FALSE`

- md:

  Dataset field(s) to include in the output data frame, or `'.all'` to
  include all metadata fields. Default: `'.all'`

- adiv:

  Alpha diversity metric(s) to use. Options are:
  `c("ace", "berger", "brillouin", "chao1", "faith", "fisher", "simpson", "inv_simpson", "margalef", "mcintosh", "menhinick", "observed", "shannon", "squares")`.
  For `"faith"`, a phylogenetic tree must be present in `biom` or
  explicitly provided via `tree=`. Set `adiv=".all"` to use all metrics.
  Multiple/abbreviated values allowed. Default: `"shannon"`

- bdiv:

  Beta diversity distance algorithm(s) to use. Options are:
  `c("aitchison", "bhattacharyya", "bray", "canberra", "chebyshev", "chord", "clark", "sorensen", "divergence", "euclidean", "generalized_unifrac", "gower", "hamming", "hellinger", "horn", "jaccard", "jensen", "jsd", "lorentzian", "manhattan", "matusita", "minkowski", "morisita", "motyka", "normalized_unifrac", "ochiai", "psym_chisq", "soergel", "squared_chisq", "squared_chord", "squared_euclidean", "topsoe", "unweighted_unifrac", "variance_adjusted_unifrac", "wave_hedges", "weighted_unifrac")`.
  For the UniFrac family, a phylogenetic tree must be present in `biom`
  or explicitly provided via `tree=`. Supports partial matching.
  Multiple values are allowed for functions which return a table or
  plot. Default: `"bray"`

- taxa:

  Which taxa to display. An integer value will show the top n most
  abundant taxa. A value 0 \<= n \< 1 will show any taxa with that mean
  abundance or greater (e.g. `0.1` implies \>= 10%). A character vector
  of taxa names will show only those named taxa. Default: `6`.

- ord:

  Method for reducing dimensionality. Options are:

  `"PCoA"` -

  :   Principal coordinate analysis;
      [`ape::pcoa()`](https://rdrr.io/pkg/ape/man/pcoa.html).

  `"UMAP"` -

  :   Uniform manifold approximation and projection;
      [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).

  `"NMDS"` -

  :   Nonmetric multidimensional scaling;
      [`vegan::metaMDS()`](https://vegandevs.github.io/vegan/reference/metaMDS.html).

  `"tSNE"` -

  :   t-distributed stochastic neighbor embedding;
      [`tsne::tsne()`](https://rdrr.io/pkg/tsne/man/tsne.html).

  Multiple/abbreviated values allowed. Default: `"PCoA"`

- weighted:

  (Deprecated - weighting is now inherent in bdiv metric name.) Take
  relative abundances into account. When `weighted=FALSE`, only
  presence/absence is considered. Multiple values allowed. Default:
  `NULL`

- normalized:

  (Deprecated - normalization is now inherent in bdiv metric name.) Only
  changes the "Weighted UniFrac" calculation. Divides result by the
  total branch weights. Default: `NULL`

- delta:

  For numeric metadata, report the absolute difference in values for the
  two samples, for instance `2` instead of `"10 vs 12"`. Default: `TRUE`

- rank:

  What rank(s) of taxa to display. E.g. `"Phylum"`, `"Genus"`, `".otu"`,
  etc. An integer vector can also be given, where `1` is the highest
  rank, `2` is the second highest, `-1` is the lowest rank, `-2` is the
  second lowest, and `0` is the OTU "rank". Run `biom$ranks` to see all
  options for a given rbiom object. Default: `-1`.

- lineage:

  Include all ranks in the name of the taxa. For instance, setting to
  `TRUE` will produce
  `Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales`.
  Otherwise the taxa name will simply be `Coriobacteriales`. You want to
  set this to TRUE when `unc = "asis"` and you have taxa names (such as
  *Incertae_Sedis*) that map to multiple higher level ranks. Default:
  `FALSE`

- unc:

  How to handle unclassified, uncultured, and similarly ambiguous taxa
  names. Options are:

  `"singly"` -

  :   Replaces them with the OTU name.

  `"grouped"` -

  :   Replaces them with a higher rank's name.

  `"drop"` -

  :   Excludes them from the result.

  `"asis"` -

  :   To not check/modify any taxa names.

  Abbreviations are allowed. Default: `"singly"`

- other:

  Sum all non-itemized taxa into an "Other" taxa. When `FALSE`, only
  returns taxa matched by the `taxa` argument. Specifying `TRUE` adds
  "Other" to the returned set. A string can also be given to imply
  `TRUE`, but with that value as the name to use instead of "Other".
  Default: `FALSE`

- sparse:

  If `TRUE`, returns a sparse matrix from the `Matrix` package,
  otherwise returns a normal R matrix object. Default: `FALSE`

- p.top:

  Only display taxa with the most significant differences in abundance.
  If `p.top` is \>= 1, then the `p.top` most significant taxa are
  displayed. If `p.top` is less than one, all taxa with an adjusted
  p-value \<= `p.top` are displayed. Recommended to be used in
  combination with the `taxa` parameter to set a lower bound on the mean
  abundance of considered taxa. Default: `Inf`

- y.transform:

  The transformation to apply to the y-axis. Visualizing differences of
  both high- and low-abundance taxa is best done with a non-linear axis.
  Options are:

  `"sqrt"` -

  :   square-root transformation

  `"log1p"` -

  :   log(y + 1) transformation

  `"none"` -

  :   no transformation

  These methods allow visualization of both high- and low-abundance taxa
  simultaneously, without complaint about 'zero' count observations.
  Default: `"sqrt"` Use `xaxis.transform` or `yaxis.transform` to pass
  custom values directly to ggplot2's `scale_*` functions.

- flip:

  Transpose the axes, so that taxa are present as rows instead of
  columns. Default: `FALSE`

- stripe:

  Shade every other x position. Default: *same as flip*

- ci:

  How to calculate min/max of the **crossbar**, **errorbar**,
  **linerange**, and **pointrange** layers. Options are: `"ci"`
  (confidence interval), `"range"`, `"sd"` (standard deviation), `"se"`
  (standard error), and `"mad"` (median absolute deviation). The center
  mark of **crossbar** and **pointrange** represents the mean, except
  for `"mad"` in which case it represents the median. Default: `"ci"`

- p.label:

  Minimum adjusted p-value to display on the plot with a bracket.

  `p.label = 0.05` -

  :   Show p-values that are \<= 0.05.

  `p.label = 0` -

  :   Don't show any p-values on the plot.

  `p.label = 1` -

  :   Show all p-values on the plot.

  If a numeric vector with more than one value is provided, they will be
  used as breaks for asterisk notation. Default: `0.05`

- level:

  The confidence level for calculating a confidence interval. Default:
  `0.95`

- caption:

  Add methodology caption beneath the plot. Default: `TRUE`

- outliers:

  Show boxplot outliers? `TRUE` to always show. `FALSE` to always hide.
  `NULL` to only hide them when overlaying a dot or strip chart.
  Default: `NULL`

- xlab.angle:

  Angle of the labels at the bottom of the plot. Options are `"auto"`,
  `'0'`, `'30'`, and `'90'`. Default: `"auto"`.

- k:

  Number of ordination dimensions to return. Either `2L` or `3L`.
  Default: `2L`

- split.by:

  Dataset field(s) that the data should be split by prior to any
  calculations. Must be categorical. Default: `NULL`

- dm:

  A `dist`-class distance matrix, as returned from
  [`bdiv_distmat()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md)
  or [`stats::dist()`](https://rdrr.io/r/stats/dist.html). Required.

- groups:

  A named vector of grouping values. The names should correspond to
  `attr(dm, 'Labels')`. Values can be either categorical or numeric.
  Required.

- df:

  The dataset (data.frame or tibble object). "Dataset fields" mentioned
  below should match column names in `df`. Required.

- regr:

  Dataset field with the x-axis (independent; predictive) values. Must
  be numeric. Default: `NULL`

- resp:

  Dataset field with the y-axis (dependent; response) values, such as
  taxa abundance or alpha diversity. Default: `attr(df, 'response')`

- stat.by:

  Dataset field with the statistical groups. Must be categorical.
  Default: `NULL`

- color.by:

  Dataset field with the group to color by. Must be categorical.
  Default: `stat.by`

- shape.by:

  Dataset field with the group for shapes. Must be categorical. Default:
  `stat.by`

- facet.by:

  Dataset field(s) to use for faceting. Must be categorical. Default:
  `NULL`

- colors:

  How to color the groups. Options are:

  `TRUE` -

  :   Automatically select colorblind-friendly colors.

  `FALSE` or `NULL` -

  :   Don't use colors.

  a palette name -

  :   Auto-select colors from this set. E.g. `"okabe"`

  character vector -

  :   Custom colors to use. E.g. `c("red", "#00FF00")`

  named character vector -

  :   Explicit mapping. E.g. `c(Male = "blue", Female = "red")`

  See "Aesthetics" section below for additional information. Default:
  `TRUE`

- shapes:

  Shapes for each group. Options are similar to `colors`'s: `TRUE`,
  `FALSE`, `NULL`, shape names (typically integers 0 - 17), or a named
  vector mapping groups to specific shape names. See "Aesthetics"
  section below for additional information. Default: `TRUE`

- patterns:

  Patterns for each group. Options are similar to `colors`'s: `TRUE`,
  `FALSE`, `NULL`, pattern names (`"brick"`, `"chevron"`, `"fish"`,
  `"grid"`, etc), or a named vector mapping groups to specific pattern
  names. See "Aesthetics" section below for additional information.
  Default: `FALSE`

- test:

  Method for computing p-values: `'wilcox'`, `'kruskal'`, `'emmeans'`,
  or `'emtrends'`. Default: `'emmeans'`

- fit:

  How to fit the trendline. `'lm'`, `'log'`, or `'gam'`. Default:
  `'gam'`

- at:

  Position(s) along the x-axis where the means or slopes should be
  evaluated. Default: `NULL`, which samples 100 evenly spaced positions
  and selects the position where the p-value is most significant.

- alt:

  Alternative hypothesis direction. Options are `'!='` (two-sided; not
  equal to `mu`), `'<'` (less than `mu`), or `'>'` (greater than `mu`).
  Default: `'!='`

- mu:

  Reference value to test against. Default: `0`

- check:

  Generate additional plots to aid in assessing data normality. Default:
  `FALSE`

- within, between:

  Dataset field(s) for intra- or inter- sample comparisons.
  Alternatively, dataset field names given elsewhere can be prefixed
  with `'=='` or `'!='` to assign them to `within` or `between`,
  respectively. Default: `NULL`

- norm:

  Normalize the incoming counts. Options are:

  - `'none'`: No transformation.

  - `'percent'`: Relative abundance (sample abundances sum to 1).

  - `'binary'`: Unweighted presence/absence (each count is either 0 or
    1).

  - `'clr'`: Centered log ratio.

  Default: `'none'`.

- pairs:

  Which combinations of samples should distances be calculated for? The
  default value (`NULL`) calculates all-vs-all. Provide a numeric or
  logical vector specifying positions in the distance matrix to
  calculate. See examples.

- power:

  Scaling factor for the magnitude of differences between communities
  (\\p\\) when `bdiv = 'minkowski'`. Ignored for other beta diversity
  metrics. Default: `1.5`

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'`. Ignored for
  other normalization methods. Default: `NULL` (emits a warning).

- alpha:

  The alpha term to use in Generalized UniFrac. How much weight to give
  to relative abundances; a value between 0 and 1, inclusive. Setting
  `alpha=1` is equivalent to Normalized UniFrac. Default: `0.5`

- seed:

  Random seed for permutations. Must be a non-negative integer. Default:
  `0`

- cpus:

  The number of CPUs to use. Set to `NULL` to use all available, or to
  `1` to disable parallel processing. Default: `NULL`

- permutations:

  Number of random permutations to use. Default: `999`

- p.adj:

  Method to use for multiple comparisons adjustment of p-values. Run
  `p.adjust.methods` for a list of available options. Default: `"fdr"`

- depths:

  Rarefaction depths to show in the plot, or `NULL` to auto-select.
  Default: `NULL`

- rline:

  Where to draw a horizontal line on the plot, intended to show a
  particular rarefaction depth. Set to `TRUE` to show an auto-selected
  rarefaction depth or `FALSE` to not show a line. Default: `NULL`

- clone:

  Create a copy of `biom` before modifying. If `FALSE`, `biom` is
  modified in place as a side-effect. See [speed
  ups](https://cmmr.github.io/rbiom/reference/speed.md) for use cases.
  Default: `TRUE`

- labels:

  Show sample names under each bar. Default: `FALSE`

- transform:

  Transformation to apply to calculated values. Options are:
  `c("none", "rank", "log", "log1p", "sqrt", "percent")`. `"rank"` is
  useful for correcting for non-normally distributions before applying
  regression statistics. Default: `"none"`

- ties:

  When `transform="rank"`, how to rank identical values. Options are:
  `c("average", "first", "last", "random", "max", "min")`. See
  [`rank()`](https://rdrr.io/r/base/rank.html) for details. Default:
  `"random"`
