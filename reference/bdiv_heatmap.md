# Display beta diversities in an all vs all grid.

Display beta diversities in an all vs all grid.

## Usage

``` r
bdiv_heatmap(
  biom,
  bdiv = "bray",
  tree = NULL,
  tracks = NULL,
  grid = "devon",
  label = TRUE,
  label_size = NULL,
  rescale = "none",
  clust = "complete",
  trees = TRUE,
  asp = 1,
  tree_height = 10,
  track_height = 10,
  legend = "right",
  title = TRUE,
  xlab.angle = "auto",
  alpha = 0.5,
  cpus = n_cpus(),
  ...
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- bdiv:

  Beta diversity distance algorithm(s) to use. Options are:
  `c("aitchison", "bhattacharyya", "bray", "canberra", "chebyshev", "chord", "clark", "sorensen", "divergence", "euclidean", "generalized_unifrac", "gower", "hamming", "hellinger", "horn", "jaccard", "jensen", "jsd", "lorentzian", "manhattan", "matusita", "minkowski", "morisita", "motyka", "normalized_unifrac", "ochiai", "psym_chisq", "soergel", "squared_chisq", "squared_chord", "squared_euclidean", "topsoe", "unweighted_unifrac", "variance_adjusted_unifrac", "wave_hedges", "weighted_unifrac")`.
  For the UniFrac family, a phylogenetic tree must be present in `biom`
  or explicitly provided via `tree=`. Supports partial matching.
  Multiple values are allowed for functions which return a table or
  plot. Default: `"bray"`

- tree:

  A `phylo` object representing the phylogenetic relationships of the
  taxa in `biom`. Only required when computing UniFrac distances.
  Default: `biom$tree`

- tracks:

  A character vector of metadata fields to display as tracks at the top
  of the plot. Or, a list as expected by the `tracks` argument of
  [`plot_heatmap()`](https://cmmr.github.io/rbiom/reference/plot_heatmap.md).
  Default: `NULL`

- grid:

  Color palette name, or a list with entries for `label`, `colors`,
  `range`, `bins`, `na.color`, and/or `guide`. See the Track Definitions
  section for details. Default: `"devon"`

- label:

  Label the matrix rows and columns. You can supply a list or logical
  vector of length two to control row labels and column labels
  separately, for example `label = c(rows = TRUE, cols = FALSE)`, or
  simply `label = c(TRUE, FALSE)`. Other valid options are `"rows"`,
  `"cols"`, `"both"`, `"bottom"`, `"right"`, and `"none"`. Default:
  `TRUE`

- label_size:

  The font size to use for the row and column labels. You can supply a
  numeric vector of length two to control row label sizes and column
  label sizes separately, for example `c(rows = 20, cols = 8)`, or
  simply `c(20, 8)`. Default: `NULL`, which computes:
  `pmax(8, pmin(20, 100 / dim(mtx)))`

- rescale:

  Rescale rows or columns to all have a common min/max. Options:
  `"none"`, `"rows"`, or `"cols"`. Default: `"none"`

- clust:

  Clustering algorithm for reordering the rows and columns by
  similarity. You can supply a list or character vector of length two to
  control the row and column clustering separately, for example
  `clust = c(rows = "complete", cols = NA)`, or simply
  `clust = c("complete", NA)`. Options are:

  `FALSE` or `NA` -

  :   Disable reordering.

  An `hclust` class object

  :   E.g. from
      [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html).

  A method name -

  :   `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"`,
      `"mcquitty"`, `"median"`, or `"centroid"`.

  Default: `"complete"`

- trees:

  Draw a dendrogram for rows (left) and columns (top). You can supply a
  list or logical vector of length two to control the row tree and
  column tree separately, for example
  `trees = c(rows = TRUE, cols = FALSE)`, or simply
  `trees = c(TRUE, FALSE)`. Other valid options are `"rows"`, `"cols"`,
  `"both"`, `"left"`, `"top"`, and `"none"`. Default: `TRUE`

- asp:

  Aspect ratio (height/width) for entire grid. Default: `1` (square)

- tree_height, track_height:

  The height of the dendrogram or annotation tracks as a percentage of
  the overall grid size. Use a numeric vector of length two to assign
  `c(top, left)` independently. Default: `10` (10% of the grid's height)

- legend:

  Where to place the legend. Options are: `"right"` or `"bottom"`.
  Default: `"right"`

- title:

  Plot title. Set to `TRUE` for a default title, `NULL` for no title, or
  any character string. Default: `TRUE`

- xlab.angle:

  Angle of the labels at the bottom of the plot. Options are `"auto"`,
  `'0'`, `'30'`, and `'90'`. Default: `"auto"`.

- alpha:

  The alpha term to use in Generalized UniFrac. How much weight to give
  to relative abundances; a value between 0 and 1, inclusive. Setting
  `alpha=1` is equivalent to Normalized UniFrac. Default: `0.5`

- cpus:

  The number of CPUs to use. Set to `NULL` to use all available, or to
  `1` to disable parallel processing. Default: `NULL`

- ...:

  Additional arguments to pass on to ggplot2::theme(). For example,
  `labs.subtitle = "Plot subtitle"`.

## Value

A `ggplot2` plot. The computed data points and ggplot command are
available as `$data` and `$code`, respectively.

## Annotation Tracks

Metadata can be displayed as colored tracks above the heatmap. Common
use cases are provided below, with more thorough documentation available
at https://cmmr.github.io/rbiom .

    ## Categorical ----------------------------
    tracks = "Body Site"
    tracks = list('Body Site' = "bright")
    tracks = list('Body Site' = c('Stool' = "blue", 'Saliva' = "green"))

    ## Numeric --------------------------------
    tracks = "Age"
    tracks = list('Age' = "reds")

    ## Multiple Tracks ------------------------
    tracks = c("Body Site", "Age")
    tracks = list('Body Site' = "bright", 'Age' = "reds")
    tracks = list(
      'Body Site' = c('Stool' = "blue", 'Saliva' = "green"),
      'Age'       = list('colors' = "reds") )

The following entries in the track definitions are understood:

- `colors` - :

  A pre-defined palette name or custom set of colors to map to.

- `range` - :

  The c(min,max) to use for scale values.

- `label` - :

  Label for this track. Defaults to the name of this list element.

- `side` - :

  Options are `"top"` (default) or `"left"`.

- `na.color` - :

  The color to use for `NA` values.

- `bins` - :

  Bin a gradient into this many bins/steps.

- `guide` - :

  A list of arguments for guide_colorbar() or guide_legend().

All built-in color palettes are colorblind-friendly.

Categorical palette names: `"okabe"`, `"carto"`, `"r4"`, `"polychrome"`,
`"tol"`, `"bright"`, `"light"`, `"muted"`, `"vibrant"`, `"tableau"`,
`"classic"`, `"alphabet"`, `"tableau20"`, `"kelly"`, and `"fishy"`.

Numeric palette names: `"reds"`, `"oranges"`, `"greens"`, `"purples"`,
`"grays"`, `"acton"`, `"bamako"`, `"batlow"`, `"bilbao"`, `"buda"`,
`"davos"`, `"devon"`, `"grayC"`, `"hawaii"`, `"imola"`, `"lajolla"`,
`"lapaz"`, `"nuuk"`, `"oslo"`, `"tokyo"`, `"turku"`, `"bam"`,
`"berlin"`, `"broc"`, `"cork"`, `"lisbon"`, `"roma"`, `"tofino"`,
`"vanimo"`, and `"vik"`.

## See also

Other beta_diversity:
[`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md),
[`bdiv_clusters()`](https://cmmr.github.io/rbiom/reference/bdiv_clusters.md),
[`bdiv_corrplot()`](https://cmmr.github.io/rbiom/reference/bdiv_corrplot.md),
[`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md),
[`bdiv_ord_table()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_table.md),
[`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md),
[`bdiv_table()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md),
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md)

Other visualization:
[`adiv_boxplot()`](https://cmmr.github.io/rbiom/reference/adiv_boxplot.md),
[`adiv_corrplot()`](https://cmmr.github.io/rbiom/reference/adiv_corrplot.md),
[`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md),
[`bdiv_corrplot()`](https://cmmr.github.io/rbiom/reference/bdiv_corrplot.md),
[`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md),
[`plot_heatmap()`](https://cmmr.github.io/rbiom/reference/plot_heatmap.md),
[`rare_corrplot()`](https://cmmr.github.io/rbiom/reference/rare_corrplot.md),
[`rare_multiplot()`](https://cmmr.github.io/rbiom/reference/rare_multiplot.md),
[`rare_stacked()`](https://cmmr.github.io/rbiom/reference/rare_stacked.md),
[`stats_boxplot()`](https://cmmr.github.io/rbiom/reference/stats_boxplot.md),
[`stats_corrplot()`](https://cmmr.github.io/rbiom/reference/stats_corrplot.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md)

## Examples

``` r
    library(rbiom)
    
    # Subset to 10 samples and rarefy them.
    hmp10 <- rarefy(hmp50[1:10])
    
    bdiv_heatmap(hmp10, tracks=c("Body Site", "Age"))

    
    bdiv_heatmap(hmp10, bdiv=c("u_unifrac", "n_unifrac"), tracks="sex")
```
