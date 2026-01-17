# Display taxa abundances as a heatmap.

Display taxa abundances as a heatmap.

## Usage

``` r
taxa_heatmap(
  biom,
  rank = -1,
  taxa = 6,
  tracks = NULL,
  grid = "bilbao",
  other = FALSE,
  unc = "singly",
  lineage = FALSE,
  label = TRUE,
  label_size = NULL,
  rescale = "none",
  trees = TRUE,
  clust = "complete",
  dist = "euclidean",
  asp = 1,
  tree_height = 10,
  track_height = 10,
  legend = "right",
  title = TRUE,
  xlab.angle = "auto",
  ...
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- rank:

  What rank(s) of taxa to display. E.g. `"Phylum"`, `"Genus"`, `".otu"`,
  etc. An integer vector can also be given, where `1` is the highest
  rank, `2` is the second highest, `-1` is the lowest rank, `-2` is the
  second lowest, and `0` is the OTU "rank". Run `biom$ranks` to see all
  options for a given rbiom object. Default: `-1`.

- taxa:

  Which taxa to display. An integer value will show the top n most
  abundant taxa. A value 0 \<= n \< 1 will show any taxa with that mean
  abundance or greater (e.g. `0.1` implies \>= 10%). A character vector
  of taxa names will show only those named taxa. Default: `6`.

- tracks:

  A character vector of metadata fields to display as tracks at the top
  of the plot. Or, a list as expected by the `tracks` argument of
  [`plot_heatmap()`](https://cmmr.github.io/rbiom/reference/plot_heatmap.md).
  Default: `NULL`

- grid:

  Color palette name, or a list as expected
  [`plot_heatmap()`](https://cmmr.github.io/rbiom/reference/plot_heatmap.md).
  Default: `"bilbao"`

- other:

  Sum all non-itemized taxa into an "Other" taxa. When `FALSE`, only
  returns taxa matched by the `taxa` argument. Specifying `TRUE` adds
  "Other" to the returned set. A string can also be given to imply
  `TRUE`, but with that value as the name to use instead of "Other".
  Default: `FALSE`

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

- lineage:

  Include all ranks in the name of the taxa. For instance, setting to
  `TRUE` will produce
  `Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales`.
  Otherwise the taxa name will simply be `Coriobacteriales`. You want to
  set this to TRUE when `unc = "asis"` and you have taxa names (such as
  *Incertae_Sedis*) that map to multiple higher level ranks. Default:
  `FALSE`

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

- trees:

  Draw a dendrogram for rows (left) and columns (top). You can supply a
  list or logical vector of length two to control the row tree and
  column tree separately, for example
  `trees = c(rows = TRUE, cols = FALSE)`, or simply
  `trees = c(TRUE, FALSE)`. Other valid options are `"rows"`, `"cols"`,
  `"both"`, `"left"`, `"top"`, and `"none"`. Default: `TRUE`

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

- dist:

  Distance algorithm to use when reordering the rows and columns by
  similarity. You can supply a list or character vector of length two to
  control the row and column clustering separately, for example
  `dist = c(rows = "euclidean", cols = "maximum")`, or simply
  `dist = c("euclidean", "maximum")`. Options are:

  A `dist` class object

  :   E.g. from [`stats::dist()`](https://rdrr.io/r/stats/dist.html) or
      [`bdiv_distmat()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md).

  A method name -

  :   `"euclidean"`, `"maximum"`, `"manhattan"`, `"canberra"`,
      `"binary"`, or `"minkowski"`.

  Default: `"euclidean"`

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

- ...:

  Additional arguments to pass on to ggplot2::theme().

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

Other taxa_abundance:
[`sample_sums()`](https://cmmr.github.io/rbiom/reference/sample_sums.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_clusters()`](https://cmmr.github.io/rbiom/reference/taxa_clusters.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md),
[`taxa_sums()`](https://cmmr.github.io/rbiom/reference/taxa_sums.md),
[`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)

Other visualization:
[`adiv_boxplot()`](https://cmmr.github.io/rbiom/reference/adiv_boxplot.md),
[`adiv_corrplot()`](https://cmmr.github.io/rbiom/reference/adiv_corrplot.md),
[`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md),
[`bdiv_corrplot()`](https://cmmr.github.io/rbiom/reference/bdiv_corrplot.md),
[`bdiv_heatmap()`](https://cmmr.github.io/rbiom/reference/bdiv_heatmap.md),
[`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md),
[`plot_heatmap()`](https://cmmr.github.io/rbiom/reference/plot_heatmap.md),
[`rare_corrplot()`](https://cmmr.github.io/rbiom/reference/rare_corrplot.md),
[`rare_multiplot()`](https://cmmr.github.io/rbiom/reference/rare_multiplot.md),
[`rare_stacked()`](https://cmmr.github.io/rbiom/reference/rare_stacked.md),
[`stats_boxplot()`](https://cmmr.github.io/rbiom/reference/stats_boxplot.md),
[`stats_corrplot()`](https://cmmr.github.io/rbiom/reference/stats_corrplot.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md)

## Examples

``` r
    library(rbiom)
    
    # Keep and rarefy the 10 most deeply sequenced samples.
    hmp10 <- rarefy(hmp50, n = 10)
    
    taxa_heatmap(hmp10, rank = "Phylum", tracks = "Body Site")

    
    taxa_heatmap(hmp10, rank = "Genus", tracks = c("sex", "bo"))

    
    taxa_heatmap(hmp10, rank = "Phylum", tracks = list(
      'Sex'       = list(colors = c(m = "#0000FF", f = "violetred")), 
      'Body Site' = list(colors = "muted", label = "Source") ))

    
```
