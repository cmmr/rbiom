# documentation_heatmap

documentation_heatmap

## Arguments

- grid:

  Color palette name, or a list with entries for `label`, `colors`,
  `range`, `bins`, `na.color`, and/or `guide`. See the Track Definitions
  section for details. Default:
  `list(label = "Grid Value", colors = "imola")`

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

- tree_height, track_height:

  The height of the dendrogram or annotation tracks as a percentage of
  the overall grid size. Use a numeric vector of length two to assign
  `c(top, left)` independently. Default: `10` (10% of the grid's height)

- asp:

  Aspect ratio (height/width) for entire grid. Default: `1` (square)

- legend:

  Where to place the legend. Options are: `"right"` or `"bottom"`.
  Default: `"right"`

- title:

  Plot title. Set to `TRUE` for a default title, `NULL` for no title, or
  any character string. Default: `TRUE`

- ...:

  Additional arguments to pass on to ggplot2::theme().

## Value

A `ggplot2` plot. The computed data points and ggplot command are
available as `$data` and `$code`, respectively.
