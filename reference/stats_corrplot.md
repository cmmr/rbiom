# Visualize regression with scatterplots and trendlines.

Visualize regression with scatterplots and trendlines.

## Usage

``` r
stats_corrplot(
  df,
  x,
  y = attr(df, "response"),
  layers = "tc",
  stat.by = NULL,
  facet.by = NULL,
  colors = TRUE,
  shapes = TRUE,
  test = "emmeans",
  fit = "gam",
  at = NULL,
  level = 0.95,
  p.adj = "fdr",
  p.top = Inf,
  alt = "!=",
  mu = 0,
  caption = TRUE,
  check = FALSE,
  ...
)
```

## Arguments

- df:

  The dataset (data.frame or tibble object). "Dataset fields" mentioned
  below should match column names in `df`. Required.

- x:

  Dataset field with the x-axis values. Equivalent to the `regr`
  argument in
  [`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md).
  Required.

- y:

  A numeric metadata column name to use for the y-axis. Default:
  `attr(df, 'response')`

- layers:

  One or more of
  `c("trend", "confidence", "point", "name", "residual")`. Single letter
  abbreviations are also accepted. For instance, `c("trend", "point")`
  is equivalent to `c("t", "p")` and `"tp"`. Default: `"tc"`

- stat.by:

  Dataset field with the statistical groups. Must be categorical.
  Default: `NULL`

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

- test:

  Method for computing p-values: `'none'`, `'emmeans'`, or `'emtrends'`.
  Default: `'emmeans'`

- fit:

  How to fit the trendline. `'lm'`, `'log'`, or `'gam'`. Default:
  `'gam'`

- at:

  Position(s) along the x-axis where the means or slopes should be
  evaluated. Default: `NULL`, which samples 100 evenly spaced positions
  and selects the position where the p-value is most significant.

- level:

  The confidence level for calculating a confidence interval. Default:
  `0.95`

- p.adj:

  Method to use for multiple comparisons adjustment of p-values. Run
  `p.adjust.methods` for a list of available options. Default: `"fdr"`

- p.top:

  Only display taxa with the most significant differences in abundance.
  If `p.top` is \>= 1, then the `p.top` most significant taxa are
  displayed. If `p.top` is less than one, all taxa with an adjusted
  p-value \<= `p.top` are displayed. Recommended to be used in
  combination with the `taxa` parameter to set a lower bound on the mean
  abundance of considered taxa. Default: `Inf`

- alt:

  Alternative hypothesis direction. Options are `'!='` (two-sided; not
  equal to `mu`), `'<'` (less than `mu`), or `'>'` (greater than `mu`).
  Default: `'!='`

- mu:

  Reference value to test against. Default: `0`

- caption:

  Add methodology caption beneath the plot. Default: `TRUE`

- check:

  Generate additional plots to aid in assessing data normality. Default:
  `FALSE`

- ...:

  Additional parameters to pass along to ggplot2 functions. Prefix a
  parameter name with a layer name to pass it to only that layer. For
  instance, `p.size = 2` ensures only the points have their size set to
  `2`.

## Value

A `ggplot2` plot. The computed data points, ggplot2 command, stats
table, and stats table commands are available as `$data`, `$code`,
`$stats`, and `$stats$code`, respectively.

## Aesthetics

All built-in color palettes are colorblind-friendly. The available
categorical palette names are: `"okabe"`, `"carto"`, `"r4"`,
`"polychrome"`, `"tol"`, `"bright"`, `"light"`, `"muted"`, `"vibrant"`,
`"tableau"`, `"classic"`, `"alphabet"`, `"tableau20"`, `"kelly"`, and
`"fishy"`.

Shapes can be given as per base R - numbers 0 through 17 for various
shapes, or the decimal value of an ascii character, e.g. a-z = 65:90;
A-Z = 97:122 to use letters instead of shapes on the plot. Character
strings may used as well.

## See also

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
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md)

## Examples

``` r
    library(rbiom)
    
    biom <- subset(hmp50, `Body Site` %in% c('Saliva', 'Stool'))
    df   <- adiv_table(rarefy(biom))
    stats_corrplot(df, "age", stat.by = "body")

    stats_corrplot(
      df       = df, 
      x        = "Age", 
      stat.by  = "Body Site", 
      facet.by = "Sex", 
      layers   = "trend" )
```
