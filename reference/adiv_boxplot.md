# Visualize alpha diversity with boxplots.

Visualize alpha diversity with boxplots.

## Usage

``` r
adiv_boxplot(
  biom,
  x = NULL,
  adiv = "Shannon",
  layers = "x",
  stat.by = x,
  facet.by = NULL,
  colors = TRUE,
  shapes = TRUE,
  patterns = FALSE,
  test = "auto",
  flip = FALSE,
  stripe = NULL,
  ci = "ci",
  level = 0.95,
  p.adj = "fdr",
  outliers = NULL,
  xlab.angle = "auto",
  p.label = 0.05,
  transform = "none",
  strata = NULL,
  caption = TRUE,
  ...
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- x:

  A categorical metadata column name to use for the x-axis. Or `NULL`,
  which groups all samples into a single category.

- adiv:

  Alpha diversity metric(s) to use. Options are:
  `c("ace", "berger", "brillouin", "chao1", "faith", "fisher", "simpson", "inv_simpson", "margalef", "mcintosh", "menhinick", "observed", "shannon", "squares")`.
  For `"faith"`, a phylogenetic tree must be present in `biom` or
  explicitly provided via `tree=`. Set `adiv=".all"` to use all metrics.
  Multiple/abbreviated values allowed. Default: `"shannon"`

- layers:

  One or more of
  `c("bar", "box" ("x"), "violin", "dot", "strip", "crossbar", "errorbar", "linerange", "pointrange")`.
  Single letter abbreviations are also accepted. For instance,
  `c("box", "dot")` is equivalent to `c("x", "d")` and `"xd"`. Default:
  `"x"`

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

- patterns:

  Patterns for each group. Options are similar to `colors`'s: `TRUE`,
  `FALSE`, `NULL`, pattern names (`"brick"`, `"chevron"`, `"fish"`,
  `"grid"`, etc), or a named vector mapping groups to specific pattern
  names. See "Aesthetics" section below for additional information.
  Default: `FALSE`

- test:

  Method for computing p-values. `'auto'` runs pairwise Mann-Whitney U
  tests (`'wilcox'`) when `x` and `stat.by` are the same field (or
  `'emmeans'` when `strata` is given), and Kruskal-Wallis tests
  (`'kruskal'`) across the `stat.by` groups within each `x` category
  otherwise. Beta diversity boxplots use PERMANOVA (`'adonis2'`),
  `'mrpp'`, or `'betadisper'` on the samples instead. `'none'` skips the
  statistics. See
  [`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md)
  and
  [`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md).
  Default: `'auto'`

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
  for `"mad"` in which case it represents the median. Intervals are
  drawn as computed, including below zero, unless the y-axis is
  transformed: bounds below zero are drawn at zero on a `"sqrt"` or
  `"log1p"` axis, and bounds at or below zero are drawn at the smallest
  positive y value on a log axis. Default: `"ci"`

- level:

  The confidence level for calculating a confidence interval. Default:
  `0.95`

- p.adj:

  Method to use for multiple comparisons adjustment of p-values. Run
  `p.adjust.methods` for a list of available options. The adjustment
  spans every p-value in the returned table (all metrics, strata, pairs,
  and positions). Default: `"fdr"`

- outliers:

  Show boxplot outliers? `TRUE` to always show. `FALSE` to always hide.
  `NULL` to only hide them when overlaying a dot or strip chart.
  Default: `NULL`

- xlab.angle:

  Angle of the labels at the bottom of the plot. Options are `"auto"`,
  `'0'`, `'30'`, and `'90'`. Default: `"auto"`.

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

- transform:

  Transformation to apply to calculated values. Options are:
  `c("none", "rank", "log", "log1p", "sqrt")`. `"rank"` is useful for
  correcting for non-normal distributions before applying regression
  statistics. Default: `"none"`

- strata:

  Dataset field with the blocking/subject variable. Permutation tests
  permute within strata; regression models include a random intercept
  per stratum. Default: `NULL`

- caption:

  Add methodology caption beneath the plot. Default: `TRUE`

- ...:

  Additional parameters to pass along to ggplot2 functions. Prefix a
  parameter name with a layer name to pass it to only that layer. For
  instance, `d.size = 2` ensures only the points on the **dot** layer
  have their size set to `2`. Unprefixed parameters go to every layer
  whose function accepts them (e.g. `scales = "free"` for the facets).
  An argument that no layer accepts is an error.

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

Patterns are added using the fillpattern R package, which must be
installed to use `patterns`. Options are `"brick"`, `"chevron"`,
`"fish"`, `"grid"`, `"herringbone"`, `"hexagon"`, `"octagon"`, `"rain"`,
`"saw"`, `"shingle"`, `"rshingle"`, `"stripe"`, and `"wave"`, optionally
abbreviated and/or suffixed with modifiers. For example, `"hex10_sm"`
for the hexagon pattern rotated 10 degrees and shrunk by 2x. See
[`fillpattern::fill_pattern()`](https://cmmr.github.io/fillpattern/reference/fill_pattern.html)
for complete documentation of options.

Shapes can be given as per base R - numbers 0 through 17 for various
shapes, or the decimal value of an ascii character, e.g. a-z = 65:90;
A-Z = 97:122 to use letters instead of shapes on the plot. Character
strings may used as well.

The **dot** and **strip** layers use the ggbeeswarm package when it is
installed, and jittered points otherwise.

## See also

Other alpha_diversity:
[`adiv_corrplot()`](https://cmmr.github.io/rbiom/reference/adiv_corrplot.md),
[`adiv_stats()`](https://cmmr.github.io/rbiom/reference/adiv_stats.md),
[`adiv_table()`](https://cmmr.github.io/rbiom/reference/adiv_table.md)

Other visualization:
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
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md)

## Examples

``` r
    library(rbiom)

    biom <- rarefy(hmp50)

    adiv_boxplot(biom, x="Body Site", stat.by="Body Site")


    adiv_boxplot(biom, x="Sex", stat.by="Body Site", adiv=c("otu", "shan"), layers = "bld")


    adiv_boxplot(biom, x="body", stat.by="sex", adiv=c('sha', 'sim'), flip=TRUE, layers="p")
#> Warning: ! Statistics could not be computed for 2 of 10 comparisons: "Mid vagina /
#>   shannon" and "Mid vagina / simpson".
#> ℹ Message from the underlying tests:
#> • all observations are in the same group



    # Each plot object includes additional information.
    fig <- adiv_boxplot(biom, x="Body Site")

    ## Computed Data Points -------------------
    fig$data
#> # A tibble: 49 × 5
#>    .sample .depth .adiv   .diversity `Body Site`   
#>    <chr>    <dbl> <fct>        <dbl> <fct>         
#>  1 HMP01     1183 shannon       1.71 Buccal mucosa 
#>  2 HMP02     1183 shannon       2.58 Buccal mucosa 
#>  3 HMP03     1183 shannon       2.92 Saliva        
#>  4 HMP04     1183 shannon       3.26 Saliva        
#>  5 HMP05     1183 shannon       1.43 Buccal mucosa 
#>  6 HMP06     1183 shannon       3.05 Saliva        
#>  7 HMP07     1183 shannon       1.19 Buccal mucosa 
#>  8 HMP08     1183 shannon       2.50 Saliva        
#>  9 HMP09     1183 shannon       3.59 Saliva        
#> 10 HMP10     1183 shannon       1.75 Anterior nares
#> # ℹ 39 more rows

    ## Statistics Table -----------------------
    fig$stats
#> # Model:    wilcox.test(.resp ~ .stat.by)
#> # A tibble: 10 × 9
#>    `Body Site`        .loc.shift .h1    .p.val  .adj.p .lower .upper    .n .stat
#>    <chr>                   <dbl> <fct>   <dbl>   <dbl>  <dbl>  <dbl> <int> <dbl>
#>  1 Anterior nares - …     -1.51  != 0  1.08e-5 4.33e-5 -1.85  -1.30     20     0
#>  2 Mid vagina - Sali…     -2.71  != 0  1.08e-5 4.33e-5 -3.03  -2.43     20     0
#>  3 Anterior nares - …     -1.01  != 0  2.17e-5 4.33e-5 -1.32  -0.761    19     0
#>  4 Buccal mucosa - S…     -1.67  != 0  2.17e-5 4.33e-5 -2.20  -0.975    20     1
#>  5 Mid vagina - Stool     -2.19  != 0  2.17e-5 4.33e-5 -2.54  -1.90     19     0
#>  6 Buccal mucosa - S…     -1.11  != 0  9.74e-4 1.62e-3 -1.70  -0.466    19     7
#>  7 Saliva - Stool          0.481 != 0  1.45e-3 1.88e-3  0.210  0.847    19    82
#>  8 Anterior nares - …      1.20  != 0  1.50e-3 1.88e-3  0.841  1.64     20    90
#>  9 Buccal mucosa - M…      1.09  != 0  2.09e-3 2.32e-3  0.397  1.67     20    89
#> 10 Anterior nares - …      0.101 != 0  6.84e-1 6.84e-1 -0.439  0.595    20    56

    ## ggplot2 Command ------------------------
    fig$code
#> # Generated by rbiom 4.0.0
#> library(rbiom)
#> library(ggplot2)
#> 
#> biom <- hmp50
#> biom <- rarefy(biom, depth = 1183L)
#> 
#> data <- adiv_table(biom, "Shannon", "Body Site")
#> 
#> # Bracket segments connecting the compared groups.
#> brackets <- tibble::tibble(
#>   .x    = c(1, 3, 1, 2, 3, 2, 4, 1, 2, 1, 3, 1, 2, 3, 2, 4, 1, 2, 4, 4, 5, 4, 5, 5, 5, 3, 3),
#>   .xend = c(4, 4, 5, 4, 5, 5, 5, 3, 3, 1, 3, 1, 2, 3, 2, 4, 1, 2, 4, 4, 5, 4, 5, 5, 5, 3, 3),
#>   .y    = c(4, 4.52, 5.04, 5.56, 6.08, 6.6, 7.12, 7.64, 8.16, 4, 4.52, 5.04, 5.56, 6.08, 6.6, 7.12, 7.64, 8.16, 
#>     4, 4.52, 5.04, 5.56, 6.08, 6.6, 7.12, 7.64, 8.16),
#>   .yend = c(4, 4.52, 5.04, 5.56, 6.08, 6.6, 7.12, 7.64, 8.16, 3.87, 4.39, 4.91, 5.43, 5.95, 6.47, 6.99, 7.51, 8.03, 
#>     3.87, 4.39, 4.91, 5.43, 5.95, 6.47, 6.99, 7.51, 8.03) )
#> 
#> # p-values from the statistics table (see `p$stats$code`), placed above the boxes.
#> labels <- tibble::tibble(
#>   .label = c("italic(p)== 4e-05", "italic(p)== 4e-05", "italic(p)== 4e-05", "italic(p)== 4e-05", "italic(p)== 4e-05", 
#>     "italic(p)== 0.002", "italic(p)== 0.002", "italic(p)== 0.002", "italic(p)== 0.002"),
#>   .x     = c(2.5, 3.5, 3, 3, 4, 3.5, 4.5, 2, 2.5),
#>   .y     = c(4.0346666666666664, 4.554666666666666, 5.0746666666666664, 5.594666666666666, 6.1146666666666665, 6.634666666666666, 
#>     7.1546666666666665, 7.6746666666666661, 8.1946666666666665) )
#> 
#> ggplot(data, aes(x = `Body Site`, y = .diversity)) +
#>   annotate(
#>     geom  = "rect", 
#>     xmin  = -Inf, 
#>     xmax  = Inf, 
#>     ymin  = 4, 
#>     ymax  = Inf, 
#>     color = NA, 
#>     fill  = "white" ) +
#>   geom_boxplot(
#>     mapping      = aes(color = `Body Site`, fill = `Body Site`), 
#>     outlier.size = 1.14, 
#>     alpha        = 0.6, 
#>     width        = 0.7 ) +
#>   geom_segment(
#>     mapping = aes(x = .x, xend = .xend, y = .y, yend = .yend), 
#>     data    = brackets ) +
#>   geom_text(
#>     mapping = aes(x = .x, y = .y, label = .label), 
#>     data    = labels, 
#>     parse   = TRUE, 
#>     size    = 3, 
#>     vjust   = 0 ) +
#>   labs(
#>     caption = "Pairwise Mann-Whitney U tests; p-values adjusted with the FDR (Benjamini-Hochberg) method.", 
#>     y       = "Shannon Diversity Index" ) +
#>   scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")) +
#>   scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")) +
#>   scale_x_discrete() +
#>   scale_y_continuous(
#>     breaks       = c(0, 1, 2, 3), 
#>     minor_breaks = c(0.5, 1.5, 2.5, 3.5), 
#>     expand       = c(0.02, 0, 0.08, 0) ) +
#>   theme_bw() +
#>   theme(
#>     text               = element_text(size = 14), 
#>     axis.text.x        = element_text(hjust = 0, vjust = 1, angle = -30), 
#>     panel.grid.major.x = element_blank(), 
#>     plot.caption       = element_text(face = "italic", size = 9) )
```
