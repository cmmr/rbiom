# Visualize the number of observations per sample.

Visualize the number of observations per sample.

## Usage

``` r
rare_stacked(
  biom,
  rline = TRUE,
  counts = TRUE,
  labels = TRUE,
  y.transform = "log10",
  ...
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- rline:

  Where to draw a horizontal line on the plot, intended to show a
  particular rarefaction depth. Set to `TRUE` to show an auto-selected
  rarefaction depth, `FALSE` to not show a line, or an integer for a
  custom position. Default: `TRUE`.

- counts:

  Display the number of samples and reads remaining after rarefying to
  `rline` reads per sample. Default: `TRUE`.

- labels:

  Show sample names under each bar. Default: `TRUE`.

- y.transform:

  Y-axis transformation. Options are `"log10"` or `"none"`. Default:
  `"log10"`. Use `xaxis.transform` or `yaxis.transform` to pass custom
  values directly to ggplot2's `scale_*` functions.

- ...:

  Additional parameters to pass along to ggplot2 functions. Prefix a
  parameter name with `r.` to ensure it gets passed to (and only to)
  [geom_hline](https://ggplot2.tidyverse.org/reference/geom_abline.html).
  For instance, `r.color = "black"` ensures only the horizontal
  rarefaction line has its color set to `"black"`.

## Value

A `ggplot2` plot. The computed data points and ggplot command are
available as `$data` and `$code`, respectively.

## See also

Other rarefaction:
[`rare_corrplot()`](https://cmmr.github.io/rbiom/reference/rare_corrplot.md),
[`rare_multiplot()`](https://cmmr.github.io/rbiom/reference/rare_multiplot.md),
[`sample_sums()`](https://cmmr.github.io/rbiom/reference/sample_sums.md)

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
[`stats_boxplot()`](https://cmmr.github.io/rbiom/reference/stats_boxplot.md),
[`stats_corrplot()`](https://cmmr.github.io/rbiom/reference/stats_corrplot.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md)

## Examples

``` r
    library(rbiom) 
    
    rare_stacked(hmp50)

    
    rare_stacked(hmp50, rline = 500, r.linewidth = 2, r.linetype = "twodash")

    
    fig <- rare_stacked(hmp50, counts = FALSE)
    fig$code
#> ggplot(data) +
#>   geom_rect(
#>     mapping = aes(xmin = .xmin, xmax = .xmax, ymin = .ymin, ymax = .ymax, fill = .group), 
#>     color   = NA ) +
#>   geom_hline(
#>     yintercept = 1183, 
#>     color      = "red", 
#>     linetype   = "dashed" ) +
#>   labs(
#>     fill = "Reads", 
#>     x    = "Sample", 
#>     y    = "Sequencing Depth\n(log10 scale)" ) +
#>   scale_x_discrete() +
#>   scale_y_continuous(
#>     breaks       = 10^(0:5), 
#>     minor_breaks = as.vector(2:9 %o% 10^(0:4)), 
#>     labels       = scales::label_number(scale_cut = scales::cut_si("")), 
#>     expand       = c(0, 0), 
#>     transform    = "log10" ) +
#>   theme_bw() +
#>   theme(
#>     text               = element_text(size = 14), 
#>     panel.grid.major.x = element_blank() )
    
```
