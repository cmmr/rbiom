# Run statistics on a distance matrix vs a categorical or numeric variable.

Run statistics on a distance matrix vs a categorical or numeric
variable.

## Usage

``` r
distmat_stats(dm, groups, test = "adonis2", seed = 0, permutations = 999)
```

## Arguments

- dm:

  A `dist`-class distance matrix, as returned from
  [`bdiv_distmat()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md)
  or [`stats::dist()`](https://rdrr.io/r/stats/dist.html). Required.

- groups:

  A named vector of grouping values. The names should correspond to
  `attr(dm, 'Labels')`. Values can be either categorical or numeric.
  Required.

- test:

  Permutational test for accessing significance. Options are:

  `"adonis2"` -

  :   Permutational MANOVA;
      [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html).

  `"mrpp"` -

  :   Multiple response permutation procedure;
      [`vegan::mrpp()`](https://vegandevs.github.io/vegan/reference/mrpp.html).

  `"none"` -

  :   Don't run any statistics.

  Abbreviations are allowed. Default: `"adonis2"`

- seed:

  Random seed for permutations. Must be a non-negative integer. Default:
  `0`

- permutations:

  Number of random permutations to use. Default: `999`

## Value

A data.frame with summary statistics from
[`vegan::permustats()`](https://vegandevs.github.io/vegan/reference/permustats.html).
The columns are:

- *.n* - :

  The size of the distance matrix.

- *.stat* - :

  The observed statistic. For mrpp, this is the overall weighted mean of
  group mean distances.

- *.z* - :

  The difference of observed statistic and mean of permutations divided
  by the standard deviation of permutations (also known as z-values).
  Evaluated from permuted values without observed statistic.

- *.p.val* - :

  Probability calculated by `test`.

R commands for reproducing the results are in `$code`.

## See also

Other beta_diversity:
[`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md),
[`bdiv_clusters()`](https://cmmr.github.io/rbiom/reference/bdiv_clusters.md),
[`bdiv_corrplot()`](https://cmmr.github.io/rbiom/reference/bdiv_corrplot.md),
[`bdiv_heatmap()`](https://cmmr.github.io/rbiom/reference/bdiv_heatmap.md),
[`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md),
[`bdiv_ord_table()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_table.md),
[`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md),
[`bdiv_table()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md)

Other stats_tables:
[`adiv_stats()`](https://cmmr.github.io/rbiom/reference/adiv_stats.md),
[`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md),
[`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md)

## Examples

``` r
    library(rbiom)
    
    hmp10        <- hmp50$clone()
    hmp10$counts <- hmp10$counts[,1:10]
    
    dm <- bdiv_distmat(hmp10, 'w_unifrac')
    
    distmat_stats(dm, groups = pull(hmp10, 'Body Site'))
#> # A tibble: 1 × 4
#>      .n .stat    .z .p.val
#>   <int> <dbl> <dbl>  <dbl>
#> 1    10  9.27  11.0  0.001
    
    distmat_stats(dm, groups = pull(hmp10, 'Age'))
#> # A tibble: 1 × 4
#>      .n .stat    .z .p.val
#>   <int> <dbl> <dbl>  <dbl>
#> 1    10  4.13  3.56  0.012
    
    # See the R code used to calculate these statistics:
    stats <- distmat_stats(dm, groups = pull(hmp10, 'Age'))
    stats$code
#> grouping <- groups[attr(dm, 'Labels')]
#> set.seed(0)
#> 
#> vegan::adonis2(formula = dm ~ grouping, permutations = 999) %>%
#>   vegan::permustats() %>%
#>   summary() %>%
#>   with(data.frame(.stat = statistic, .z = z, .p.val = p)) %>%
#>   tryCatch(
#>     error   = function (e) data.frame(.stat=NA, .z=NA, .p.val=NA), 
#>     warning = function (w) data.frame(.stat=NA, .z=NA, .p.val=NA) ) %>%
#>   data.frame(row.names = NULL, .n = attr(dm, 'Size'), .)
```
