# Test alpha diversity for associations with metadata.

A convenience wrapper for
[`adiv_table()`](https://cmmr.github.io/rbiom/reference/adiv_table.md) +
[`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md).

## Usage

``` r
adiv_stats(
  biom,
  regr = NULL,
  stat.by = NULL,
  adiv = "Shannon",
  split.by = NULL,
  transform = "none",
  test = "emmeans",
  fit = "gam",
  at = NULL,
  level = 0.95,
  alt = "!=",
  mu = 0,
  p.adj = "fdr"
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- regr:

  Dataset field with the x-axis (independent; predictive) values. Must
  be numeric. Default: `NULL`

- stat.by:

  Dataset field with the statistical groups. Must be categorical.
  Default: `NULL`

- adiv:

  Alpha diversity metric(s) to use. Options are:
  `c("ace", "berger", "brillouin", "chao1", "faith", "fisher", "simpson", "inv_simpson", "margalef", "mcintosh", "menhinick", "observed", "shannon", "squares")`.
  For `"faith"`, a phylogenetic tree must be present in `biom` or
  explicitly provided via `tree=`. Set `adiv=".all"` to use all metrics.
  Multiple/abbreviated values allowed. Default: `"shannon"`

- split.by:

  Dataset field(s) that the data should be split by prior to any
  calculations. Must be categorical. Default: `NULL`

- transform:

  Transformation to apply. Options are:
  `c("none", "rank", "log", "log1p", "sqrt", "percent")`. `"rank"` is
  useful for correcting for non-normally distributions before applying
  regression statistics. Default: `"none"`

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

- level:

  The confidence level for calculating a confidence interval. Default:
  `0.95`

- alt:

  Alternative hypothesis direction. Options are `'!='` (two-sided; not
  equal to `mu`), `'<'` (less than `mu`), or `'>'` (greater than `mu`).
  Default: `'!='`

- mu:

  Reference value to test against. Default: `0`

- p.adj:

  Method to use for multiple comparisons adjustment of p-values. Run
  `p.adjust.methods` for a list of available options. Default: `"fdr"`

## Value

A tibble data.frame with fields from the table below. This tibble object
provides the `$code` operator to print the R code used to generate the
statistics.

|              |                                                                                                                |
|--------------|----------------------------------------------------------------------------------------------------------------|
| **Field**    | **Description**                                                                                                |
| .mean        | Estimated marginal mean. See [`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html). |
| .mean.diff   | Difference in means.                                                                                           |
| .slope       | Trendline slope. See [`emmeans::emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.html).       |
| .slope.diff  | Difference in slopes.                                                                                          |
| .h1          | Alternate hypothesis.                                                                                          |
| .p.val       | Probability that null hypothesis is correct.                                                                   |
| .adj.p       | `.p.val` after adjusting for multiple comparisons.                                                             |
| .effect.size | Effect size. See [`emmeans::eff_size()`](https://rvlenth.github.io/emmeans/reference/eff_size.html).           |
| .lower       | Confidence interval lower bound.                                                                               |
| .upper       | Confidence interval upper bound.                                                                               |
| .se          | Standard error.                                                                                                |
| .n           | Number of samples.                                                                                             |
| .df          | Degrees of freedom.                                                                                            |
| .stat        | Wilcoxon or Kruskal-Wallis rank sum statistic.                                                                 |
| .t.ratio     | `.mean` / `.se`                                                                                                |
| .r.sqr       | Percent of variation explained by the model.                                                                   |
| .adj.r       | `.r.sqr`, taking degrees of freedom into account.                                                              |
| .aic         | Akaike Information Criterion (predictive models).                                                              |
| .bic         | Bayesian Information Criterion (descriptive models).                                                           |
| .loglik      | Log-likelihood goodness-of-fit score.                                                                          |
| .fit.p       | P-value for observing this fit by chance.                                                                      |

## See also

Other alpha_diversity:
[`adiv_boxplot()`](https://cmmr.github.io/rbiom/reference/adiv_boxplot.md),
[`adiv_corrplot()`](https://cmmr.github.io/rbiom/reference/adiv_corrplot.md),
[`adiv_table()`](https://cmmr.github.io/rbiom/reference/adiv_table.md)

Other stats_tables:
[`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md),
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md),
[`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md)

## Examples

``` r
    library(rbiom) 
    
    biom <- rarefy(hmp50)
      
    adiv_stats(biom, stat.by = "Sex")[,1:6]
#> # Model:    gam(.diversity ~ Sex, method = "REML")
#> # A tibble: 1 × 6
#>   Sex           .mean.diff .h1    .p.val  .adj.p .effect.size
#>   <chr>              <dbl> <fct>   <dbl>   <dbl>        <dbl>
#> 1 Female - Male     -0.770 != 0  0.00932 0.00932       -0.795
      
    adiv_stats(biom, stat.by = "Sex", split.by = "Body Site")[,1:6]
#> # Model:    gam(.diversity ~ Sex, method = "REML")
#> # A tibble: 5 × 6
#>   `Body Site`    Sex           .mean.diff .h1   .p.val .adj.p
#>   <fct>          <chr>              <dbl> <fct>  <dbl>  <dbl>
#> 1 Saliva         Female - Male    -0.257  != 0   0.193  0.461
#> 2 Buccal mucosa  Female - Male    -0.544  != 0   0.248  0.461
#> 3 Stool          Female - Male    -0.197  != 0   0.346  0.461
#> 4 Anterior nares Female - Male    -0.0671 != 0   0.736  0.736
#> 5 Mid vagina     NA               NA      NA    NA     NA    
    
    adiv_stats(biom, stat.by = "Body Site", test = "kruskal")
#> # Model:    kruskal.test(.diversity ~ `Body Site`)
#> # A tibble: 1 × 6
#>   .stat .h1         .p.val       .adj.p    .n   .df
#>   <dbl> <fct>        <dbl>        <dbl> <int> <int>
#> 1  38.9 > 0   0.0000000736 0.0000000736    49     4
```
