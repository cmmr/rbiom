# Run non-parametric statistics on a data.frame.

A simple interface to lower-level statistics functions, including
[`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html),
[`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html),
[`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html),
and
[`emmeans::emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.html).

## Usage

``` r
stats_table(
  df,
  regr = NULL,
  resp = attr(df, "response"),
  stat.by = NULL,
  split.by = NULL,
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

- split.by:

  Dataset field(s) that the data should be split by prior to any
  calculations. Must be categorical. Default: `NULL`

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

Other stats_tables:
[`adiv_stats()`](https://cmmr.github.io/rbiom/reference/adiv_stats.md),
[`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md),
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md)

## Examples

``` r
    library(rbiom)
    
    biom <- rarefy(hmp50)
    
    df <- taxa_table(biom, rank = "Family")
    stats_table(df, stat.by = "Body Site")[,1:6]
#> # Model:    gam(.abundance ~ `Body Site`, method = "REML")
#> # A tibble: 10 × 6
#>    `Body Site`                    .mean.diff .h1   .p.val .adj.p .effect.size
#>    <chr>                               <dbl> <fct>  <dbl>  <dbl>        <dbl>
#>  1 Mid vagina - Saliva                 83.4  != 0   0.140  0.506       0.271 
#>  2 Anterior nares - Mid vagina        -85.5  != 0   0.161  0.506      -0.257 
#>  3 Buccal mucosa - Saliva              60.3  != 0   0.163  0.506       0.257 
#>  4 Anterior nares - Buccal mucosa     -62.5  != 0   0.203  0.506      -0.234 
#>  5 Mid vagina - Stool                  73.9  != 0   0.255  0.510       0.215 
#>  6 Buccal mucosa - Stool               50.9  != 0   0.331  0.552       0.183 
#>  7 Buccal mucosa - Mid vagina         -23.1  != 0   0.732  0.873      -0.0627
#>  8 Anterior nares - Stool             -11.6  != 0   0.785  0.873      -0.0512
#>  9 Saliva - Stool                      -9.47 != 0   0.786  0.873      -0.0511
#> 10 Anterior nares - Saliva             -2.12 != 0   0.947  0.947      -0.0121
    
    df <- adiv_table(biom)
    stats_table(df, stat.by = "Sex", split.by = "Body Site")[,1:7]
#> # Model:    gam(.diversity ~ Sex, method = "REML")
#> # A tibble: 5 × 7
#>   `Body Site`    Sex           .mean.diff .h1   .p.val .adj.p .effect.size
#>   <fct>          <chr>              <dbl> <fct>  <dbl>  <dbl>        <dbl>
#> 1 Saliva         Female - Male    -0.257  != 0   0.193  0.461       -0.899
#> 2 Buccal mucosa  Female - Male    -0.544  != 0   0.248  0.461       -0.788
#> 3 Stool          Female - Male    -0.197  != 0   0.346  0.461       -0.678
#> 4 Anterior nares Female - Male    -0.0671 != 0   0.736  0.736       -0.221
#> 5 Mid vagina     NA               NA      NA    NA     NA           NA    
```
