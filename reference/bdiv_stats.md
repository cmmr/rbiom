# Test beta diversity for associations with metadata.

A convenience wrapper for
[`bdiv_table()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md) +
[`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md).

## Usage

``` r
bdiv_stats(
  biom,
  regr = NULL,
  stat.by = NULL,
  bdiv = "bray",
  weighted = NULL,
  tree = NULL,
  within = NULL,
  between = NULL,
  split.by = NULL,
  transform = "none",
  test = "emmeans",
  fit = "gam",
  at = NULL,
  level = 0.95,
  alt = "!=",
  mu = 0,
  p.adj = "fdr",
  alpha = 0.5,
  cpus = n_cpus()
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

- bdiv:

  Beta diversity distance algorithm(s) to use. Options are:
  `c("aitchison", "bhattacharyya", "bray", "canberra", "chebyshev", "chord", "clark", "sorensen", "divergence", "euclidean", "generalized_unifrac", "gower", "hamming", "hellinger", "horn", "jaccard", "jensen", "jsd", "lorentzian", "manhattan", "matusita", "minkowski", "morisita", "motyka", "normalized_unifrac", "ochiai", "psym_chisq", "soergel", "squared_chisq", "squared_chord", "squared_euclidean", "topsoe", "unweighted_unifrac", "variance_adjusted_unifrac", "wave_hedges", "weighted_unifrac")`.
  For the UniFrac family, a phylogenetic tree must be present in `biom`
  or explicitly provided via `tree=`. Supports partial matching.
  Multiple values are allowed for functions which return a table or
  plot. Default: `"bray"`

- weighted:

  (Deprecated - weighting is now inherent in bdiv metric name.) Take
  relative abundances into account. When `weighted=FALSE`, only
  presence/absence is considered. Multiple values allowed. Default:
  `NULL`

- tree:

  A `phylo` object representing the phylogenetic relationships of the
  taxa in `biom`. Only required when computing UniFrac distances.
  Default: `biom$tree`

- within, between:

  Dataset field(s) for intra- or inter- sample comparisons.
  Alternatively, dataset field names given elsewhere can be prefixed
  with `'=='` or `'!='` to assign them to `within` or `between`,
  respectively. Default: `NULL`

- split.by:

  Dataset field(s) that the data should be split by prior to any
  calculations. Must be categorical. Default: `NULL`

- transform:

  Transformation to apply to calculated values. Options are:
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

- alpha:

  The alpha term to use in Generalized UniFrac. How much weight to give
  to relative abundances; a value between 0 and 1, inclusive. Setting
  `alpha=1` is equivalent to Normalized UniFrac. Default: `0.5`

- cpus:

  The number of CPUs to use. Set to `NULL` to use all available, or to
  `1` to disable parallel processing. Default: `NULL`

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

Other beta_diversity:
[`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md),
[`bdiv_clusters()`](https://cmmr.github.io/rbiom/reference/bdiv_clusters.md),
[`bdiv_corrplot()`](https://cmmr.github.io/rbiom/reference/bdiv_corrplot.md),
[`bdiv_heatmap()`](https://cmmr.github.io/rbiom/reference/bdiv_heatmap.md),
[`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md),
[`bdiv_ord_table()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_table.md),
[`bdiv_table()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md),
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md)

Other stats_tables:
[`adiv_stats()`](https://cmmr.github.io/rbiom/reference/adiv_stats.md),
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md),
[`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md)

## Examples

``` r
    library(rbiom)
    
    biom <- rarefy(hmp50)
      
    bdiv_stats(biom, stat.by = "Sex", bdiv = c("bray", "w_unifrac"))[,1:7]
#> # Model:    gam(.distance ~ Sex, method = "REML")
#> # A tibble: 6 × 7
#>   .bdiv     Sex                     .mean.diff .h1    .p.val .adj.p .effect.size
#>   <fct>     <chr>                        <dbl> <fct>   <dbl>  <dbl>        <dbl>
#> 1 bray      Female - Male              0.0569  != 0  0.00642 0.0308       0.247 
#> 2 bray      Male - Female vs Male     -0.0468  != 0  0.0103  0.0308      -0.224 
#> 3 w_unifrac Female - Female vs Male   -0.0262  != 0  0.0611  0.122       -0.119 
#> 4 w_unifrac Male - Female vs Male     -0.0188  != 0  0.290   0.436       -0.0923
#> 5 bray      Female - Female vs Male    0.0100  != 0  0.472   0.566        0.0458
#> 6 w_unifrac Female - Male             -0.00741 != 0  0.720   0.720       -0.0323
    
    biom <- subset(biom, `Body Site` %in% c('Saliva', 'Stool', 'Buccal mucosa'))
    bdiv_stats(biom, stat.by = "Body Site", split.by = "==Sex")[,1:6]
#> # Model:    gam(.distance ~ `Body Site`, method = "REML")
#> # A tibble: 30 × 6
#>    Sex    `Body Site`                         .mean.diff .h1     .p.val   .adj.p
#>    <fct>  <chr>                                    <dbl> <fct>    <dbl>    <dbl>
#>  1 Female Buccal mucosa - Buccal mucosa vs S…     -0.792 != 0  8.67e-36 1.50e-34
#>  2 Female Buccal mucosa - Saliva vs Stool         -0.790 != 0  1.00e-35 1.50e-34
#>  3 Male   Saliva - Saliva vs Stool                -0.511 != 0  3.25e-23 2.47e-22
#>  4 Male   Saliva - Buccal mucosa vs Stool         -0.511 != 0  3.29e-23 2.47e-22
#>  5 Female Buccal mucosa - Buccal mucosa vs S…     -0.590 != 0  2.01e-22 1.21e-21
#>  6 Female Saliva - Buccal mucosa vs Stool         -0.517 != 0  3.04e-21 1.50e-20
#>  7 Female Saliva - Saliva vs Stool                -0.515 != 0  3.50e-21 1.50e-20
#>  8 Female Buccal mucosa vs Saliva - Buccal m…     -0.202 != 0  4.73e-20 1.77e-19
#>  9 Female Buccal mucosa vs Saliva - Saliva v…     -0.200 != 0  7.30e-20 2.43e-19
#> 10 Male   Stool - Saliva vs Stool                 -0.403 != 0  2.59e-18 7.20e-18
#> # ℹ 20 more rows
```
