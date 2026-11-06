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
  strata = NULL,
  covariates = NULL,
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

- strata:

  Dataset field with a blocking variable, such as a subject identifier
  for repeated measures. Must be categorical. When given, regression
  models (`test = "emmeans"` or `"emtrends"`) include a random intercept
  per stratum, fitted with
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) as
  `s(strata, bs = "re")` (a linear mixed model for `fit = "lm"` and
  `"log"`). Not supported by the rank tests. Default: `NULL`

- covariates:

  Dataset field(s) to adjust for, such as age or antibiotic use. Each is
  added to the regression model as a linear term (numeric fields) or a
  factor (categorical fields), and the reported means, slopes, and
  differences are adjusted for them:
  [`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html)
  holds numeric covariates at their mean and averages over the levels of
  categorical ones with equal weights. Only for `test = "emmeans"` or
  `"emtrends"`; the rank tests can't adjust for covariates. Rows with an
  `NA` in a covariate are dropped. A covariate that has the same value
  in every observation of a comparison (e.g. within a pair of groups)
  gives an `NA` row and a warning naming it. Default: `NULL`

- transform:

  Transformation to apply to the diversity values. Options are:
  `c("none", "rank", "log", "log1p", "sqrt")`. `"rank"` is useful for
  correcting for non-normally distributions before applying regression
  statistics. Default: `"none"`

- test:

  Method for computing p-values: `'wilcox'`, `'kruskal'`, `'emmeans'`,
  or `'emtrends'`. Default: `'emmeans'`

- fit:

  How to fit the trendline. `'lm'`, `'log'`, or `'gam'`. Default:
  `'gam'`

- at:

  Position(s) along `regr` where group means or local slopes are
  evaluated. Only used by `test = "emmeans"` with `stat.by`, and by
  `test = "emtrends"` with `fit = "gam"` (the other combinations test a
  single, global slope). Default: `NULL`, which evaluates at one
  prespecified position: the median of `regr` (computed within each
  `split.by` stratum). A vector evaluates at each position; the returned
  table then has one row per position, and the p-values are adjusted
  across all of them.

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
  `p.adjust.methods` for a list of available options. The adjustment
  spans the entire returned table: every combination of `split.by`
  strata, `stat.by` pairs, and `at` positions. Default: `"fdr"`

## Value

A tibble data.frame with the fields below. The R code used to generate
the statistics is in `$code`.

|  |  |
|----|----|
| **Field** | **Description** |
| .mean | Estimated marginal mean. See [`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html). |
| .mean.diff | Difference in means (first group minus second). |
| .slope | Trendline slope. See [`emmeans::emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.html). |
| .slope.diff | Difference in slopes (first group minus second). |
| .loc.shift | Hodges-Lehmann estimate of the location shift between two groups (wilcox). |
| .pseudomedian | Hodges-Lehmann estimate of the location (one-sample wilcox). |
| .h1 | Alternative hypothesis. |
| .p.val | The probability, under the null hypothesis, of a test statistic at least as extreme as the one observed. |
| .adj.p | `.p.val` after adjusting for multiple comparisons across the whole table. |
| .effect.size | Cohen's d effect size. See [`emmeans::eff_size()`](https://rvlenth.github.io/emmeans/reference/eff_size.html). |
| .lower | Confidence interval lower bound. |
| .upper | Confidence interval upper bound. |
| .se | Standard error. |
| .n | Number of observations used. |
| .df | Degrees of freedom. |
| .stat | Wilcoxon or Kruskal-Wallis rank sum statistic. |
| .t.ratio | t statistic: estimate divided by `.se`. |
| .edf, .ref.df | Effective and reference degrees of freedom of the smooth term (gam). |
| .f | F statistic for the smooth term (gam). |
| .r.sqr | R-squared (lm); adjusted R-squared (gam). |
| .adj.r | Adjusted R-squared (lm). |
| .dev.expl | Proportion of deviance explained (gam). |
| .aic | Akaike Information Criterion. |
| .bic | Bayesian Information Criterion. |
| .loglik | Log-likelihood of the fitted model. |

## Tests

Which quantity is tested depends on `test` and on whether `regr` and
`stat.by` are given. Rows with an `NA` in any of the `regr`, `resp`,
`stat.by`, `split.by`, `strata`, or `covariates` fields are dropped
(with a message).

|  |  |  |  |
|----|----|----|----|
| `regr` | `stat.by` | `test = "emmeans"` | `test = "emtrends"` |
| no | no | mean of `resp` vs `mu` | (not allowed) |
| no | yes | pairwise differences in group means (`.mean.diff`) | (not allowed) |
| yes | no | `fit = "lm"`/`"log"`: the slope's t-test (`.slope`). `fit = "gam"`: the smooth term's approximate significance (`.edf`, `.f`) | the slope at `at` (`.slope`; global for `"lm"`/`"log"`, local for `"gam"`) |
| yes | yes | pairwise differences in group means at `at` (`.mean.diff`) | pairwise differences in slopes (`.slope.diff`; global for `"lm"`/`"log"`, local at `at` for `"gam"`) |

For `fit = "log"`, the slope is with respect to `log(regr)`. `mu` and
`alt` apply to the mean or slope being tested; they are ignored for
pairwise comparisons (always `!= 0`) and for the smooth-term test.
Pairwise comparisons fit a separate model to each pair of groups.
Without `regr`, `fit` is irrelevant: the model is `resp ~ stat.by`.
`covariates` are added to every model (in the formula, as `.cov1`,
`.cov2`, ...; the table header names them), so that each of the
quantities above is adjusted for them. For example,
`stat.by = "Treatment", covariates = "Age"` compares the groups' means
at the average age.

`test = "wilcox"` runs
[`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)
(one-sample against `mu`, or two-sample for each pair of `stat.by`
groups) with `exact = NULL`, so R chooses an exact or
normal-approximation p-value based on the sample sizes and ties.
`test = "kruskal"` runs
[`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html)
across all `stat.by` groups.

## See also

Other alpha_diversity:
[`adiv_boxplot()`](https://cmmr.github.io/rbiom/reference/adiv_boxplot.md),
[`adiv_corrplot()`](https://cmmr.github.io/rbiom/reference/adiv_corrplot.md),
[`adiv_table()`](https://cmmr.github.io/rbiom/reference/adiv_table.md)

Other stats_tables:
[`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md),
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md),
[`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md),
[`taxa_da()`](https://cmmr.github.io/rbiom/reference/taxa_da.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md)

## Examples

``` r
    library(rbiom)

    biom <- rarefy(hmp50)

    adiv_stats(biom, stat.by = "Sex")[,1:6]
#> # Model:    lm(.resp ~ .stat.by)
#> # A tibble: 1 × 6
#>   Sex           .mean.diff .h1    .p.val  .adj.p .lower
#>   <chr>              <dbl> <fct>   <dbl>   <dbl>  <dbl>
#> 1 Female - Male     -0.770 != 0  0.00932 0.00932  -1.34

    adiv_stats(biom, stat.by = "Sex", split.by = "Body Site")[,1:6]
#> Warning: ! Statistics could not be computed for 1 of 5 comparisons: "Mid vagina / NA".
#> ℹ Message from the underlying tests:
#> • fewer than two groups in this stratum
#> # Model:    lm(.resp ~ .stat.by)
#> # A tibble: 5 × 6
#>   `Body Site`    Sex           .mean.diff .h1   .p.val .adj.p
#>   <fct>          <chr>              <dbl> <fct>  <dbl>  <dbl>
#> 1 Saliva         Female - Male    -0.257  != 0   0.193  0.461
#> 2 Buccal mucosa  Female - Male    -0.544  != 0   0.248  0.461
#> 3 Stool          Female - Male    -0.197  != 0   0.346  0.461
#> 4 Anterior nares Female - Male    -0.0671 != 0   0.736  0.736
#> 5 Mid vagina     NA               NA      NA    NA     NA    

    adiv_stats(biom, stat.by = "Body Site", test = "kruskal")
#> # Model:    kruskal.test(.resp ~ .stat.by)
#> # A tibble: 1 × 5
#>   .stat   .df       .p.val       .adj.p    .n
#>   <dbl> <int>        <dbl>        <dbl> <int>
#> 1  38.9     4 0.0000000736 0.0000000736    49

    adiv_stats(biom, regr = "Age", fit = "lm")
#> # Model:    lm(.resp ~ .regr)
#> # A tibble: 1 × 15
#>   .slope .h1   .p.val .adj.p  .lower .upper    .se    .n   .df .t.ratio  .r.sqr
#>    <dbl> <fct>  <dbl>  <dbl>   <dbl>  <dbl>  <dbl> <int> <dbl>    <dbl>   <dbl>
#> 1 0.0159 != 0   0.625  0.625 -0.0490 0.0807 0.0322    49    47    0.492 0.00512
#> # ℹ 4 more variables: .adj.r <dbl>, .aic <dbl>, .bic <dbl>, .loglik <dbl>
```
