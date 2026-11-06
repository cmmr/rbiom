# Test taxa abundances for associations with metadata.

A convenience wrapper for
[`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md) +
[`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md).

## Usage

``` r
taxa_stats(
  biom,
  regr = NULL,
  stat.by = NULL,
  rank = -1,
  taxa = 6,
  lineage = FALSE,
  unc = "singly",
  other = FALSE,
  split.by = NULL,
  strata = NULL,
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

- lineage:

  Include all ranks in the name of the taxa. For instance, setting to
  `TRUE` will produce
  `Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales`.
  Otherwise the taxa name will simply be `Coriobacteriales`. You want to
  set this to TRUE when `unc = "asis"` and you have taxa names (such as
  *Incertae_Sedis*) that map to multiple higher level ranks. Default:
  `FALSE`

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

- other:

  Sum all non-itemized taxa into an "Other" taxa. When `FALSE`, only
  returns taxa matched by the `taxa` argument. Specifying `TRUE` adds
  "Other" to the returned set. A string can also be given to imply
  `TRUE`, but with that value as the name to use instead of "Other".
  Default: `FALSE`

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

- transform:

  Transformation to apply to calculated values. Options are:
  `c("none", "rank", "log", "log1p", "sqrt")`. `"rank"` is useful for
  correcting for non-normal distributions before applying regression
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
`stat.by`, `split.by`, or `strata` fields are dropped (with a message).

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

`test = "wilcox"` runs
[`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)
(one-sample against `mu`, or two-sample for each pair of `stat.by`
groups) with `exact = NULL`, so R chooses an exact or
normal-approximation p-value based on the sample sizes and ties.
`test = "kruskal"` runs
[`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html)
across all `stat.by` groups.

## See also

Other taxa_abundance:
[`sample_sums()`](https://cmmr.github.io/rbiom/reference/sample_sums.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_clusters()`](https://cmmr.github.io/rbiom/reference/taxa_clusters.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md),
[`taxa_sums()`](https://cmmr.github.io/rbiom/reference/taxa_sums.md),
[`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)

Other stats_tables:
[`adiv_stats()`](https://cmmr.github.io/rbiom/reference/adiv_stats.md),
[`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md),
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md),
[`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md)

## Examples

``` r
    library(rbiom)

    biom <- rarefy(hmp50)

    taxa_stats(biom, stat.by = "Body Site", rank = "Family")[,1:6]
#> Warning: ! Statistics could not be computed for 1 of 60 comparisons: "Lactobacillaceae /
#>   Saliva - Stool".
#> # Model:    lm(.resp ~ .stat.by)
#> # A tibble: 60 × 6
#>    .taxa            `Body Site`                .mean.diff .h1     .p.val  .adj.p
#>    <fct>            <chr>                           <dbl> <fct>    <dbl>   <dbl>
#>  1 Bacteroidaceae   Mid vagina - Stool              -611. != 0  4.59e-10 7.01e-9
#>  2 Bacteroidaceae   Anterior nares - Stool          -610. != 0  4.74e-10 7.01e-9
#>  3 Bacteroidaceae   Buccal mucosa - Stool           -610. != 0  4.75e-10 7.01e-9
#>  4 Bacteroidaceae   Saliva - Stool                  -610. != 0  4.75e-10 7.01e-9
#>  5 Streptococcaceae Buccal mucosa - Mid vagina       788. != 0  1.29e- 8 1.53e-7
#>  6 Streptococcaceae Anterior nares - Buccal m…      -767. != 0  2.15e- 8 2.07e-7
#>  7 Lactobacillaceae Mid vagina - Saliva             1049. != 0  3.14e- 8 2.07e-7
#>  8 Lactobacillaceae Buccal mucosa - Mid vagina     -1049. != 0  3.14e- 8 2.07e-7
#>  9 Lactobacillaceae Anterior nares - Mid vagi…     -1048. != 0  3.15e- 8 2.07e-7
#> 10 Streptococcaceae Buccal mucosa - Stool            789. != 0  4.88e- 8 2.88e-7
#> # ℹ 50 more rows
```
