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
#> # Model:    gam(.abundance ~ `Body Site`, method = "REML")
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
