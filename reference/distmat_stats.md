# Run permutation statistics on a distance matrix.

Tests whether the distances in `dm` are structured by `groups`, using
permutation tests. Samples, not distances, are the unit of permutation,
so the n(n-1)/2 pairwise distances are never treated as independent
observations.

## Usage

``` r
distmat_stats(
  dm,
  groups,
  test = "adonis2",
  seed = 0,
  permutations = 999,
  strata = NULL
)
```

## Arguments

- dm:

  A `dist`-class distance matrix, as returned from
  [`bdiv_distmat()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md)
  or [`stats::dist()`](https://rdrr.io/r/stats/dist.html). Required.

- groups:

  A vector of grouping values, one per sample. Either named with the
  sample names in `attr(dm, 'Labels')`, or unnamed and in the same
  order. Categorical values define groups; numeric values are used as a
  covariate (`test = "adonis2"`) or correlated with the distances
  (`test = "mantel"`). Samples with an `NA` value are dropped (with a
  message). Required.

- test:

  The permutation test to run. Options are:

  `"adonis2"` -

  :   PERMANOVA;
      [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html).
      Accepts categorical or numeric `groups`.

  `"mrpp"` -

  :   Multiple response permutation procedure;
      [`vegan::mrpp()`](https://vegandevs.github.io/vegan/reference/mrpp.html).
      Categorical `groups` only.

  `"betadisper"` -

  :   Homogeneity of multivariate dispersion, i.e. whether groups differ
      in their spread rather than their location;
      [`vegan::betadisper()`](https://vegandevs.github.io/vegan/reference/betadisper.html)
      followed by
      [`vegan::permutest()`](https://vegandevs.github.io/vegan/reference/anova.cca.html).
      Categorical `groups` only. A significant PERMANOVA together with a
      significant dispersion test means that group differences may be
      due to differences in spread.

  `"mantel"` -

  :   Mantel test of the correlation (Spearman) between `dm` and the
      pairwise absolute differences in a numeric `groups`;
      [`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html).
      Numeric `groups` only.

  `"none"` -

  :   Don't run any statistics; returns `NULL`.

  Abbreviations are allowed. Default: `"adonis2"`

- seed:

  Random seed for permutations, rarefaction, random tie-breaking, and
  stochastic ordinations. Results are reproducible for a given seed, and
  the caller's random number generator state is left unchanged. Must be
  a non-negative integer. Default: `0`

- permutations:

  Number of random permutations to use. Default: `999`

- strata:

  A blocking variable for repeated measures, e.g. a subject identifier.
  Given like `groups` (one value per sample, named or in
  `attr(dm, 'Labels')` order). How permutations are restricted depends
  on whether `groups` varies within strata; see Details. Default: `NULL`
  (samples are exchangeable).

## Value

A tibble data.frame with one row and the columns:

- `.n` - :

  The number of samples in the distance matrix.

- `.stat` - :

  The observed test statistic: the pseudo-F ratio (adonis2, betadisper),
  the weighted mean within-group distance delta (mrpp), or Spearman's r
  (mantel).

- `.r.sqr` - :

  adonis2 only: the proportion of the total sum of squares explained by
  `groups`.

- `.z` - :

  mrpp only: the observed statistic minus the mean of the permuted
  statistics, divided by their standard deviation.

- `.p.val` - :

  The permutation p-value: the probability, under the null hypothesis,
  of a statistic at least as extreme as the observed one, estimated from
  the permutations (plus the observed data). Its resolution is limited
  by `permutations` (or by the number of distinct arrangements of whole
  strata; see Details).

The R code used to calculate these statistics is in `$code`. The
`groups` (and `strata`) values are embedded in that code so that it can
be re-run as is.

## Details

**Permutation schemes with `strata`.** Samples from the same stratum
(e.g. repeated samples of one subject) are not independent, so they
can't be shuffled freely. Which restriction is valid depends on
`groups`:

- **Within-strata variable** (`groups` varies within at least one
  stratum, e.g. time point or a within-subject treatment): labels are
  permuted only within each stratum, using
  `permute::how(nperm = permutations, blocks = strata)` with the vegan
  function for `test`.

- **Between-strata variable** (`groups` is constant within every
  stratum, e.g. delivery mode with subject as the stratum): permuting
  within strata would never change the labels (every permutation would
  reproduce the observed statistic, giving p = 1). Instead, the
  stratum-level values are permuted among whole strata (strata may have
  different sizes), and the statistic is recomputed on all samples. The
  p-value is `(1 + #{permuted >= observed}) / (1 + permutations)`. When
  there are no more distinct arrangements of the values among strata
  than `permutations`, all of them are enumerated instead and the
  p-value is exact: the proportion of arrangements (including the
  observed one) with a statistic at least as large. The statistics are
  computed by rbiom with the same formulas as vegan: the PERMANOVA
  pseudo-F and R-squared from the Gower-centred matrix of squared
  distances, and Spearman's r between the distances and pairwise
  differences in `groups`. The effective sample size is the number of
  strata, so a warning is given when the null distribution is coarse
  (e.g. a group present in only one stratum). Only `test = "adonis2"`
  and `test = "mantel"` support this scheme; `"mrpp"` and `"betadisper"`
  give an error.

The scheme used, and the number of permutations or enumerated
arrangements, is reported in the table header.

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
[`taxa_da()`](https://cmmr.github.io/rbiom/reference/taxa_da.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md)

## Examples

``` r
    library(rbiom)

    hmp10 <- hmp50[1:10]

    dm <- bdiv_distmat(hmp10, 'w_unifrac')

    distmat_stats(dm, groups = pull(hmp10, 'Body Site'))
#> # Test:     PERMANOVA (vegan::adonis2) ~ groups; 999 permutations.
#> # A tibble: 1 × 4
#>      .n .stat .r.sqr .p.val
#>   <int> <dbl>  <dbl>  <dbl>
#> 1    10  9.27  0.726  0.001

    distmat_stats(dm, groups = pull(hmp10, 'Age'), test = "mantel")
#> # Test:     Mantel (vegan::mantel) ~ groups; 999 permutations.
#> # A tibble: 1 × 3
#>      .n .stat .p.val
#>   <int> <dbl>  <dbl>
#> 1    10 0.279  0.081

    # See the R code used to calculate these statistics:
    stats <- distmat_stats(dm, groups = pull(hmp10, 'Age'), test = "mantel")
    stats$code
#> # Generated by rbiom 4.0.0.9001
#> library(rbiom)
#> 
#> biom <- hmp50
#> biom <- biom[1:10]
#> 
#> data <- bdiv_distmat(biom, "w_unifrac")
#> 
#> dm     <- data
#> groups <- c(HMP01 = 22, HMP02 = 24, HMP03 = 28, HMP04 = 25, HMP05 = 27, HMP06 = 32, HMP07 = 26, HMP08 = 27, HMP09 = 33, 
#>   HMP10 = 22)
#> 
#> grouping <- groups[attr(dm, 'Labels')]
#> set.seed(0)
#> 
#> ptest <- vegan::mantel(dm, stats::dist(grouping), method = 'spearman', permutations = 999)
#> stats <- tibble::tibble(
#>   .n     = attr(dm, 'Size'),
#>   .stat  = ptest$statistic,
#>   .p.val = ptest$signif )

    # Repeated measures: delivery mode is constant within each subject,
    # so whole subjects are permuted (all 220 arrangements of the 12
    # infants are enumerated, giving an exact p-value).
    b  <- rarefy(subset(babies, `Age (days)` <= 10))
    dm <- bdiv_distmat(b, 'bray')
    distmat_stats(dm, pull(b, 'Delivery mode'), strata = pull(b, 'Subject ID'))
#> # Test:     PERMANOVA (adonis2 pseudo-F) ~ groups; whole strata permuted
#> #   (between-strata variable; all 220 arrangements, exact).
#> # A tibble: 1 × 4
#>      .n .stat .r.sqr .p.val
#>   <int> <dbl>  <dbl>  <dbl>
#> 1   120  4.81 0.0392  0.255
```
