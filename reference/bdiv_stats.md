# Test beta diversity for associations with metadata.

Pairwise distances are not independent observations, so beta diversity
is tested with permutation tests in which samples, not distances, are
permuted: PERMANOVA and friends for a categorical `stat.by`, or a Mantel
test for a numeric `regr`. See
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md)
for the tests.

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
  strata = NULL,
  test = if (is.null(regr)) "adonis2" else "mantel",
  seed = 0,
  permutations = 999,
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

  Dataset field with a numeric variable to correlate the distances with
  (Mantel test). When both `regr` and `stat.by` are given, the Mantel
  test is run within each `stat.by` group. Default: `NULL`

- stat.by:

  Dataset field with the groups to compare. Must be categorical. Any
  `==`/`!=` prefix is ignored. Default: `NULL`

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

  Accepted for compatibility with
  [`bdiv_table()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md)
  and
  [`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md),
  where they select which sample *pairs* are displayed. The tests here
  are computed on samples, so these arguments have no effect.

- split.by:

  Dataset field(s) that the data should be split by prior to any
  calculations. Must be categorical. Default: `NULL`

- strata:

  Dataset field with a blocking variable, such as a subject identifier
  for repeated measures. Must be categorical. When `stat.by` (or `regr`)
  varies within strata, permutations are restricted to within each
  stratum; when it is constant within every stratum (a between-subject
  variable), whole strata are permuted instead. See the Details of
  [`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md).
  `test = "mrpp"` and `"betadisper"` are not available for a
  between-strata variable. Default: `NULL`

- test:

  The permutation test. With `stat.by`: `"adonis2"` (PERMANOVA; the
  default), `"mrpp"`, or `"betadisper"` (homogeneity of dispersion).
  With `regr`: `"mantel"` (the default) or `"adonis2"` (PERMANOVA with a
  numeric covariate). See
  [`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md).

- seed:

  Random seed for permutations, rarefaction, random tie-breaking, and
  stochastic ordinations. Results are reproducible for a given seed, and
  the caller's random number generator state is left unchanged. Must be
  a non-negative integer. Default: `0`

- permutations:

  Number of random permutations to use. Default: `999`

- p.adj:

  Method to use for multiple comparisons adjustment of p-values. Run
  `p.adjust.methods` for a list of available options. The adjustment
  spans the entire returned table. Default: `"fdr"`

- alpha:

  The alpha term to use in Generalized UniFrac. How much weight to give
  to relative abundances; a value between 0 and 1, inclusive. Setting
  `alpha=1` is equivalent to Normalized UniFrac. Default: `0.5`

- cpus:

  The number of CPUs to use. Set to `NULL` to use all available, or to
  `1` to disable parallel processing. Does not affect results. Default:
  [`ecodive::n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.html)
  (all available cores)

## Value

A tibble data.frame with one row per test. With a categorical `stat.by`
and more than two groups, each `split.by` stratum (and `bdiv` metric)
gets an omnibus test across all groups plus one pairwise test per pair
of groups, computed on that pair's samples. The columns are:

- `split.by` fields, `.bdiv` - :

  The stratum. `.bdiv` is only present when several `bdiv` metrics were
  requested.

- `stat.by` field - :

  The groups compared, joined by `" vs "`. The omnibus row lists all
  groups.

- `.groups` - :

  The number of groups in the comparison: 2 for pairwise rows, more for
  the omnibus row.

- `.n`, `.stat`, `.r.sqr`, `.z`, `.p.val` - :

  See
  [`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md).

- `.adj.p` - :

  `.p.val` adjusted for multiple comparisons across the whole table.

With a numeric `regr`, one Mantel test per stratum, with `.n`, `.stat`
(Spearman's r), `.p.val`, and `.adj.p`. The R code used to compute the
table is in `$code`.

## Details

With `strata` (e.g. repeated samples per subject), the permutation
scheme depends on the variable being tested. A variable that varies
within strata (e.g. time point) is permuted within each stratum. A
variable that is constant within every stratum (e.g. delivery mode, with
subject as the stratum) is permuted at the stratum level: whole strata
exchange their values, so the effective sample size is the number of
strata. When the number of distinct arrangements is no larger than
`permutations`, all of them are enumerated and the p-value is exact. The
scheme is chosen separately for each row of the table (a pair of groups
can be a between-strata comparison even when the omnibus test is not)
and is reported in the table header. See
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md)
for details.

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

    bdiv_stats(biom, stat.by = "Sex", bdiv = c("bray", "w_unifrac"))
#> # Test:     PERMANOVA ~ Sex; 999 permutations.
#> # A tibble: 2 × 8
#>   .bdiv            Sex            .groups    .n .stat .r.sqr .p.val .adj.p
#>   <chr>            <chr>            <int> <int> <dbl>  <dbl>  <dbl>  <dbl>
#> 1 weighted_unifrac Female vs Male       2    49  2.19 0.0444  0.075  0.088
#> 2 bray             Female vs Male       2    49  1.62 0.0332  0.088  0.088

    biom <- subset(biom, `Body Site` %in% c('Saliva', 'Stool', 'Buccal mucosa'))
    bdiv_stats(biom, stat.by = "Body Site", split.by = "Sex")
#> # Test:     PERMANOVA ~ `Body Site`; 999 permutations.
#> # A tibble: 8 × 8
#>   Sex    `Body Site`                    .groups    .n .stat .r.sqr .p.val .adj.p
#>   <fct>  <chr>                            <int> <int> <dbl>  <dbl>  <dbl>  <dbl>
#> 1 Female Buccal mucosa vs Saliva vs St…       3    15 14.0   0.701  0.001  0.004
#> 2 Male   Buccal mucosa vs Saliva vs St…       3    14 10.1   0.647  0.001  0.004
#> 3 Female Saliva vs Stool                      2    10 10.6   0.569  0.007  0.016
#> 4 Male   Buccal mucosa vs Stool               2     9 11.8   0.627  0.011  0.016
#> 5 Female Buccal mucosa vs Stool               2    10 16.7   0.675  0.012  0.016
#> 6 Male   Saliva vs Stool                      2     9 11.7   0.625  0.013  0.016
#> 7 Female Buccal mucosa vs Saliva              2    10 17.9   0.691  0.016  0.016
#> 8 Male   Buccal mucosa vs Saliva              2    10  6.37  0.443  0.016  0.016

    bdiv_stats(biom, regr = "Age")
#> # Test:     Mantel (Spearman) ~ Age; 999 permutations.
#> # A tibble: 1 × 4
#>      .n   .stat .p.val .adj.p
#>   <int>   <dbl>  <dbl>  <dbl>
#> 1    29 0.00817  0.399  0.399

    # Repeated measures: delivery mode is constant within each subject,
    # so whole subjects are permuted.
    b <- rarefy(subset(babies, `Age (days)` <= 30))
    bdiv_stats(b, stat.by = "Delivery mode", strata = "Subject ID")
#> # Test:     PERMANOVA ~ `Delivery mode`; whole `Subject ID` strata permuted
#> #   (between-strata variable; all 220 arrangements, exact).
#> # A tibble: 1 × 7
#>   `Delivery mode`     .groups    .n .stat .r.sqr .p.val .adj.p
#>   <chr>                 <int> <int> <dbl>  <dbl>  <dbl>  <dbl>
#> 1 Cesarean vs Vaginal       2   360  9.17 0.0250    0.5    0.5
```
