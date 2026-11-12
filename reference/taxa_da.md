# Differential abundance of taxa.

Finds the taxa whose abundance differs between groups (`stat.by`) or
changes with a numeric variable (`regr`). Unlike
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md),
which tests each taxon's abundance on its own, these methods account for
the compositional nature of sequencing data: a sample's counts are a
sample from the community, so an increase in one taxon lowers the
observed proportions of all the others. The method is LinDA (Zhou et al.
2022), built into rbiom.

## Usage

``` r
taxa_da(
  biom,
  stat.by = NULL,
  regr = NULL,
  rank = -1,
  taxa = NULL,
  method = "linda",
  covariates = NULL,
  split.by = NULL,
  strata = NULL,
  prev.cut = 0.1,
  level = 0.95,
  p.adj = "fdr",
  lineage = FALSE,
  unc = "singly"
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- stat.by:

  Dataset field with the groups to compare. Must be categorical. Each
  pair of groups is compared, in the order of the field's levels: with
  levels `A`, `B`, and `C`, the comparisons are `"B - A"`, `"C - A"`,
  and `"C - B"`, and a positive log fold change means the taxon is more
  abundant in the first-named group. Give either `stat.by` or `regr`,
  not both. Default: `NULL`

- regr:

  Dataset field with a numeric variable. The log fold change is then per
  unit increase in `regr`. Give either `stat.by` or `regr`, not both.
  Default: `NULL`

- rank:

  What rank(s) of taxa to test. E.g. `"Phylum"`, `"Genus"`, `".otu"`,
  etc. An integer vector can also be given, where `1` is the highest
  rank, `-1` is the lowest rank, and `0` is the OTU "rank". Each rank is
  analyzed separately. Default: `-1`

- taxa:

  Which taxa to report. An integer value reports the top n most abundant
  taxa, a value 0 \<= n \< 1 reports the taxa with that mean relative
  abundance or greater, and a character vector reports the named taxa.
  Every taxon that passes `prev.cut` is still part of the analysis,
  since LinDA estimates the compositional bias from all of them, and
  p-values are adjusted for all of them. Default: `NULL` (all tested
  taxa).

- method:

  The differential abundance method. Currently only `"linda"`. See
  Details. Default: `"linda"`

- covariates:

  Dataset field(s) to adjust for, such as age or antibiotic use. They
  are included in every model, but their own effects are not reported.
  Default: `NULL`

- split.by:

  Dataset field(s) that the samples should be split by. Each subset is
  analyzed separately. Must be categorical. Default: `NULL`

- strata:

  Dataset field with a blocking variable, such as a subject identifier
  for repeated measures. Must be categorical. Adds a random intercept
  per stratum to every model, fitted as linear mixed models with
  [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html)
  (requires the lmerTest package). Default: `NULL`

- prev.cut:

  Prevalence filter: taxa found in fewer than this proportion of the
  samples are not tested. Default: `0.1`

- level:

  The confidence level for the confidence intervals (`.lower` and
  `.upper`). Default: `0.95`

- p.adj:

  Method to use for multiple comparisons adjustment of p-values. Run
  `p.adjust.methods` for a list of available options. As is standard for
  differential abundance, p-values are adjusted across taxa separately
  for each comparison (and each rank and `split.by` subset). Default:
  `"fdr"`

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

## Value

A tibble data.frame with one row per taxon and comparison, and the
columns:

- `.rank` - :

  The rank. Only present when several ranks were requested.

- `split.by` fields - :

  The subset of samples.

- `stat.by` field - :

  The comparison, e.g. `"Stool - Saliva"` (positive fold changes: more
  abundant in Stool). Absent when `regr` is given.

- `.taxa` - :

  The taxon.

- `.lfc` - :

  The bias-corrected log2 fold change: the difference in log2 abundance
  between the two groups, or the change per unit of `regr`.

- `.se` - :

  The standard error of `.lfc`.

- `.lower`, `.upper` - :

  The confidence interval of `.lfc` at `level`.

- `.t.ratio` - :

  The test statistic, `.lfc / .se`.

- `.df` - :

  Degrees of freedom of the t-test.

- `.p.val` - :

  The p-value.

- `.adj.p` - :

  `.p.val` adjusted for multiple comparisons across taxa (see `p.adj`).

- `.n` - :

  The number of samples analyzed.

Rows are sorted by `.p.val`. The R code that reproduces the table is in
`$code`.

## Details

LinDA expects unrarefied counts: it models differences in sequencing
depth itself, and rarefying discards data.

**Counts.** The counts are summed to `rank` with
[`taxa_matrix()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)
(using `unc` and `lineage`) and every sample with a value for all of the
fields used is included. Samples with an `NA` in `stat.by`, `regr`,
`covariates`, `split.by`, or `strata` are dropped (with a message).

**LinDA** (Zhou et al. 2022) fits one linear model per taxon to the
centered log2-ratios of the counts. Compositionality shifts every
taxon's coefficient by the same unknown amount; assuming that most taxa
don't change, that shift is the mode of the coefficients across taxa,
and subtracting it gives the log2 fold changes. rbiom's implementation
follows the published algorithm and the defaults of
`MicrobiomeStat::linda()`:

1.  Taxa present in fewer than `prev.cut` of the samples are removed,
    and then any samples without counts.

2.  Each taxon's relative abundances are winsorized at its 97th
    percentile, and converted back to (rounded) counts.

3.  Zeros: when the sequencing depth (log of the sample sums) is
    associated with the model's variables (any p \<= 0.1), each zero is
    imputed as the sample's depth divided by the largest depth among the
    samples where the taxon is absent. Otherwise a pseudo-count of 0.5
    is added to every count.

4.  Centered log2-ratio transformation (each sample's log2 counts minus
    their mean).

5.  A linear model per taxon, `clr ~ stat.by + covariates` (or `regr`),
    fitted with [`stats::lm()`](https://rdrr.io/r/stats/lm.html), or
    with
    [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html) and
    a random intercept per stratum when `strata` is given. Numeric
    variables are standardized for fitting; `regr`'s fold changes are
    reported per unit.

6.  Each comparison is a coefficient (`regr`, or a group versus the
    first group) or a difference of two coefficients (two other groups).
    Its bias is its mode across taxa, estimated by the mean-shift
    algorithm with a Gaussian kernel (the same estimator as
    `modeest::mlv(method = "meanshift")`), and is subtracted.

7.  The standard error comes from the coefficients' covariance, and the
    p-value from a t-test with the model's residual degrees of freedom
    (Satterthwaite's for mixed models).

For a numeric `regr` and for comparisons with the first group, this
gives the results of `MicrobiomeStat::linda()`. For a comparison between
two other groups, `MicrobiomeStat::linda.wald.test()` would subtract the
difference of the two coefficients' modes; rbiom takes the mode of the
difference instead, which assumes that most taxa don't differ between
those two groups (rather than between each of them and the first group),
and doesn't depend on which group comes first.

One deliberate difference from `MicrobiomeStat::linda()` (version 1.4):
when it reports that sequencing depth is associated with the model, it
still adds the pseudo-count, because it compares `"Imputation"` with
`"imputation"`; rbiom imputes, as the paper describes.

## References

Zhou H, He K, Chen J, Zhang X (2022). LinDA: linear models for
differential abundance analysis of microbiome compositional data.
*Genome Biology* 23:95.
[doi:10.1186/s13059-022-02655-5](https://doi.org/10.1186/s13059-022-02655-5)

## See also

[`taxa_forest()`](https://cmmr.github.io/rbiom/reference/taxa_forest.md)
to plot the results.

Other taxa_abundance:
[`sample_sums()`](https://cmmr.github.io/rbiom/reference/sample_sums.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_clusters()`](https://cmmr.github.io/rbiom/reference/taxa_clusters.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_forest()`](https://cmmr.github.io/rbiom/reference/taxa_forest.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md),
[`taxa_sums()`](https://cmmr.github.io/rbiom/reference/taxa_sums.md),
[`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)

Other stats_tables:
[`adiv_stats()`](https://cmmr.github.io/rbiom/reference/adiv_stats.md),
[`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md),
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md),
[`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md)

## Examples

``` r
    library(rbiom)

    biom <- subset(hmp50, `Body Site` %in% c("Saliva", "Stool"))

    da <- taxa_da(biom, stat.by = "Body Site", rank = "Genus")
    da
#> # Method:   LinDA (Zhou et al. 2022)
#> # Model:    log abundance ~ `Body Site`, per taxon
#> # p.adj:    FDR (Benjamini-Hochberg), across taxa within each comparison
#> # A tibble: 147 × 11
#>    `Body Site`  .taxa  .lfc   .se .lower .upper .t.ratio   .df   .p.val   .adj.p
#>    <fct>        <chr> <dbl> <dbl>  <dbl>  <dbl>    <dbl> <dbl>    <dbl>    <dbl>
#>  1 Stool - Sal… Veil… -9.40 0.393 -10.2   -8.57    -23.9    18 4.35e-15 6.40e-13
#>  2 Stool - Sal… Orib… -7.83 0.342  -8.55  -7.11    -22.9    18 9.41e-15 6.92e-13
#>  3 Stool - Sal… Prev… -7.69 0.394  -8.52  -6.86    -19.5    18 1.44e-13 7.05e-12
#>  4 Stool - Sal… Capn… -6.82 0.356  -7.57  -6.07    -19.1    18 2.06e-13 7.56e-12
#>  5 Stool - Sal… Gran… -6.46 0.347  -7.19  -5.73    -18.6    18 3.38e-13 9.92e-12
#>  6 Stool - Sal… Stre… -8.73 0.518  -9.82  -7.64    -16.8    18 1.83e-12 4.49e-11
#>  7 Stool - Sal… Alis…  7.29 0.451   6.34   8.24     16.2    18 3.71e-12 7.80e-11
#>  8 Stool - Sal… Berg… -5.29 0.335  -5.99  -4.58    -15.8    18 5.54e-12 1.01e-10
#>  9 Stool - Sal… Geme… -5.62 0.359  -6.38  -4.87    -15.7    18 6.17e-12 1.01e-10
#> 10 Stool - Sal… Acti… -7.02 0.454  -7.97  -6.06    -15.5    18 7.81e-12 1.08e-10
#> # ℹ 137 more rows
#> # ℹ 1 more variable: .n <int>

    # Adjusting for a covariate.
    taxa_da(biom, regr = "Age", covariates = "Sex", rank = "Phylum")
#> # Method:   LinDA (Zhou et al. 2022)
#> # Model:    log abundance ~ Age + Sex, per taxon
#> # p.adj:    FDR (Benjamini-Hochberg), across taxa within each comparison
#> # A tibble: 11 × 10
#>    .taxa           .lfc    .se  .lower .upper .t.ratio   .df .p.val .adj.p    .n
#>    <chr>          <dbl>  <dbl>   <dbl>  <dbl>    <dbl> <dbl>  <dbl>  <dbl> <int>
#>  1 Synergiste…  1.66e-1 0.0681  0.0227  0.310  2.44e+0    17 0.0258  0.284    20
#>  2 Spirochaet…  1.82e-1 0.0942 -0.0164  0.381  1.94e+0    17 0.0696  0.383    20
#>  3 Fusobacter…  1.95e-1 0.210  -0.247   0.637  9.32e-1    17 0.364   0.999    20
#>  4 Firmicutes   3.67e-2 0.0653 -0.101   0.174  5.62e-1    17 0.582   0.999    20
#>  5 Actinobact…  1.00e-1 0.182  -0.284   0.484  5.50e-1    17 0.590   0.999    20
#>  6 Bacteroide…  3.57e-2 0.139  -0.258   0.330  2.56e-1    17 0.801   0.999    20
#>  7 Sacchariba… -2.30e-2 0.0923 -0.218   0.172 -2.49e-1    17 0.806   0.999    20
#>  8 Verrucomic… -3.93e-2 0.175  -0.408   0.329 -2.25e-1    17 0.825   0.999    20
#>  9 Gracilibac… -1.15e-2 0.0753 -0.170   0.147 -1.52e-1    17 0.881   0.999    20
#> 10 Proteobact… -8.33e-3 0.139  -0.301   0.284 -6.01e-2    17 0.953   0.999    20
#> 11 Tenericutes  1.45e-4 0.175  -0.370   0.370  8.24e-4    17 0.999   0.999    20

    # See the R code that reproduces the table.
    da$code
#> # Generated by rbiom 4.0.0.9001
#> library(rbiom)
#> 
#> biom <- hmp50
#> biom <- subset(biom, `Body Site` %in% c("Saliva", "Stool"))
#> 
#> # Differential abundance with LinDA (Zhou et al. 2022), as implemented in
#> # rbiom; see ?taxa_da for the algorithm. For each rank, the counts from
#> # taxa_matrix() are modeled per taxon as centered log2-ratios with
#> #   log abundance ~ `Body Site`
#> # and each coefficient is corrected by its mode across taxa.
#> stats <- taxa_da(biom, "Body Site", rank = "Genus")
```
