# Summarize the taxa observations in each sample.

Summarize the taxa observations in each sample.

## Usage

``` r
sample_sums(biom, rank = -1, sort = NULL, unc = "singly")

sample_apply(biom, FUN, rank = -1, sort = NULL, unc = "singly", ...)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- rank:

  What rank(s) of taxa to display. E.g. `"Phylum"`, `"Genus"`, `".otu"`,
  etc. An integer vector can also be given, where `1` is the highest
  rank, `2` is the second highest, `-1` is the lowest rank, `-2` is the
  second lowest, and `0` is the OTU "rank". Run `biom$ranks` to see all
  options for a given rbiom object. Default: `-1`.

- sort:

  Sort the result. Options: `NULL` - don't sort; `"asc"` - in ascending
  order (smallest to largest); `"desc"` - in descending order (largest
  to smallest). Ignored when the result is not a simple numeric vector.
  Default: `NULL`

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

- FUN:

  The function to apply to each column of
  [`taxa_matrix()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md).

- ...:

  Optional arguments to `FUN`.

## Value

For `sample_sums`, A named numeric vector of the number of observations
in each sample. For `sample_apply`, a named vector or list with the
results of `FUN`. The names are the taxa IDs.

## See also

Other samples:
[`pull.rbiom()`](https://cmmr.github.io/rbiom/reference/pull.rbiom.md)

Other rarefaction:
[`rare_corrplot()`](https://cmmr.github.io/rbiom/reference/rare_corrplot.md),
[`rare_multiplot()`](https://cmmr.github.io/rbiom/reference/rare_multiplot.md),
[`rare_stacked()`](https://cmmr.github.io/rbiom/reference/rare_stacked.md)

Other taxa_abundance:
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_clusters()`](https://cmmr.github.io/rbiom/reference/taxa_clusters.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md),
[`taxa_sums()`](https://cmmr.github.io/rbiom/reference/taxa_sums.md),
[`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)

## Examples

``` r
    library(rbiom)
    library(ggplot2)
    
    sample_sums(hmp50, sort = 'asc') %>% head()
#> HMP36 HMP24 HMP03 HMP02 HMP42 HMP17 
#>   182  1183  1353  1371  1489  1579 
    
    # Unique OTUs and "cultured" classes per sample
    nnz <- function (x) sum(x > 0) # number of non-zeroes
    sample_apply(hmp50, nnz, 'otu') %>% head()
#> HMP01 HMP02 HMP03 HMP04 HMP05 HMP06 
#>    49    75    75    83    67   105 
    sample_apply(hmp50, nnz, 'class', unc = 'drop') %>% head()
#> HMP01 HMP02 HMP03 HMP04 HMP05 HMP06 
#>    10    13    12    13    12    15 
    
    # Number of reads in each sample's most abundant family
    sample_apply(hmp50, base::max, 'f', sort = 'desc') %>% head()
#> HMP44 HMP25 HMP11 HMP21 HMP34 HMP46 
#> 16220  9581  6308  5786  4645  4050 
    
    ggplot() + geom_histogram(aes(x=sample_sums(hmp50)), bins = 20)
```
