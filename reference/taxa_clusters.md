# Cluster samples by taxa abundances (k-means).

Runs [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) on the
samples x taxa abundance matrix, with `nstart = 10` random starts so
that the partition is stable, and with the random number generator
seeded by `seed` so that it is reproducible.

## Usage

``` r
taxa_clusters(biom, rank = ".otu", k = 5, seed = 0, ...)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- rank:

  Which taxa rank to use. E.g. `"Phylum"`, `"Genus"`, `".otu"`, etc. An
  integer can also be given, where `1` is the highest rank, `2` is the
  second highest, `-1` is the lowest rank, `-2` is the second lowest,
  and `0` is the OTU "rank". Run `biom$ranks` to see all options for a
  given rbiom object. Default: `.otu`.

- k:

  Number of clusters. Default: `5L`

- seed:

  Random seed for the k-means starting centers. Default: `0`

- ...:

  Passed on to [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html).

## Value

A numeric factor assigning samples to clusters.

## See also

Other taxa_abundance:
[`sample_sums()`](https://cmmr.github.io/rbiom/reference/sample_sums.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md),
[`taxa_sums()`](https://cmmr.github.io/rbiom/reference/taxa_sums.md),
[`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)

Other clustering:
[`bdiv_clusters()`](https://cmmr.github.io/rbiom/reference/bdiv_clusters.md)

## Examples

``` r
    library(rbiom)

    biom <- rarefy(hmp50)
    biom$metadata$otu_cluster <- taxa_clusters(biom)

    pull(biom, 'otu_cluster')[1:10]
#> HMP01 HMP02 HMP03 HMP04 HMP05 HMP06 HMP07 HMP08 HMP09 HMP10 
#>     3     3     2     2     3     2     3     2     2     4 
#> Levels: 1 2 3 4 5

    bdiv_ord_plot(biom, layers = "p", stat.by = "otu_cluster")
```
