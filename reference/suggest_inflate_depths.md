# Suggest Inflation Depths

Estimates the optimal sequencing depth for each sample in a matrix by
leveraging the global abundance distribution structure.

## Usage

``` r
suggest_inflate_depths(biom, adjust = 1.5)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- adjust:

  Numeric. Bandwidth adjustment for the kernel density estimation.
  Default: `1.5`.

## Value

A named integer vector of recommended depths for each sample.

## Exact Depths

When `depth = NULL`,
[`biom_inflate()`](https://cmmr.github.io/rbiom/reference/biom_inflate.md)
calls this function to find the original sequencing depth of each
sample. Values on any scale are first divided by their sample's total,
so proportions, percentages, and counts divided by a constant are all
handled alike.

If every non-zero value in a sample is a whole multiple of the sample's
smallest non-zero value, that value is taken to be one read, and the
depth is the sample's total divided by it. This recovers the original
counts exactly from counts divided by a constant, and from proportions
of integer counts whenever the sample contains a singleton. Other
samples get the heuristic below.

## The Singleton Peak Heuristic

The underlying assumption is that in typical microbiome datasets, the
most frequent count value (the mode of the abundance distribution) is 1
(a singleton).

The algorithm works as follows:

1.  **Log-Transformation:** Non-zero relative abundances are
    log10-transformed.

2.  **Global Consensus:** To overcome sparsity in individual samples,
    distributions are centered by their medians and aggregated across
    all samples.

3.  **Peak Detection:** Kernel Density Estimation (KDE) is used to
    identify the peak (mode) of this aggregated distribution.

4.  **Scaling:** A scaling factor is calculated for each sample that
    shifts this peak to correspond to an integer count of 1.

This approach effectively "shoehorns" relative abundance data into
integer formats required by diversity metrics (like rarefaction or
Chao1) by maximizing the number of singletons in the resulting matrix.

Samples with fewer than two non-zero features carry no information about
the distribution's shape; they are assigned the median depth of the
other samples (or, if there are too few informative samples for the
heuristic, every sample is scaled so that its smallest non-zero value
becomes 1). Counts that are already integers are returned as their
current sums.

**Accuracy:** the heuristic is an estimate, not a recovery. On `hmp50`
proportions perturbed so that no exact depth can be found, its depths
ranged from 0.06 to 2.9 times the true depths (median 0.72), so the
inflated counts, and any rarefaction depth chosen from them, can differ
substantially from the original data. Rounded proportions (e.g. to 6
decimal places) are in this situation too. Pass `depth` to
[`biom_inflate()`](https://cmmr.github.io/rbiom/reference/biom_inflate.md)
when the original depths are known.

## See also

[`biom_inflate()`](https://cmmr.github.io/rbiom/reference/biom_inflate.md)
which uses this heuristic when `depth = NULL`.

## Examples

``` r
    library(rbiom)

    depths <- suggest_inflate_depths(biom_relativize(hmp50))
    head(depths)
#> HMP01 HMP02 HMP03 HMP04 HMP05 HMP06 
#>  1660  1371  1353  1895  3939  4150 
```
