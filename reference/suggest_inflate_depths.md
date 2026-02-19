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

## The Singleton Peak Heuristic

When `depth = NULL`,
[`biom_inflate()`](https://cmmr.github.io/rbiom/reference/biom_inflate.md)
calls this function to estimate the original sequencing depth for each
sample. The underlying assumption is that in typical microbiome
datasets, the most frequent count value (the mode of the abundance
distribution) is 1 (a singleton).

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

## See also

[`biom_inflate()`](https://cmmr.github.io/rbiom/reference/biom_inflate.md)
which uses this heuristic when `depth = NULL`.

## Examples

``` r
    library(rbiom)
    
    depths <- suggest_inflate_depths(hmp50)
    head(depths)
#> HMP01 HMP02 HMP03 HMP04 HMP05 HMP06 
#>     1     1     1     1     1     0 
```
