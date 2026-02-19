# Relativize Counts to Proportions

This function normalizes the data by dividing each observation by the
total library size of its sample. The resulting values represent the
proportion (0 to 1) of the sample composed of that specific feature.

This is a common transformation for microbiome data, as it accounts for
differences in sequencing depth across samples, allowing for comparison
of community composition.

## Usage

``` r
biom_relativize(biom, clone = TRUE)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- clone:

  Create a copy of `biom` before modifying. If `FALSE`, `biom` is
  modified in place as a side-effect. See [speed
  ups](https://cmmr.github.io/rbiom/reference/speed.md) for use cases.
  Default: `TRUE`

## Value

An [rbiom
object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md).

## Details

Convert absolute counts to relative abundances (proportions) where each
sample sums to 1.

## See also

Other transformations:
[`biom_inflate()`](https://cmmr.github.io/rbiom/reference/biom_inflate.md),
[`biom_rescale()`](https://cmmr.github.io/rbiom/reference/biom_rescale.md),
[`modify_metadata`](https://cmmr.github.io/rbiom/reference/modify_metadata.md),
[`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md),
[`slice_metadata`](https://cmmr.github.io/rbiom/reference/slice_metadata.md),
[`subset()`](https://cmmr.github.io/rbiom/reference/subset.md),
[`with()`](https://cmmr.github.io/rbiom/reference/with.md)

## Examples

``` r
    library(rbiom)
    
    biom <- hmp50[1:5]
    
    # Raw counts sum to different library sizes
    sample_sums(biom)
#> HMP01 HMP02 HMP03 HMP04 HMP05 
#>  1660  1371  1353  1895  3939 
    
    # Relativized counts sum to 1
    biom_rel <- biom_relativize(biom)
    sample_sums(biom_rel)
#> HMP01 HMP02 HMP03 HMP04 HMP05 
#>     1     1     1     1     1 
```
