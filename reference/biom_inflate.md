# Inflate Relative Abundances to Integer Counts

Scaling a matrix of proportions (or counts) to a new target depth,
rounding to integers while preserving the original total abundance sum
exactly.

## Usage

``` r
biom_inflate(biom, depth = NULL, clone = TRUE)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- depth:

  The target library size (sum) for each sample. Must be an integer
  greater than 0. If `NULL` (the default), the depth is estimated
  per-sample using the "Singleton Peak Heuristic". See
  [`suggest_inflate_depths()`](https://cmmr.github.io/rbiom/reference/suggest_inflate_depths.md)
  for algorithm details.

- clone:

  Create a copy of `biom` before modifying. If `FALSE`, `biom` is
  modified in place as a side-effect. See [speed
  ups](https://cmmr.github.io/rbiom/reference/speed.md) for use cases.
  Default: `TRUE`

## Value

An [rbiom
object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md).

## Rounding (Largest Remainder Method)

To ensure the sum of the resulting counts equals the target `depth`
exactly (avoiding drift caused by simple rounding), this function uses
the Largest Remainder Method (also known as the Hare-Niemeyer method).

It assigns the integer part of the scaled value to each feature, and
then distributes the remaining counts to the features with the largest
fractional parts.

## See also

[`suggest_inflate_depths()`](https://cmmr.github.io/rbiom/reference/suggest_inflate_depths.md)
for details on how target depths are estimated when `depth = NULL`.

Other transformations:
[`biom_relativize()`](https://cmmr.github.io/rbiom/reference/biom_relativize.md),
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
    sample_sums(biom)
#> HMP01 HMP02 HMP03 HMP04 HMP05 
#>  1660  1371  1353  1895  3939 
    
    biom <- biom_relativize(biom)
    sample_sums(biom)
#> HMP01 HMP02 HMP03 HMP04 HMP05 
#>     1     1     1     1     1 
    
    biom <- biom_inflate(biom)
    sample_sums(biom)
#> HMP01 HMP02 HMP03 HMP04 HMP05 
#>  1267  1395  1032  1157  4007 
```
