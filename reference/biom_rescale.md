# Rescale Counts to a Specific Range

This function performs a min-max scaling on each sample independently.

It is useful for normalization techniques that require data to be within
a specific bounded range, or for visualization purposes where
maintaining the relative distances between values is important but the
absolute magnitude needs adjustment.

## Usage

``` r
biom_rescale(biom, range = c(0, 1), clone = TRUE)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- range:

  Numeric vector of length 2. Target min and max. Default: `c(0, 1)`.

- clone:

  Create a copy of `biom` before modifying. If `FALSE`, `biom` is
  modified in place as a side-effect. See [speed
  ups](https://cmmr.github.io/rbiom/reference/speed.md) for use cases.
  Default: `TRUE`

## Value

An [rbiom
object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md).

## Details

Linearly rescale each sample's values to lie between a specified minimum
and maximum.

## Note

If `range` starts at a non-zero value (e.g., `c(1, 10)`), the sparsity
of the matrix will be destroyed because all zero counts will be shifted
to the minimum value. This can significantly increase memory usage for
large datasets.

## Mathematical Transformation

The rescaling is performed in two steps:

1.  **Normalize:** Divide values by the maximum value in that sample,
    scaling them to a `[0, 1]` range relative to the sample's peak.

2.  **Scale and Shift:** Apply the target range using the formula:
    \$\$x\_{new} = x\_{norm} \times (max - min) + min\$\$

## See also

Other transformations:
[`biom_inflate()`](https://cmmr.github.io/rbiom/reference/biom_inflate.md),
[`biom_relativize()`](https://cmmr.github.io/rbiom/reference/biom_relativize.md),
[`modify_metadata`](https://cmmr.github.io/rbiom/reference/modify_metadata.md),
[`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md),
[`slice_metadata`](https://cmmr.github.io/rbiom/reference/slice_metadata.md),
[`subset()`](https://cmmr.github.io/rbiom/reference/subset.md),
[`with()`](https://cmmr.github.io/rbiom/reference/with.md)

## Examples

``` r
    library(rbiom)
    
    biom <- hmp50[1:5]
    
    # Original range
    range(as.matrix(biom))
#> [1]    0 2672
    
    # Rescaled to 0-1
    biom_01 <- biom_rescale(biom)
    range(as.matrix(biom_01))
#> [1] 0 1
    
    # Rescaled to 0-100 (Percentages)
    biom_100 <- biom_rescale(biom, range = c(0, 100))
    range(as.matrix(biom_100))
#> [1]   0 100
    
```
