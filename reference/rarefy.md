# Rarefy Counts to a Constant Depth

This function reduces the number of observations (reads) in each sample
to a fixed integer value (`depth`). Samples with fewer observations than
the specified depth are discarded.

Rarefaction is a common technique in microbiome analysis used to account
for uneven sequencing effort across samples. By standardizing the
library size, it allows for fair comparisons of alpha and beta diversity
metrics.

## Usage

``` r
rarefy(
  biom,
  depth = NULL,
  seed = 0L,
  inflate = FALSE,
  clone = TRUE,
  cpus = n_cpus()
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- depth:

  The number of observations to keep per sample. Must be an integer
  greater than 0.

  - If `NULL` (the default), a depth is automatically selected that
    retains at least 10% of the dataset's total abundance while
    maximizing the number of samples kept. See
    [`suggest_rarefy_depth()`](https://cmmr.github.io/rbiom/reference/suggest_rarefy_depth.md)
    for the specific heuristic used.

  - Samples with total counts less than `depth` will be dropped from the
    result.

- seed:

  Random seed for permutations. Must be a non-negative integer. Default:
  `0`

- inflate:

  Logical. Handling for non-integer data (e.g. relative abundances).

  - `FALSE` (Default): The function will error if non-integers are
    detected. Rarefaction requires discrete counts (integers).

  - `TRUE`: The function will automatically rescale (inflate)
    non-integers to integers using
    [`biom_inflate()`](https://cmmr.github.io/rbiom/reference/biom_inflate.md)
    before rarefying. This is useful for 'shoehorning' metagenomic
    relative abundance data into diversity functions that strictly
    require integers.

- clone:

  Create a copy of `biom` before modifying. If `FALSE`, `biom` is
  modified in place as a side-effect. See [speed
  ups](https://cmmr.github.io/rbiom/reference/speed.md) for use cases.
  Default: `TRUE`

- cpus:

  The number of CPUs to use. Set to `NULL` to use all available, or to
  `1` to disable parallel processing. Default: `NULL`

## Value

An [rbiom
object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md).

## Details

Normalizes the library sizes of a dataset by randomly sub-sampling
observations from each sample to a specific depth.

## See also

[`suggest_rarefy_depth()`](https://cmmr.github.io/rbiom/reference/suggest_rarefy_depth.md)
for details on the default depth selection.

Other transformations:
[`biom_inflate()`](https://cmmr.github.io/rbiom/reference/biom_inflate.md),
[`biom_relativize()`](https://cmmr.github.io/rbiom/reference/biom_relativize.md),
[`biom_rescale()`](https://cmmr.github.io/rbiom/reference/biom_rescale.md),
[`modify_metadata`](https://cmmr.github.io/rbiom/reference/modify_metadata.md),
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
    
    # Rarefy to the lowest sample depth 
    # (All samples are kept, but counts are reduced)
    biom_rare <- rarefy(biom, depth = min(sample_sums(biom)))
    sample_sums(biom_rare)
#> HMP01 HMP02 HMP03 HMP04 HMP05 
#>  1353  1353  1353  1353  1353 
    
    # Auto-select depth (may drop samples with low coverage)
    biom_auto <- rarefy(biom)
    sample_sums(biom_auto)
#> HMP01 HMP02 HMP03 HMP04 HMP05 
#>  1353  1353  1353  1353  1353 
```
