# Rarefy OTU counts.

Sub-sample OTU observations such that all samples have an equal number.

## Usage

``` r
rarefy(
  biom,
  depth = 0.1,
  n = NULL,
  seed = 0L,
  upsample = NULL,
  clone = TRUE,
  cpus = NULL
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- depth:

  How many observations to keep per sample. When `0 < depth < 1`, it is
  taken as the minimum percentage of the dataset's observations to keep.
  Ignored when `n` is specified. Default: `0.1`

- n:

  The number of samples to keep. When `0 < n < 1`, it is taken as the
  percentage of samples to keep. If negative, that number or percentage
  of samples is dropped. If `0`, all samples are kept. If `NULL`,
  `depth` is used instead. Default: `NULL`

- seed:

  An integer seed for randomizing which observations to keep or drop. If
  you need to create different random rarefactions of the same data, set
  the seed to a different number each time.

- upsample:

  If the count data is in percentages, provide an integer value here to
  scale each sample's observations to integers that sum to this value.
  Generally not recommended, but can be used to 'shoehorn' metagenomic
  abundance estimates into rbiom's functions that were designed for
  amplicon datasets. When invoked, `depth`, `n`, and `seed` are ignored.
  The default, `NULL`, will throw an error if the counts are not all
  integers.

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

## See also

Other transformations:
[`matrix_ops`](https://cmmr.github.io/rbiom/reference/matrix_ops.md),
[`modify_metadata`](https://cmmr.github.io/rbiom/reference/modify_metadata.md),
[`slice_metadata`](https://cmmr.github.io/rbiom/reference/slice_metadata.md),
[`subset()`](https://cmmr.github.io/rbiom/reference/subset.md),
[`with()`](https://cmmr.github.io/rbiom/reference/with.md)

## Examples

``` r
    library(rbiom)
    
    sample_sums(hmp50) %>% head()
#> HMP01 HMP02 HMP03 HMP04 HMP05 HMP06 
#>  1660  1371  1353  1895  3939  4150 
    
    biom <- rarefy(hmp50)
    sample_sums(biom) %>% head()
#> HMP01 HMP02 HMP03 HMP04 HMP05 HMP06 
#>  1183  1183  1183  1183  1183  1183 
```
