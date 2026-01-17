# Transform a counts matrix.

A collection of transformations that operate directly on matrices.  
  
Note: `rarefy_cols()`, `rescale_rows()`, and `rescale_cols()` are
deprecated.

## Usage

``` r
mtx_rarefy(
  mtx,
  margin = 2L,
  depth = 0.1,
  n = NULL,
  seed = 0L,
  upsample = NULL,
  cpus = NULL
)

mtx_percent(mtx, margin = 2L)

mtx_rescale(mtx, margin = 2L, range = c(0, 1))

rarefy_cols(mtx, depth = 0.1, n = NULL, seed = 0L, cpus = NULL)

rescale_rows(mtx)

rescale_cols(mtx)
```

## Arguments

- mtx:

  A matrix-like object.

- margin:

  Apply the transformation to the matrix's rows (`margin=1L`) or columns
  (`margin=2L`). Instead of `1L` and `2L`, you may also use `'rows'` and
  `'cols'`. Default: `2L` (column-wise, aka sample-wise for otu tables)

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

  Random seed for permutations. Must be a non-negative integer. Default:
  `0`

- upsample:

  If the count data is in percentages, provide an integer value here to
  scale each sample's observations to integers that sum to this value.
  Generally not recommended, but can be used to 'shoehorn' metagenomic
  abundance estimates into rbiom's functions that were designed for
  amplicon datasets. When invoked, `depth`, `n`, and `seed` are ignored.
  The default, `NULL`, will throw an error if the counts are not all
  integers.

- cpus:

  The number of CPUs to use. Set to `NULL` to use all available, or to
  `1` to disable parallel processing. Default: `NULL`

- range:

  When rescaling, what should the minimum and maximum values be?
  Default: `c(0, 1)`

## Value

The transformed matrix. If `mtx` was a sparse matrix from the `Matrix`
package, then the result will also be a sparse matrix, otherwise the
result will be a base R matrix.

## See also

Other transformations:
[`modify_metadata`](https://cmmr.github.io/rbiom/reference/modify_metadata.md),
[`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md),
[`slice_metadata`](https://cmmr.github.io/rbiom/reference/slice_metadata.md),
[`subset()`](https://cmmr.github.io/rbiom/reference/subset.md),
[`with()`](https://cmmr.github.io/rbiom/reference/with.md)

## Examples

``` r
    library(rbiom)
    
    # mtx_rarefy --------------------------------------
    biom <- hmp50$clone()
    sample_sums(biom) %>% head(10)
#> HMP01 HMP02 HMP03 HMP04 HMP05 HMP06 HMP07 HMP08 HMP09 HMP10 
#>  1660  1371  1353  1895  3939  4150  3283  1695  2069  2509 

    biom$counts %<>% mtx_rarefy(depth=1000)
    sample_sums(biom) %>% head(10)
#> HMP01 HMP02 HMP03 HMP04 HMP05 HMP06 HMP07 HMP08 HMP09 HMP10 
#>  1000  1000  1000  1000  1000  1000  1000  1000  1000  1000 
    
    
    # rescaling ----------------------------------------
    mtx <- matrix(sample(1:20), nrow=4)
    mtx
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   10    7    3    4    8
#> [2,]   12   15    5   13   16
#> [3,]    2   20    6   19   14
#> [4,]   18    9   17    1   11
    
    colSums(mtx)
#> [1] 42 51 31 37 49
    
    colSums(mtx_rarefy(mtx))
#> [1] 31 31 31 31 31
    
    colSums(mtx_percent(mtx))
#> [1] 1 1 1 1 1
    
    apply(mtx_rescale(mtx), 2L, max)
#> [1] 1 1 1 1 1
```
