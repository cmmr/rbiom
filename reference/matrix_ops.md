# Deprecated matrix transformations

A collection of transformations that operate directly on matrices.

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

  A numeric matrix or sparse matrix of counts.

- margin:

  Apply the transformation to the matrix's rows (`margin=1L`) or columns
  (`margin=2L`). Instead of `1L` and `2L`, you may also use `'rows'` and
  `'cols'`. Default: `2L` (column-wise, aka sample-wise for otu tables)

- depth:

  How many observations to keep per sample.

- n:

  Deprecated. The number of samples to keep. This argument is ignored in
  the current version.

- seed:

  An integer seed for randomizing which observations to keep or drop.

- upsample:

  If the count data is in percentages, provide an integer value here to
  scale each sample's observations to integers that sum to this value.
  Maps to `inflate` in the new syntax.

- cpus:

  The number of CPUs to use. Set to `NULL` to use all available, or to
  `1` to disable parallel processing. Default: `NULL`

- range:

  When rescaling, what should the minimum and maximum values be?
  Default: `c(0, 1)`

## Value

A A numeric matrix or sparse matrix, depending on the input type, with
the same dimensions as `mtx`.

## Examples

``` r
    mtx <- matrix(1:12, ncol = 3, dimnames = list(paste0("OTU", 1:4), paste0("Sample", 1:3)))
    
    mtx
#>      Sample1 Sample2 Sample3
#> OTU1       1       5       9
#> OTU2       2       6      10
#> OTU3       3       7      11
#> OTU4       4       8      12
    
    suppressWarnings({
    
      mtx_rarefy(mtx)
      
      rarefy_cols(mtx)
      
      mtx_percent(mtx)
      
      mtx_rescale(mtx)
      
      rescale_rows(mtx)
      
      rescale_cols(mtx)
      
    })
#>      Sample1   Sample2   Sample3
#> OTU1     0.1 0.1923077 0.2142857
#> OTU2     0.2 0.2307692 0.2380952
#> OTU3     0.3 0.2692308 0.2619048
#> OTU4     0.4 0.3076923 0.2857143
    
```
