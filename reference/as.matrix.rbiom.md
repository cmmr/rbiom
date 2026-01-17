# Convert an rbiom object to a simple count matrix.

Identical to running `as.matrix(biom$counts)`.

## Usage

``` r
# S3 method for class 'rbiom'
as.matrix(x, ...)
```

## Arguments

- x:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), such
  as from
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- ...:

  Not used.

## Value

A base R matrix with OTUs as rows and samples as columns.

## See also

Other conversion:
[`as.list.rbiom()`](https://cmmr.github.io/rbiom/reference/as.list.rbiom.md)

## Examples

``` r
    library(rbiom)
    
    as.matrix(hmp50)[1:5,1:5]
#>          HMP01 HMP02 HMP03 HMP04 HMP05
#> Unc01yki     0     0     0     0     0
#> Unc53100  1083   543   301   223  2672
#> LtbAci52     0     0     0     0     0
#> CnbTube3     0     0     0     0     0
#> Unc02qsf     0     0     0     0     0
```
