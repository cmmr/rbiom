# Get a glimpse of your metadata.

Get a glimpse of your metadata.

## Usage

``` r
# S3 method for class 'rbiom'
glimpse(x, width = NULL, ...)
```

## Arguments

- x:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), such
  as from
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- width:

  Width of output. See
  [`pillar::glimpse()`](https://pillar.r-lib.org/reference/glimpse.html)
  documentation. Default: `NULL`

- ...:

  Not used.

## Value

The original `biom`, invisibly.

## See also

Other metadata:
[`bdply()`](https://cmmr.github.io/rbiom/reference/bdply.md)

## Examples

``` r
    library(rbiom)
    
    glimpse(hmp50)
#> Rows: 50
#> Columns: 5
#> $ .sample     <chr> "HMP01", "HMP02", "HMP03", "HMP04", "HMP05", "HMP06", "HMP…
#> $ Age         <dbl> 22, 24, 28, 25, 27, 32, 26, 27, 33, 22, 24, 35, 24, 32, 25…
#> $ BMI         <dbl> 20, 23, 26, 23, 24, 25, 22, 26, 32, 20, 23, 26, 21, 26, 21…
#> $ `Body Site` <fct> Buccal mucosa, Buccal mucosa, Saliva, Saliva, Buccal mucos…
#> $ Sex         <fct> Female, Male, Male, Male, Female, Male, Male, Female, Male…
    
```
