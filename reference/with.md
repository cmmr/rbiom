# Evaluate expressions on metadata.

`with()` will return the result of your expression.
[`within()`](https://rdrr.io/r/base/with.html) will return an rbiom
object.

## Usage

``` r
# S3 method for class 'rbiom'
with(data, expr, ...)

# S3 method for class 'rbiom'
within(data, expr, clone = TRUE, ...)
```

## Arguments

- data:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), such
  as from
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- expr:

  Passed on to [`base::with()`](https://rdrr.io/r/base/with.html) or
  [`base::within()`](https://rdrr.io/r/base/with.html).

- ...:

  Not used.

- clone:

  Create a copy of `biom` before modifying. If `FALSE`, `biom` is
  modified in place as a side-effect. See [speed
  ups](https://cmmr.github.io/rbiom/reference/speed.md) for use cases.
  Default: `TRUE`

## Value

See description.

## See also

Other transformations:
[`matrix_ops`](https://cmmr.github.io/rbiom/reference/matrix_ops.md),
[`modify_metadata`](https://cmmr.github.io/rbiom/reference/modify_metadata.md),
[`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md),
[`slice_metadata`](https://cmmr.github.io/rbiom/reference/slice_metadata.md),
[`subset()`](https://cmmr.github.io/rbiom/reference/subset.md)

## Examples

``` r
    library(rbiom) 
    
    with(hmp50, table(`Body Site`, Sex))
#>                 Sex
#> Body Site        Female Male
#>   Anterior nares      5    5
#>   Buccal mucosa       5    5
#>   Mid vagina         10    0
#>   Saliva              5    5
#>   Stool               5    5
    
    biom <- within(hmp50, {
      age_bin = cut(Age, 5)
      bmi_bin = cut(BMI, 5)
    })
    biom$metadata
#> # A tibble: 50 × 7
#>    .sample   Age   BMI `Body Site`    Sex    bmi_bin     age_bin    
#>  * <chr>   <dbl> <dbl> <fct>          <fct>  <fct>       <fct>      
#>  1 HMP01      22    20 Buccal mucosa  Female (19,21.6]   (21,24.8]  
#>  2 HMP02      24    23 Buccal mucosa  Male   (21.6,24.2] (21,24.8]  
#>  3 HMP03      28    26 Saliva         Male   (24.2,26.8] (24.8,28.6]
#>  4 HMP04      25    23 Saliva         Male   (21.6,24.2] (24.8,28.6]
#>  5 HMP05      27    24 Buccal mucosa  Female (21.6,24.2] (24.8,28.6]
#>  6 HMP06      32    25 Saliva         Male   (24.2,26.8] (28.6,32.4]
#>  7 HMP07      26    22 Buccal mucosa  Male   (21.6,24.2] (24.8,28.6]
#>  8 HMP08      27    26 Saliva         Female (24.2,26.8] (24.8,28.6]
#>  9 HMP09      33    32 Saliva         Male   (29.4,32]   (32.4,36.2]
#> 10 HMP10      22    20 Anterior nares Female (19,21.6]   (21,24.8]  
#> # ℹ 40 more rows
```
