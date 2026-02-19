# Create, modify, and delete metadata fields.

mutate() creates new fields in `$metadata` that are functions of
existing metadata fields. It can also modify (if the name is the same as
an existing field) and delete fields (by setting their value to NULL).

## Usage

``` r
# S3 method for class 'rbiom'
mutate(.data, ..., clone = TRUE)

# S3 method for class 'rbiom'
rename(.data, ..., clone = TRUE)
```

## Arguments

- .data:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), such
  as from
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- ...:

  Passed on to
  [`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
  or
  [`dplyr::rename()`](https://dplyr.tidyverse.org/reference/rename.html).

- clone:

  Create a copy of `biom` before modifying. If `FALSE`, `biom` is
  modified in place as a side-effect. See [speed
  ups](https://cmmr.github.io/rbiom/reference/speed.md) for use cases.
  Default: `TRUE`

## Value

An [rbiom
object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md).

## See also

Other transformations:
[`biom_inflate()`](https://cmmr.github.io/rbiom/reference/biom_inflate.md),
[`biom_relativize()`](https://cmmr.github.io/rbiom/reference/biom_relativize.md),
[`biom_rescale()`](https://cmmr.github.io/rbiom/reference/biom_rescale.md),
[`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md),
[`slice_metadata`](https://cmmr.github.io/rbiom/reference/slice_metadata.md),
[`subset()`](https://cmmr.github.io/rbiom/reference/subset.md),
[`with()`](https://cmmr.github.io/rbiom/reference/with.md)

## Examples

``` r
    library(rbiom) 
    
    biom <- slice_max(hmp50, BMI, n = 6)
    biom$metadata
#> # A tibble: 7 × 5
#>   .sample   Age   BMI `Body Site`    Sex   
#> * <chr>   <dbl> <dbl> <fct>          <fct> 
#> 1 HMP09      33    32 Saliva         Male  
#> 2 HMP25      33    32 Anterior nares Male  
#> 3 HMP37      38    31 Mid vagina     Female
#> 4 HMP27      24    30 Saliva         Female
#> 5 HMP43      24    30 Mid vagina     Female
#> 6 HMP36      24    29 Stool          Male  
#> 7 HMP48      24    29 Anterior nares Male  
    
    # Add a new field to the metadata
    biom <- mutate(biom, Obsese = BMI >= 30)
    biom$metadata
#> # A tibble: 7 × 6
#>   .sample   Age   BMI `Body Site`    Sex    Obsese
#> * <chr>   <dbl> <dbl> <fct>          <fct>  <lgl> 
#> 1 HMP09      33    32 Saliva         Male   TRUE  
#> 2 HMP25      33    32 Anterior nares Male   TRUE  
#> 3 HMP37      38    31 Mid vagina     Female TRUE  
#> 4 HMP27      24    30 Saliva         Female TRUE  
#> 5 HMP43      24    30 Mid vagina     Female TRUE  
#> 6 HMP36      24    29 Stool          Male   FALSE 
#> 7 HMP48      24    29 Anterior nares Male   FALSE 
    
    # Rename a metadata field
    biom <- rename(biom, 'Age (years)' = "Age")
    biom$metadata
#> # A tibble: 7 × 6
#>   .sample `Age (years)`   BMI `Body Site`    Sex    Obsese
#> * <chr>           <dbl> <dbl> <fct>          <fct>  <lgl> 
#> 1 HMP09              33    32 Saliva         Male   TRUE  
#> 2 HMP25              33    32 Anterior nares Male   TRUE  
#> 3 HMP37              38    31 Mid vagina     Female TRUE  
#> 4 HMP27              24    30 Saliva         Female TRUE  
#> 5 HMP43              24    30 Mid vagina     Female TRUE  
#> 6 HMP36              24    29 Stool          Male   FALSE 
#> 7 HMP48              24    29 Anterior nares Male   FALSE 
```
