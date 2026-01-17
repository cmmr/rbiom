# Map sample names to metadata field values.

Map sample names to metadata field values.

## Usage

``` r
# S3 method for class 'rbiom'
pull(.data, var = -1, name = ".sample", ...)
```

## Arguments

- .data:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), such
  as from
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- var:

  The metadata field name specified as:

  - The metadata field name to retrieve. Can be abbreviated.

  - A positive integer, giving the position counting from the left.

  - A negative integer, giving the position counting from the right.

  Default: `-1`

- name:

  The column to be used as names for a named vector. Specified in a
  similar manner as var. Default: `".sample"`

- ...:

  Not used.

## Value

A vector of metadata values, named with sample names.

## See also

[`taxa_map()`](https://cmmr.github.io/rbiom/reference/taxa_map.md)

Other samples:
[`sample_sums()`](https://cmmr.github.io/rbiom/reference/sample_sums.md)

## Examples

``` r
    library(rbiom)
    
    pull(hmp50, 'Age') %>% head()
#> HMP01 HMP02 HMP03 HMP04 HMP05 HMP06 
#>    22    24    28    25    27    32 
    
    pull(hmp50, 'bod') %>% head(4)
#>         HMP01         HMP02         HMP03         HMP04 
#> Buccal mucosa Buccal mucosa        Saliva        Saliva 
#> Levels: Anterior nares Buccal mucosa Mid vagina Saliva Stool
    
```
