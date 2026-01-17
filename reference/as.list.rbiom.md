# Convert an rbiom object to a base R list.

Convert an rbiom object to a base R list.

## Usage

``` r
# S3 method for class 'rbiom'
as.list(x, ...)
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

A list with names
`c('counts', 'metadata', 'taxonomy', 'tree', 'sequences', 'id', 'comment', 'date', 'generated_by')`.

## See also

Other conversion:
[`as.matrix.rbiom()`](https://cmmr.github.io/rbiom/reference/as.matrix.rbiom.md)
