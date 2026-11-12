# List the color palettes.

The palettes that rbiom's `colors`, `grid`, and track `colors` arguments
accept by name. Prefix a name with `"-"` to reverse its colors, e.g.
`"-bilbao"`.

## Usage

``` r
list_palettes(kind = NULL)
```

## Arguments

- kind:

  Which kinds to list: any of `"discrete"`, `"sequential"`, and
  `"diverging"`. `NULL` lists them all. Default: `NULL`

## Value

A tibble with one row per palette: `name`, `kind`, and `colors` (a list
column of hex codes).

## Details

The discrete palettes are for categories. `"okabe"`, `"tol"`,
`"bright"`, `"vibrant"`, `"muted"`, `"light"`, and `"carto"` are
designed for color vision deficiencies. The sequential and diverging
palettes are for numbers: heat map grids and numeric tracks.

## See also

[`fn_schema()`](https://cmmr.github.io/rbiom/reference/fn_schema.md)

## Examples

``` r
    library(rbiom)

    list_palettes()
#> # A tibble: 44 × 3
#>    name    kind     colors    
#>    <chr>   <chr>    <list>    
#>  1 bright  discrete <chr [7]> 
#>  2 vibrant discrete <chr [7]> 
#>  3 okabe   discrete <chr [8]> 
#>  4 r4      discrete <chr [8]> 
#>  5 light   discrete <chr [8]> 
#>  6 muted   discrete <chr [9]> 
#>  7 fishy   discrete <chr [9]> 
#>  8 tableau discrete <chr [10]>
#>  9 classic discrete <chr [10]>
#> 10 carto   discrete <chr [12]>
#> # ℹ 34 more rows

    list_palettes("diverging")$name
#> [1] "bam"    "berlin" "broc"   "cork"   "lisbon" "roma"   "tofino" "vanimo"
#> [9] "vik"   
```
