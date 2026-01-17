# Apply a function to each subset of an rbiom object.

`blply()` and `bdply()` let you divide your biom dataset into smaller
pieces, run a function on those smaller rbiom objects, and return the
results as a data.frame or list.

## Usage

``` r
bdply(biom, vars, FUN, ..., iters = list(), prefix = FALSE)

blply(biom, vars, FUN, ..., iters = list(), prefix = FALSE)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- vars:

  A character vector of metadata fields. Each unique combination of
  values in these columns will be used to create a subsetted rbiom
  object to pass to `FUN.` If `NULL`, `biom` will be passed to `FUN`
  unaltered. Unambiguous abbreviations of metadata fields are also
  accepted.

- FUN:

  The function to execute on each subset of `biom`. For `bdply()`, the
  returned value will be coerced to a data.frame. For `blply()`, any
  returned value is unmodified.

- ...:

  Additional arguments to pass on to `FUN`.

- iters:

  A named list of values to pass to `FUN`. Unlike `...`, these will be
  iterated over in all combinations. Default:
  [`list()`](https://rdrr.io/r/base/list.html)

- prefix:

  When `TRUE`, prefixes the names in in `iters` with a '.' in the final
  data.frame or 'split_labels' attribute. Default: `FALSE`

## Value

For `bdply()`, a tibble data.frame comprising the accumulated outputs of
`FUN`, along with the columns specified by `vars` and `iters`. For
`blply()`, a named list that has details about `vars` and `iters` in
`attr(,'split_labels')`.

## Details

You can also specify additional variables for your function to iterate
over in unique combinations.

Calls [`plyr::ddply()`](https://rdrr.io/pkg/plyr/man/ddply.html) or
[`plyr::dlply()`](https://rdrr.io/pkg/plyr/man/dlply.html) internally.

## See also

Other metadata:
[`glimpse.rbiom()`](https://cmmr.github.io/rbiom/reference/glimpse.rbiom.md)

Other biom:
[`biom_merge()`](https://cmmr.github.io/rbiom/reference/biom_merge.md)

## Examples

``` r
    library(rbiom)
    
    bdply(hmp50, "Sex", `$`, 'n_samples')
#> # A tibble: 2 × 2
#>   Sex       V1
#>   <fct>  <int>
#> 1 Female    30
#> 2 Male      20
    
    blply(hmp50, "Sex", `$`, 'n_samples') %>% unlist()
#> Female   Male 
#>     30     20 
    
    bdply(hmp50, c("Body Site", "Sex"), function (b) {
      adm <- adiv_matrix(b)[,c("shannon", "simpson")]
      apply(adm, 2L, mean)
    })
#> # A tibble: 9 × 4
#>   `Body Site`    Sex    shannon simpson
#>   <fct>          <fct>    <dbl>   <dbl>
#> 1 Anterior nares Female   1.43    0.681
#> 2 Anterior nares Male     1.51    0.665
#> 3 Buccal mucosa  Female   1.17    0.408
#> 4 Buccal mucosa  Male     1.71    0.602
#> 5 Mid vagina     Female   0.407   0.167
#> 6 Saliva         Female   2.93    0.893
#> 7 Saliva         Male     3.17    0.913
#> 8 Stool          Female   2.43    0.850
#> 9 Stool          Male     2.51    0.835
    
    iters <- list(d = c("bray", "euclid"))
    bdply(hmp50, "Sex", iters = iters, function (b, d) {
      r <- range(bdiv_distmat(biom = b, bdiv = d))
      round(data.frame(min = r[[1]], max = r[[2]]))
    })
#> # A tibble: 4 × 4
#>   Sex    d        min   max
#>   <fct>  <chr>  <dbl> <dbl>
#> 1 Female bray       0     1
#> 2 Female euclid     0     1
#> 3 Male   bray       0     1
#> 4 Male   euclid     0     1
```
