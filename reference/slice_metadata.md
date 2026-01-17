# Subset to a specific number of samples.

Subset to a specific number of samples.

## Usage

``` r
# S3 method for class 'rbiom'
slice(.data, ..., .by = NULL, .preserve = FALSE, clone = TRUE)

# S3 method for class 'rbiom'
slice_head(.data, n, prop, by = NULL, clone = TRUE, ...)

# S3 method for class 'rbiom'
slice_tail(.data, n, prop, by = NULL, clone = TRUE, ...)

# S3 method for class 'rbiom'
slice_min(
  .data,
  order_by,
  n,
  prop,
  by = NULL,
  with_ties = TRUE,
  na_rm = FALSE,
  clone = TRUE,
  ...
)

# S3 method for class 'rbiom'
slice_max(
  .data,
  order_by,
  n,
  prop,
  by = NULL,
  with_ties = TRUE,
  na_rm = FALSE,
  clone = TRUE,
  ...
)

# S3 method for class 'rbiom'
slice_sample(
  .data,
  n,
  prop,
  by = NULL,
  weight_by = NULL,
  replace = FALSE,
  clone = TRUE,
  ...
)
```

## Arguments

- .data:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), such
  as from
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- ...:

  For [`slice()`](https://dplyr.tidyverse.org/reference/slice.html),
  integer row indexes. For other `slice_*()` functions, not used. See
  [`dplyr::slice()`](https://dplyr.tidyverse.org/reference/slice.html).

- .by, by:

  **\[experimental\]**

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optionally, a selection of columns to group by for just this
  operation, functioning as an alternative to
  [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
  For details and examples, see
  [?dplyr_by](https://dplyr.tidyverse.org/reference/dplyr_by.html).

- .preserve:

  Relevant when the `.data` input is grouped. If `.preserve = FALSE`
  (the default), the grouping structure is recalculated based on the
  resulting data, otherwise the grouping is kept as is.

- clone:

  Create a copy of `biom` before modifying. If `FALSE`, `biom` is
  modified in place as a side-effect. See [speed
  ups](https://cmmr.github.io/rbiom/reference/speed.md) for use cases.
  Default: `TRUE`

- n, prop:

  Provide either `n`, the number of rows, or `prop`, the proportion of
  rows to select. If neither are supplied, `n = 1` will be used. If `n`
  is greater than the number of rows in the group (or `prop > 1`), the
  result will be silently truncated to the group size. `prop` will be
  rounded towards zero to generate an integer number of rows.

  A negative value of `n` or `prop` will be subtracted from the group
  size. For example, `n = -2` with a group of 5 rows will select 5 - 2 =
  3 rows; `prop = -0.25` with 8 rows will select 8 \* (1 - 0.25) = 6
  rows.

- order_by:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Variable or function of variables to order by. To order by multiple
  variables, wrap them in a data frame or tibble.

- with_ties:

  Should ties be kept together? The default, `TRUE`, may return more
  rows than you request. Use `FALSE` to ignore ties, and return the
  first `n` rows.

- na_rm:

  Should missing values in `order_by` be removed from the result? If
  `FALSE`, `NA` values are sorted to the end (like in
  [`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html)), so
  they will only be included if there are insufficient non-missing
  values to reach `n`/`prop`.

- weight_by:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Sampling weights. This must evaluate to a vector of non-negative
  numbers the same length as the input. Weights are automatically
  standardised to sum to 1.

- replace:

  Should sampling be performed with (`TRUE`) or without (`FALSE`, the
  default) replacement.

## Value

An [rbiom
object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md).

## See also

Other transformations:
[`matrix_ops`](https://cmmr.github.io/rbiom/reference/matrix_ops.md),
[`modify_metadata`](https://cmmr.github.io/rbiom/reference/modify_metadata.md),
[`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md),
[`subset()`](https://cmmr.github.io/rbiom/reference/subset.md),
[`with()`](https://cmmr.github.io/rbiom/reference/with.md)

## Examples

``` r
    library(rbiom)
    
    # The last 3 samples in the metadata table.
    biom <- slice_tail(hmp50, n = 3)
    biom$metadata
#> # A tibble: 3 × 5
#>   .sample   Age   BMI `Body Site`    Sex   
#> * <chr>   <dbl> <dbl> <fct>          <fct> 
#> 1 HMP48      24    29 Anterior nares Male  
#> 2 HMP49      21    22 Mid vagina     Female
#> 3 HMP50      31    20 Mid vagina     Female
    
    # The 3 oldest subjects sampled.
    biom <- slice_max(hmp50, Age, n = 3)
    biom$metadata
#> # A tibble: 3 × 5
#>   .sample   Age   BMI `Body Site`    Sex   
#> * <chr>   <dbl> <dbl> <fct>          <fct> 
#> 1 HMP34      40    27 Anterior nares Female
#> 2 HMP37      38    31 Mid vagina     Female
#> 3 HMP44      38    19 Mid vagina     Female
    
    # Pick 3 samples at random.
    biom <- slice_sample(hmp50, n = 3)
    biom$metadata
#> # A tibble: 3 × 5
#>   .sample   Age   BMI `Body Site`   Sex   
#> * <chr>   <dbl> <dbl> <fct>         <fct> 
#> 1 HMP43      24    30 Mid vagina    Female
#> 2 HMP07      26    22 Buccal mucosa Male  
#> 3 HMP29      36    25 Saliva        Female
```
