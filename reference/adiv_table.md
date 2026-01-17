# Calculate the alpha diversity of each sample.

Calculate the alpha diversity of each sample.

## Usage

``` r
adiv_table(
  biom,
  adiv = "shannon",
  md = ".all",
  tree = NULL,
  transform = "none",
  ties = "random",
  seed = 0,
  cpus = NULL
)

adiv_matrix(
  biom,
  adiv = c("observed", "shannon", "simpson"),
  tree = NULL,
  transform = "none",
  ties = "random",
  seed = 0,
  cpus = NULL
)

adiv_vector(
  biom,
  adiv = "shannon",
  tree = NULL,
  transform = "none",
  ties = "random",
  seed = 0,
  cpus = NULL
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- adiv:

  Alpha diversity metric(s) to use. Options are:
  `c("ace", "berger", "brillouin", "chao1", "faith", "fisher", "simpson", "inv_simpson", "margalef", "mcintosh", "menhinick", "observed", "shannon", "squares")`.
  For `"faith"`, a phylogenetic tree must be present in `biom` or
  explicitly provided via `tree=`. Set `adiv=".all"` to use all metrics.
  Multiple/abbreviated values allowed. Default: `"shannon"`

- md:

  Dataset field(s) to include in the output data frame, or `'.all'` to
  include all metadata fields. Default: `'.all'`

- tree:

  A `phylo` object representing the phylogenetic relationships of the
  taxa in `biom`. Only required when computing UniFrac distances.
  Default: `biom$tree`

- transform:

  Transformation to apply. Options are:
  `c("none", "rank", "log", "log1p", "sqrt", "percent")`. `"rank"` is
  useful for correcting for non-normally distributions before applying
  regression statistics. Default: `"none"`

- ties:

  When `transform="rank"`, how to rank identical values. Options are:
  `c("average", "first", "last", "random", "max", "min")`. See
  [`rank()`](https://rdrr.io/r/base/rank.html) for details. Default:
  `"random"`

- seed:

  Random seed for permutations. Must be a non-negative integer. Default:
  `0`

- cpus:

  The number of CPUs to use. Set to `NULL` to use all available, or to
  `1` to disable parallel processing. Default: `NULL`

## Value

- `adiv_vector()` - :

  A named numeric vector.

- `adiv_matrix()` - :

  A matrix of samples x metric. The first column, 'depth', is never
  transformed.

- `adiv_table()` - :

  A tibble data.frame of alpha diversity values. Each combination of
  sample/`adiv` has its own row. Column names are **.sample**,
  **.depth**, **.adiv**, and **.diversity**, followed by any metadata
  fields requested by `md`.

## See also

[`sample_sums()`](https://cmmr.github.io/rbiom/reference/sample_sums.md)
for sample depths.

Other alpha_diversity:
[`adiv_boxplot()`](https://cmmr.github.io/rbiom/reference/adiv_boxplot.md),
[`adiv_corrplot()`](https://cmmr.github.io/rbiom/reference/adiv_corrplot.md),
[`adiv_stats()`](https://cmmr.github.io/rbiom/reference/adiv_stats.md)

## Examples

``` r
    library(rbiom)
    
    biom <- hmp50[1:5]
    
    adiv_table(biom)
#> # A tibble: 5 × 8
#>   .sample .depth .adiv   .diversity   Age   BMI `Body Site`   Sex   
#>   <chr>    <dbl> <fct>        <dbl> <dbl> <dbl> <fct>         <fct> 
#> 1 HMP01     1660 shannon       1.74    22    20 Buccal mucosa Female
#> 2 HMP02     1371 shannon       2.59    24    23 Buccal mucosa Male  
#> 3 HMP03     1353 shannon       2.95    28    26 Saliva        Male  
#> 4 HMP04     1895 shannon       3.26    25    23 Saliva        Male  
#> 5 HMP05     3939 shannon       1.46    27    24 Buccal mucosa Female
    
    biom <- rarefy(biom)
    adiv_table(biom, md = NULL)
#> # A tibble: 5 × 4
#>   .sample .depth .adiv   .diversity
#>   <fct>    <dbl> <fct>        <dbl>
#> 1 HMP01     1353 shannon       1.75
#> 2 HMP02     1353 shannon       2.59
#> 3 HMP03     1353 shannon       2.95
#> 4 HMP04     1353 shannon       3.24
#> 5 HMP05     1353 shannon       1.43
    
    adiv_vector(biom, 'faith')
    
    adiv_matrix(biom)
#>       depth observed  shannon   simpson
#> HMP01  1353       49 1.745200 0.5724740
#> HMP02  1353       75 2.587377 0.8125427
#> HMP03  1353       75 2.950982 0.8936622
#> HMP04  1353       73 3.235844 0.9315730
#> HMP05  1353       44 1.430414 0.5248395
#> attr(,"cmd")
#> [1] "adiv_matrix(biom, c(\"observed\", \"shannon\", \"simpson\"))"
```
