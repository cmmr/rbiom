# Run ordinations on a distance matrix.

Run ordinations on a distance matrix.

## Usage

``` r
distmat_ord_table(dm, ord = "PCoA", k = 2L, ...)
```

## Arguments

- dm:

  A `dist`-class distance matrix, as returned from
  [`bdiv_distmat()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md)
  or [`stats::dist()`](https://rdrr.io/r/stats/dist.html). Required.

- ord:

  Method for reducing dimensionality. Options are:

  `"PCoA"` -

  :   Principal coordinate analysis;
      [`ape::pcoa()`](https://rdrr.io/pkg/ape/man/pcoa.html).

  `"UMAP"` -

  :   Uniform manifold approximation and projection;
      [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).

  `"NMDS"` -

  :   Nonmetric multidimensional scaling;
      [`vegan::metaMDS()`](https://vegandevs.github.io/vegan/reference/metaMDS.html).

  `"tSNE"` -

  :   t-distributed stochastic neighbor embedding;
      [`tsne::tsne()`](https://rdrr.io/pkg/tsne/man/tsne.html).

  Multiple/abbreviated values allowed. Default: `"PCoA"`

- k:

  Number of ordination dimensions to return. Either `2L` or `3L`.
  Default: `2L`

- ...:

  Additional arguments for `ord`.

## Value

A data.frame with columns `.sample`, `.ord`, `.x`, `.y`, and
(optionally) `.z`.

## See also

Other ordination:
[`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md),
[`bdiv_ord_table()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_table.md)

## Examples

``` r
    library(rbiom) 
    
    dm  <- bdiv_distmat(hmp50, "bray")
    ord <- distmat_ord_table(dm, "PCoA")
    head(ord)
#> # Ordination: pcoa(D = dm)
#> # A tibble:   6 × 4
#>   .ord  .sample     .x       .y
#>   <fct> <chr>    <dbl>    <dbl>
#> 1 PCoA  HMP01   -0.387 -0.00895
#> 2 PCoA  HMP02   -0.461 -0.0135 
#> 3 PCoA  HMP03   -0.407 -0.0143 
#> 4 PCoA  HMP04   -0.378 -0.0123 
#> 5 PCoA  HMP05   -0.452 -0.0114 
#> 6 PCoA  HMP06   -0.411 -0.0123 
    
```
