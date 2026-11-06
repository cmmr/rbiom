# Run ordinations on a distance matrix.

Run ordinations on a distance matrix.

## Usage

``` r
distmat_ord_table(dm, ord = "PCoA", k = 2L, seed = 0, ...)
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
      [`stats::cmdscale()`](https://rdrr.io/r/stats/cmdscale.html).

  `"UMAP"` -

  :   Uniform manifold approximation and projection;
      [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).

  `"NMDS"` -

  :   Nonmetric multidimensional scaling;
      [`vegan::metaMDS()`](https://vegandevs.github.io/vegan/reference/metaMDS.html).

  `"tSNE"` -

  :   t-distributed stochastic neighbor embedding;
      [`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html).

  Multiple/abbreviated values allowed. Default: `"PCoA"`

- k:

  Number of ordination dimensions to return. Either `2L` or `3L`.
  Default: `2L`

- seed:

  Random seed for the stochastic ordinations (NMDS, tSNE, UMAP), so that
  they are reproducible. PCoA is deterministic. Must be a non-negative
  integer. Default: `0`

- ...:

  Additional arguments for the `ord` function. For tSNE, `perplexity`
  defaults to 30, or to `(n - 1) / 3` for fewer than 91 samples.

## Value

A tibble data.frame with columns `.ord`, `.sample`, `.x`, `.y`, and
(optionally) `.z`. The attribute `'axes'` is a tibble with one row per
axis (`.ord`, `.axis`, `.label`, `.var`): for PCoA, `.var` is the
percent of the variance (the positive eigenvalues) explained by that
axis and `.label` reads e.g. `"PCoA 1 (24.3%)"`; for the other
ordinations `.var` is `NA`.

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
#> # Ordination: stats::cmdscale(dm, k = 2, eig = TRUE)
#> # A tibble:   6 × 4
#>   .ord  .sample     .x       .y
#>   <fct> <chr>    <dbl>    <dbl>
#> 1 PCoA  HMP01   -0.358  0.00236
#> 2 PCoA  HMP02   -0.420  0.0104 
#> 3 PCoA  HMP03   -0.382  0.0135 
#> 4 PCoA  HMP04   -0.388  0.00485
#> 5 PCoA  HMP05   -0.415 -0.0115 
#> 6 PCoA  HMP06   -0.414 -0.00703

    attr(ord, 'axes')
#> # A tibble: 2 × 4
#>   .ord  .axis .label          .var
#>   <fct> <chr> <chr>          <dbl>
#> 1 PCoA  .x    PCoA 1 (21.6%)  21.6
#> 2 PCoA  .y    PCoA 2 (14.7%)  14.7
```
