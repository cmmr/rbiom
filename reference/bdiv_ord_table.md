# Calculate PCoA and other ordinations, including taxa biplots and statistics.

The biplot parameters (`taxa`, `unc`, `p.top`, and `p.adj`) only only
have an effect when `rank` is not `NULL`.

## Usage

``` r
bdiv_ord_table(
  biom,
  bdiv = "bray",
  ord = "PCoA",
  weighted = NULL,
  md = NULL,
  k = 2,
  stat.by = NULL,
  split.by = NULL,
  tree = NULL,
  test = "adonis2",
  seed = 0,
  permutations = 999,
  rank = NULL,
  taxa = 6,
  p.top = Inf,
  p.adj = "fdr",
  unc = "singly",
  alpha = 0.5,
  cpus = NULL,
  ...
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- bdiv:

  Beta diversity distance algorithm(s) to use. Options are:
  `c("aitchison", "bhattacharyya", "bray", "canberra", "chebyshev", "chord", "clark", "sorensen", "divergence", "euclidean", "generalized_unifrac", "gower", "hamming", "hellinger", "horn", "jaccard", "jensen", "jsd", "lorentzian", "manhattan", "matusita", "minkowski", "morisita", "motyka", "normalized_unifrac", "ochiai", "psym_chisq", "soergel", "squared_chisq", "squared_chord", "squared_euclidean", "topsoe", "unweighted_unifrac", "variance_adjusted_unifrac", "wave_hedges", "weighted_unifrac")`.
  For the UniFrac family, a phylogenetic tree must be present in `biom`
  or explicitly provided via `tree=`. Supports partial matching.
  Multiple values are allowed for functions which return a table or
  plot. Default: `"bray"`

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

- weighted:

  (Deprecated - weighting is now inherent in bdiv metric name.) Take
  relative abundances into account. When `weighted=FALSE`, only
  presence/absence is considered. Multiple values allowed. Default:
  `NULL`

- md:

  Dataset field(s) to include in the output data frame, or `'.all'` to
  include all metadata fields. Default: `'.all'`

- k:

  Number of ordination dimensions to return. Either `2L` or `3L`.
  Default: `2L`

- stat.by:

  The categorical or numeric metadata field over which statistics should
  be calculated. Required.

- split.by:

  Dataset field(s) that the data should be split by prior to any
  calculations. Must be categorical. Default: `NULL`

- tree:

  A `phylo` object representing the phylogenetic relationships of the
  taxa in `biom`. Only required when computing UniFrac distances.
  Default: `biom$tree`

- test:

  Permutational test for accessing significance. Options are:

  `"adonis2"` -

  :   Permutational MANOVA;
      [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html).

  `"mrpp"` -

  :   Multiple response permutation procedure;
      [`vegan::mrpp()`](https://vegandevs.github.io/vegan/reference/mrpp.html).

  `"none"` -

  :   Don't run any statistics.

  Abbreviations are allowed. Default: `"adonis2"`

- seed:

  Random seed for permutations. Must be a non-negative integer. Default:
  `0`

- permutations:

  Number of random permutations to use. Default: `999`

- rank:

  What rank(s) of taxa to compute biplot coordinates and statistics for,
  or `NULL` to disable. E.g. `"Phylum"`, `"Genus"`, `".otu"`, etc. An
  integer vector can also be given, where `1` is the highest rank, `2`
  is the second highest, `-1` is the lowest rank, `-2` is the second
  lowest, and `0` is the OTU "rank". Run `biom$ranks` to see all options
  for a given rbiom object. Default: `NULL`.

- taxa:

  Which taxa to display. An integer value will show the top n most
  abundant taxa. A value 0 \<= n \< 1 will show any taxa with that mean
  abundance or greater (e.g. `0.1` implies \>= 10%). A character vector
  of taxa names will show only those named taxa. Default: `6`.

- p.top:

  Only display taxa with the most significant differences in abundance.
  If `p.top` is \>= 1, then the `p.top` most significant taxa are
  displayed. If `p.top` is less than one, all taxa with an adjusted
  p-value \<= `p.top` are displayed. Recommended to be used in
  combination with the `taxa` parameter to set a lower bound on the mean
  abundance of considered taxa. Default: `Inf`

- p.adj:

  Method to use for multiple comparisons adjustment of p-values. Run
  `p.adjust.methods` for a list of available options. Default: `"fdr"`

- unc:

  How to handle unclassified, uncultured, and similarly ambiguous taxa
  names. Options are:

  `"singly"` -

  :   Replaces them with the OTU name.

  `"grouped"` -

  :   Replaces them with a higher rank's name.

  `"drop"` -

  :   Excludes them from the result.

  `"asis"` -

  :   To not check/modify any taxa names.

  Abbreviations are allowed. Default: `"singly"`

- alpha:

  The alpha term to use in Generalized UniFrac. How much weight to give
  to relative abundances; a value between 0 and 1, inclusive. Setting
  `alpha=1` is equivalent to Normalized UniFrac. Default: `0.5`

- cpus:

  The number of CPUs to use. Set to `NULL` to use all available, or to
  `1` to disable parallel processing. Default: `NULL`

- ...:

  Additional arguments to pass on to
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html),
  [`ape::pcoa()`](https://rdrr.io/pkg/ape/man/pcoa.html),
  [`vegan::metaMDS()`](https://vegandevs.github.io/vegan/reference/metaMDS.html),
  or [`tsne::tsne()`](https://rdrr.io/pkg/tsne/man/tsne.html).

## Value

A data.frame with columns `.sample`, `.bdiv`, `.ord`, `.x`, `.y`, and
(optionally) `.z`. Any columns given by `md`, `split.by`, and `stat.by`
are included as well. If `stat.by` is given, then `$stats` and
`$stats$code)` are set. If `rank` is given, then `$taxa_coords`,
`$taxa_stats`, and `$taxa_stats$code` are set.

## See also

Other beta_diversity:
[`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md),
[`bdiv_clusters()`](https://cmmr.github.io/rbiom/reference/bdiv_clusters.md),
[`bdiv_corrplot()`](https://cmmr.github.io/rbiom/reference/bdiv_corrplot.md),
[`bdiv_heatmap()`](https://cmmr.github.io/rbiom/reference/bdiv_heatmap.md),
[`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md),
[`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md),
[`bdiv_table()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md),
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md)

Other ordination:
[`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md),
[`distmat_ord_table()`](https://cmmr.github.io/rbiom/reference/distmat_ord_table.md)

## Examples

``` r
    library(rbiom)
    
    ord <- bdiv_ord_table(hmp50, "bray", "pcoa", stat.by="Body Site", rank="g")
    head(ord)
#> # A tibble: 6 × 6
#>   .bdiv .ord  .sample     .x       .y `Body Site`  
#>   <chr> <fct> <chr>    <dbl>    <dbl> <fct>        
#> 1 bray  PCoA  HMP01   -0.387 -0.00895 Buccal mucosa
#> 2 bray  PCoA  HMP02   -0.461 -0.0135  Buccal mucosa
#> 3 bray  PCoA  HMP03   -0.407 -0.0143  Saliva       
#> 4 bray  PCoA  HMP04   -0.378 -0.0123  Saliva       
#> 5 bray  PCoA  HMP05   -0.452 -0.0114  Buccal mucosa
#> 6 bray  PCoA  HMP06   -0.411 -0.0123  Saliva       
    
    ord$stats
#> # Test:     adonis2 ~ `Body Site`. 999 permutations.
#> # A tibble: 1 × 6
#>   .bdiv    .n .stat    .z .p.val .adj.p
#>   <chr> <int> <dbl> <dbl>  <dbl>  <dbl>
#> 1 bray     50  19.2  72.2  0.001  0.001
    
    ord$taxa_stats
#> # Test:     adonis2 ~ taxa. 999 permutations.
#> # A tibble: 6 × 8
#>   .bdiv .rank .taxa                .n .stat    .z .p.val .adj.p
#>   <chr> <fct> <fct>             <int> <dbl> <dbl>  <dbl>  <dbl>
#> 1 bray  Genus Lactobacillus        50  8.62  16.8  0.001  0.001
#> 2 bray  Genus Streptococcus        50 10.1   19.5  0.001  0.001
#> 3 bray  Genus Bacteroides          50  8.72  16.7  0.001  0.001
#> 4 bray  Genus Corynebacterium 1    50  8.53  16.0  0.001  0.001
#> 5 bray  Genus Haemophilus          50  6.56  12.3  0.001  0.001
#> 6 bray  Genus Propionibacterium    50  5.84  11.3  0.001  0.001
    
    
```
