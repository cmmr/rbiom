# Ordinate samples and taxa on a 2D plane based on beta diversity distances.

Ordinate samples and taxa on a 2D plane based on beta diversity
distances.

## Usage

``` r
bdiv_ord_plot(
  biom,
  bdiv = "bray",
  ord = "PCoA",
  layers = "petm",
  stat.by = NULL,
  facet.by = NULL,
  colors = TRUE,
  shapes = TRUE,
  tree = NULL,
  test = "adonis2",
  seed = 0,
  permutations = 999,
  rank = -1,
  taxa = 4,
  p.top = Inf,
  p.adj = "fdr",
  unc = "singly",
  caption = TRUE,
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

- layers:

  One or more of
  `c("point", "spider", "ellipse", "name", "mean", "taxon", "arrow")`.
  The first four are sample-centric; the last three are taxa-centric.
  Single letter abbreviations are also accepted. For instance,
  `c("point", "ellipse")` is equivalent to `c("p", "e")` and `"pe"`.
  Default: `"pe"`

- stat.by:

  The categorical or numeric metadata field over which statistics should
  be calculated. Required.

- facet.by:

  Dataset field(s) to use for faceting. Must be categorical. Default:
  `NULL`

- colors:

  How to color the groups. Options are:

  `TRUE` -

  :   Automatically select colorblind-friendly colors.

  `FALSE` or `NULL` -

  :   Don't use colors.

  a palette name -

  :   Auto-select colors from this set. E.g. `"okabe"`

  character vector -

  :   Custom colors to use. E.g. `c("red", "#00FF00")`

  named character vector -

  :   Explicit mapping. E.g. `c(Male = "blue", Female = "red")`

  See "Aesthetics" section below for additional information. Default:
  `TRUE`

- shapes:

  Shapes for each group. Options are similar to `colors`'s: `TRUE`,
  `FALSE`, `NULL`, shape names (typically integers 0 - 17), or a named
  vector mapping groups to specific shape names. See "Aesthetics"
  section below for additional information. Default: `TRUE`

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

  What rank(s) of taxa to display. E.g. `"Phylum"`, `"Genus"`, `".otu"`,
  etc. An integer vector can also be given, where `1` is the highest
  rank, `2` is the second highest, `-1` is the lowest rank, `-2` is the
  second lowest, and `0` is the OTU "rank". Run `biom$ranks` to see all
  options for a given rbiom object. Default: `-1`.

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

- caption:

  Add methodology caption beneath the plot. Default: `TRUE`

- alpha:

  The alpha term to use in Generalized UniFrac. How much weight to give
  to relative abundances; a value between 0 and 1, inclusive. Setting
  `alpha=1` is equivalent to Normalized UniFrac. Default: `0.5`

- cpus:

  The number of CPUs to use. Set to `NULL` to use all available, or to
  `1` to disable parallel processing. Default: `NULL`

- ...:

  Parameters for layer geoms (e.g.
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)).
  Prefixing parameter names with a layer name ensures that a particular
  parameter is passed to, and only to, that layer. For instance,
  `point.size = 2` or `p.size = 2` ensures only the points have their
  size set to `2`. Points can also be controlled with the `pt.` prefix.

## Value

A `ggplot2` plot. The computed sample coordinates and ggplot command are
available as `$data` and `$code` respectively. If `stat.by` is given,
then `$stats` and `$stats$code` are set. If `rank` is given, then
`$data$taxa_coords`, `$taxa_stats`, and `$taxa_stats$code` are set.

## See also

Other beta_diversity:
[`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md),
[`bdiv_clusters()`](https://cmmr.github.io/rbiom/reference/bdiv_clusters.md),
[`bdiv_corrplot()`](https://cmmr.github.io/rbiom/reference/bdiv_corrplot.md),
[`bdiv_heatmap()`](https://cmmr.github.io/rbiom/reference/bdiv_heatmap.md),
[`bdiv_ord_table()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_table.md),
[`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md),
[`bdiv_table()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md),
[`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md)

Other ordination:
[`bdiv_ord_table()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_table.md),
[`distmat_ord_table()`](https://cmmr.github.io/rbiom/reference/distmat_ord_table.md)

Other visualization:
[`adiv_boxplot()`](https://cmmr.github.io/rbiom/reference/adiv_boxplot.md),
[`adiv_corrplot()`](https://cmmr.github.io/rbiom/reference/adiv_corrplot.md),
[`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md),
[`bdiv_corrplot()`](https://cmmr.github.io/rbiom/reference/bdiv_corrplot.md),
[`bdiv_heatmap()`](https://cmmr.github.io/rbiom/reference/bdiv_heatmap.md),
[`plot_heatmap()`](https://cmmr.github.io/rbiom/reference/plot_heatmap.md),
[`rare_corrplot()`](https://cmmr.github.io/rbiom/reference/rare_corrplot.md),
[`rare_multiplot()`](https://cmmr.github.io/rbiom/reference/rare_multiplot.md),
[`rare_stacked()`](https://cmmr.github.io/rbiom/reference/rare_stacked.md),
[`stats_boxplot()`](https://cmmr.github.io/rbiom/reference/stats_boxplot.md),
[`stats_corrplot()`](https://cmmr.github.io/rbiom/reference/stats_corrplot.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md)

## Examples

``` r
    library(rbiom)
    
    biom <- rarefy(hmp50)
    
    bdiv_ord_plot(biom, layers="pemt", stat.by="Body Site", rank="g")
#> Warning: Probable convergence failure

    
```
