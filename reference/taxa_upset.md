# Display taxa shared between groups as an UpSet plot.

Draws the result of
[`taxa_shared()`](https://cmmr.github.io/rbiom/reference/taxa_shared.md):
a bar for each combination of groups showing how many taxa belong to
exactly that combination, above a matrix of dots that says which groups
the combination is made of. Unlike a Venn diagram, it stays readable
with more than three groups.

## Usage

``` r
taxa_upset(
  biom,
  stat.by,
  rank = -1,
  detect = 0,
  prev = 0,
  unc = "singly",
  lineage = FALSE,
  sets = 15,
  totals = TRUE,
  ...
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- stat.by:

  Dataset field with the groups to compare. Must be categorical, with at
  least two groups. Samples with `NA` in this field are left out.
  Required.

- rank:

  What rank(s) of taxa to display. E.g. `"Phylum"`, `"Genus"`, `".otu"`,
  etc. An integer vector can also be given, where `1` is the highest
  rank, `2` is the second highest, `-1` is the lowest rank, `-2` is the
  second lowest, and `0` is the OTU "rank". Run `biom$ranks` to see all
  options for a given rbiom object. Default: `-1`.

- detect:

  Detection threshold. A taxon is present in a sample when its relative
  abundance there is greater than `detect`, a fraction from 0 to 1
  (`0.001` is 0.1%). The default counts any observation. Default: `0`

- prev:

  Prevalence threshold. A taxon belongs to a group when it is present in
  at least this fraction of the group's samples, and in at least one of
  them. So `0` needs a single sample, `0.9` needs 90% of the group's
  samples, and `1` needs all of them. Default: `0`

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

- lineage:

  Include all ranks in the name of the taxa. For instance, setting to
  `TRUE` will produce
  `Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales`.
  Otherwise the taxa name will simply be `Coriobacteriales`. You want to
  set this to TRUE when `unc = "asis"` and you have taxa names (such as
  *Incertae_Sedis*) that map to multiple higher level ranks. Default:
  `FALSE`

- sets:

  The maximum number of combinations (bars) to show. The largest are
  kept. Default: `15`

- totals:

  Show a bar for each group with the total number of taxa that belong to
  it. Default: `TRUE`

- ...:

  Parameters for the underlying ggplot2 functions, prefixed with the
  name of the part of the figure they apply to: `bar.` (the bars for
  each combination;
  [`ggplot2::geom_col()`](https://ggplot2.tidyverse.org/reference/geom_bar.html)),
  `label.` (the counts above them;
  [`ggplot2::geom_text()`](https://ggplot2.tidyverse.org/reference/geom_text.html)),
  `dot.` and `empty.` (the filled and unfilled dots;
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)),
  `line.` (the lines joining dots;
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html)),
  `total.` (the bars for each group;
  [`ggplot2::geom_col()`](https://ggplot2.tidyverse.org/reference/geom_bar.html)),
  and `theme_rbiom.`
  ([`theme_rbiom()`](https://cmmr.github.io/rbiom/reference/theme_rbiom.md),
  for all panels). For example, `bar.fill = "steelblue"`. An argument
  that no part accepts is an error.

## Value

A `ggplot2` plot (a patchwork of panels). The
[`taxa_shared()`](https://cmmr.github.io/rbiom/reference/taxa_shared.md)
table is available as `$data` and `$stats`, and the code that rebuilds
the figure as `$code`.

## Details

Whether a rare taxon is observed in a sample depends on how deeply the
sample was sequenced, so give this function rarefied counts (see
[`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md)).

## See also

Other taxa_abundance:
[`sample_sums()`](https://cmmr.github.io/rbiom/reference/sample_sums.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_clusters()`](https://cmmr.github.io/rbiom/reference/taxa_clusters.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_da()`](https://cmmr.github.io/rbiom/reference/taxa_da.md),
[`taxa_forest()`](https://cmmr.github.io/rbiom/reference/taxa_forest.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_shared()`](https://cmmr.github.io/rbiom/reference/taxa_shared.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md),
[`taxa_sums()`](https://cmmr.github.io/rbiom/reference/taxa_sums.md),
[`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md),
[`taxa_top()`](https://cmmr.github.io/rbiom/reference/taxa_top.md)

Other visualization:
[`adiv_boxplot()`](https://cmmr.github.io/rbiom/reference/adiv_boxplot.md),
[`adiv_corrplot()`](https://cmmr.github.io/rbiom/reference/adiv_corrplot.md),
[`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md),
[`bdiv_corrplot()`](https://cmmr.github.io/rbiom/reference/bdiv_corrplot.md),
[`bdiv_heatmap()`](https://cmmr.github.io/rbiom/reference/bdiv_heatmap.md),
[`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md),
[`plot_heatmap()`](https://cmmr.github.io/rbiom/reference/plot_heatmap.md),
[`rare_corrplot()`](https://cmmr.github.io/rbiom/reference/rare_corrplot.md),
[`rare_multiplot()`](https://cmmr.github.io/rbiom/reference/rare_multiplot.md),
[`rare_stacked()`](https://cmmr.github.io/rbiom/reference/rare_stacked.md),
[`stats_boxplot()`](https://cmmr.github.io/rbiom/reference/stats_boxplot.md),
[`stats_corrplot()`](https://cmmr.github.io/rbiom/reference/stats_corrplot.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_forest()`](https://cmmr.github.io/rbiom/reference/taxa_forest.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md),
[`theme_rbiom()`](https://cmmr.github.io/rbiom/reference/theme_rbiom.md)

## Examples

``` r
    library(rbiom)

    biom <- rarefy(hmp50)

    taxa_upset(biom, stat.by = "Body Site")


    # Core genera: present in at least half of a body site's samples.
    p <- taxa_upset(biom, stat.by = "Body Site", prev = 0.5, bar.fill = "steelblue")
    p


    # The taxa behind the bars.
    p$stats
#> # A tibble: 83 × 13
#>    .taxa    `Anterior nares` `Buccal mucosa` `Mid vagina` Saliva Stool .n_groups
#>    <fct>    <lgl>            <lgl>           <lgl>        <lgl>  <lgl>     <int>
#>  1 Alistip… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  2 Anaeros… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  3 Christe… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  4 Dorea    FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  5 Erysipe… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  6 Eubacte… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  7 Faecali… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  8 Fusicat… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  9 Lachnoc… FALSE            FALSE           FALSE        FALSE  TRUE          1
#> 10 Lachnos… FALSE            FALSE           FALSE        FALSE  TRUE          1
#> # ℹ 73 more rows
#> # ℹ 6 more variables: .set <fct>, `.prev.Anterior nares` <dbl>,
#> #   `.prev.Buccal mucosa` <dbl>, `.prev.Mid vagina` <dbl>, .prev.Saliva <dbl>,
#> #   .prev.Stool <dbl>
```
