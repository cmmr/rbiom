# Taxa shared between groups, and the core microbiome.

Lists which taxa are found in which groups of samples. With a prevalence
cut-off (`prev`) the taxa listed for a group are its core microbiome.
[`taxa_upset()`](https://cmmr.github.io/rbiom/reference/taxa_upset.md)
plots the result.

## Usage

``` r
taxa_shared(
  biom,
  stat.by = NULL,
  rank = -1,
  detect = 0,
  prev = 0,
  unc = "singly",
  lineage = FALSE
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- stat.by:

  Dataset field with the groups to compare. Must be categorical. Samples
  with `NA` in this field are left out. `NULL` puts every sample in a
  single group named `"All"`, which gives a plain core microbiome list.
  Default: `NULL`

- rank:

  Which taxa rank to use. E.g. `"Phylum"`, `"Genus"`, `".otu"`, etc. An
  integer can also be given, where `1` is the highest rank, `2` is the
  second highest, `-1` is the lowest rank, `-2` is the second lowest,
  and `0` is the OTU "rank". Run `biom$ranks` to see all options for a
  given rbiom object. Default: `-1`.

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

## Value

A tibble data frame with one row per taxon that belongs to at least one
group, and columns:

- `.taxa`:

  The taxon's name.

- One column per group:

  Named after the levels of `stat.by`, in level order: `TRUE` when the
  taxon belongs to that group.

- `.n_groups`:

  How many groups the taxon belongs to.

- `.set`:

  The combination of groups the taxon belongs to, e.g.
  `"Saliva & Stool"`.

- `.prev.<group>`:

  One column per group: the fraction of the group's samples in which the
  taxon is present.

Rows are sorted by `.set` (the largest set first), then by name.

The number of taxa in each set is attached as `attr(, "sets")`: a tibble
with `.set`, `.n` (taxa belonging to exactly that combination of
groups), `.n_groups`, and one logical column per group, sorted from the
largest set to the smallest.

The result carries `$code`: the R code that regenerates it, including
the steps that produced `biom`.

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
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md),
[`taxa_sums()`](https://cmmr.github.io/rbiom/reference/taxa_sums.md),
[`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md),
[`taxa_top()`](https://cmmr.github.io/rbiom/reference/taxa_top.md),
[`taxa_upset()`](https://cmmr.github.io/rbiom/reference/taxa_upset.md)

## Examples

``` r
    library(rbiom)

    biom <- rarefy(hmp50)

    # Genera found at each body site.
    shared <- taxa_shared(biom, stat.by = "Body Site")
    shared
#> # A tibble: 219 × 13
#>    .taxa    `Anterior nares` `Buccal mucosa` `Mid vagina` Saliva Stool .n_groups
#>    <fct>    <lgl>            <lgl>           <lgl>        <lgl>  <lgl>     <int>
#>  1 Acidami… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  2 Akkerma… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  3 Bilophi… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  4 Blautia  FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  5 Butyric… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  6 Butyriv… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  7 Candida… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  8 Clostri… FALSE            FALSE           FALSE        FALSE  TRUE          1
#>  9 Collins… FALSE            FALSE           FALSE        FALSE  TRUE          1
#> 10 Coproba… FALSE            FALSE           FALSE        FALSE  TRUE          1
#> # ℹ 209 more rows
#> # ℹ 6 more variables: .set <fct>, `.prev.Anterior nares` <dbl>,
#> #   `.prev.Buccal mucosa` <dbl>, `.prev.Mid vagina` <dbl>, .prev.Saliva <dbl>,
#> #   .prev.Stool <dbl>

    attr(shared, "sets")
#> # A tibble: 24 × 8
#>    .set        .n .n_groups `Anterior nares` `Buccal mucosa` `Mid vagina` Saliva
#>    <fct>    <int>     <int> <lgl>            <lgl>           <lgl>        <lgl> 
#>  1 Stool       63         1 FALSE            FALSE           FALSE        FALSE 
#>  2 Anterio…    25         1 TRUE             FALSE           FALSE        FALSE 
#>  3 Buccal …    19         2 FALSE            TRUE            FALSE        TRUE  
#>  4 Mid vag…    16         1 FALSE            FALSE           TRUE         FALSE 
#>  5 Saliva      15         1 FALSE            FALSE           FALSE        TRUE  
#>  6 Anterio…    14         3 TRUE             TRUE            FALSE        TRUE  
#>  7 Anterio…    14         2 TRUE             FALSE           FALSE        FALSE 
#>  8 Anterio…     7         5 TRUE             TRUE            TRUE         TRUE  
#>  9 Anterio…     7         4 TRUE             TRUE            TRUE         TRUE  
#> 10 Anterio…     6         2 TRUE             FALSE           TRUE         FALSE 
#> # ℹ 14 more rows
#> # ℹ 1 more variable: Stool <lgl>

    # Core phyla: above 0.1% in at least 90% of a body site's samples.
    taxa_shared(biom, "Body Site", rank = "Phylum", detect = 0.001, prev = 0.9)
#> # A tibble: 5 × 13
#>   .taxa     `Anterior nares` `Buccal mucosa` `Mid vagina` Saliva Stool .n_groups
#>   <fct>     <lgl>            <lgl>           <lgl>        <lgl>  <lgl>     <int>
#> 1 Firmicut… TRUE             TRUE            TRUE         TRUE   TRUE          5
#> 2 Bacteroi… FALSE            TRUE            FALSE        TRUE   TRUE          3
#> 3 Actinoba… TRUE             FALSE           FALSE        TRUE   FALSE         2
#> 4 Proteoba… FALSE            TRUE            FALSE        TRUE   FALSE         2
#> 5 Fusobact… FALSE            FALSE           FALSE        TRUE   FALSE         1
#> # ℹ 6 more variables: .set <fct>, `.prev.Anterior nares` <dbl>,
#> #   `.prev.Buccal mucosa` <dbl>, `.prev.Mid vagina` <dbl>, .prev.Saliva <dbl>,
#> #   .prev.Stool <dbl>

    # Core genera across all samples.
    taxa_shared(biom, prev = 0.5)
#> # A tibble: 6 × 5
#>   .taxa          All   .n_groups .set  .prev.All
#>   <fct>          <lgl>     <int> <fct>     <dbl>
#> 1 Alloprevotella TRUE          1 All       0.510
#> 2 Gemella        TRUE          1 All       0.551
#> 3 Granulicatella TRUE          1 All       0.510
#> 4 Haemophilus    TRUE          1 All       0.551
#> 5 Prevotella     TRUE          1 All       0.592
#> 6 Streptococcus  TRUE          1 All       0.714
```
