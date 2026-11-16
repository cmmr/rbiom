# The most abundant taxa overall, per group, or per sample.

[`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)
with `taxa = 5` keeps the five taxa that are most abundant across all
samples. `taxa_top()` ranks the taxa separately within each group of
samples (or each sample), so the taxa listed can differ from one group
to the next.

## Usage

``` r
taxa_top(
  biom,
  rank = -1,
  n = 5,
  by = NULL,
  unc = "singly",
  lineage = FALSE,
  transform = "percent"
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- rank:

  What rank(s) of taxa to display. E.g. `"Phylum"`, `"Genus"`, `".otu"`,
  etc. An integer vector can also be given, where `1` is the highest
  rank, `2` is the second highest, `-1` is the lowest rank, `-2` is the
  second lowest, and `0` is the OTU "rank". Run `biom$ranks` to see all
  options for a given rbiom object. Default: `-1`.

- n:

  How many taxa to list for each group. Default: `5`

- by:

  Dataset field(s) that define the groups of samples to rank taxa
  within. Must be categorical. Several fields give one group per
  combination of their values. Use `".sample"` to rank taxa within each
  sample. `NULL` ranks taxa across all samples together. Samples with
  `NA` in a `by` field are left out. Default: `NULL`

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

- transform:

  Units for the abundances in the result. Options are:
  `c("none", "rank", "log", "log1p", "sqrt", "percent")`. `"percent"`
  gives relative abundances as fractions (`0.25` is 25%), and `"none"`
  gives counts. See
  [`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)
  for the others. The choice does not change which taxa are listed.
  Default: `"percent"`

## Value

A tibble data frame with one row per listed taxon and columns:

- `.rank`:

  The taxonomic rank.

- `by` field(s):

  The group. Not present when `by = NULL`.

- `.top`:

  The taxon's position in its group: `1` for the most abundant, up to
  `n`.

- `.taxa`:

  The taxon's name.

- `.mean`:

  Mean abundance over the group's samples, in `transform` units.

- `.prevalence`:

  Fraction of the group's samples in which the taxon was observed.

When `by` includes `".sample"`, each group is one sample, and the table
has `.abundance` (in `transform` units) in place of `.mean` and
`.prevalence`.

The result carries `$code`: the R code that regenerates it, including
the steps that produced `biom`.

## Details

Taxa are ranked by their mean relative abundance in the group's samples,
as
[`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)
ranks them for its `taxa` argument. Ties are broken alphabetically by
taxon name, so the result never depends on a random seed. Taxa that are
absent from a group are not listed, so a group can have fewer than `n`
rows.

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
[`taxa_upset()`](https://cmmr.github.io/rbiom/reference/taxa_upset.md)

## Examples

``` r
    library(rbiom)

    biom <- rarefy(hmp50)

    # The five most abundant genera overall.
    taxa_top(biom)
#> # A tibble: 5 × 5
#>   .rank  .top .taxa              .mean .prevalence
#>   <fct> <int> <fct>              <dbl>       <dbl>
#> 1 Genus     1 Lactobacillus     0.181        0.286
#> 2 Genus     2 Streptococcus     0.168        0.714
#> 3 Genus     3 Bacteroides       0.0956       0.449
#> 4 Genus     4 Corynebacterium 1 0.0836       0.245
#> 5 Genus     5 Haemophilus       0.0474       0.551

    # The three most abundant phyla at each body site.
    taxa_top(biom, rank = "Phylum", n = 3, by = "Body Site")
#> # A tibble: 15 × 6
#>    .rank  `Body Site`     .top .taxa            .mean .prevalence
#>    <fct>  <fct>          <int> <fct>            <dbl>       <dbl>
#>  1 Phylum Anterior nares     1 Actinobacteria 0.657         1    
#>  2 Phylum Anterior nares     2 Firmicutes     0.270         1    
#>  3 Phylum Anterior nares     3 Proteobacteria 0.0540        0.8  
#>  4 Phylum Buccal mucosa      1 Firmicutes     0.788         1    
#>  5 Phylum Buccal mucosa      2 Proteobacteria 0.125         1    
#>  6 Phylum Buccal mucosa      3 Bacteroidetes  0.0554        1    
#>  7 Phylum Mid vagina         1 Firmicutes     0.974         1    
#>  8 Phylum Mid vagina         2 Bacteroidetes  0.0107        0.8  
#>  9 Phylum Mid vagina         3 Actinobacteria 0.00668       0.8  
#> 10 Phylum Saliva             1 Firmicutes     0.446         1    
#> 11 Phylum Saliva             2 Proteobacteria 0.261         1    
#> 12 Phylum Saliva             3 Bacteroidetes  0.181         1    
#> 13 Phylum Stool              1 Bacteroidetes  0.669         1    
#> 14 Phylum Stool              2 Firmicutes     0.287         1    
#> 15 Phylum Stool              3 Proteobacteria 0.0299        0.889

    # The dominant genus in each sample.
    taxa_top(biom, n = 1, by = ".sample")
#> # A tibble: 49 × 5
#>    .rank .sample  .top .taxa             .abundance
#>    <fct> <fct>   <int> <fct>                  <dbl>
#>  1 Genus HMP01       1 Streptococcus          0.657
#>  2 Genus HMP02       1 Streptococcus          0.461
#>  3 Genus HMP03       1 Streptococcus          0.249
#>  4 Genus HMP04       1 Fusobacterium          0.133
#>  5 Genus HMP05       1 Streptococcus          0.683
#>  6 Genus HMP06       1 Haemophilus            0.205
#>  7 Genus HMP07       1 Streptococcus          0.761
#>  8 Genus HMP08       1 Haemophilus            0.347
#>  9 Genus HMP09       1 Haemophilus            0.142
#> 10 Genus HMP10       1 Corynebacterium 1      0.440
#> # ℹ 39 more rows
```
