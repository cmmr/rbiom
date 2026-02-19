# Taxa abundances per sample.

- `taxa_matrix()` - :

  Accepts a single `rank` and returns a matrix.

- `taxa_table()` - :

  Can accept more than one `rank` and returns a tibble data.frame.

## Usage

``` r
taxa_table(
  biom,
  rank = -1,
  taxa = 6,
  lineage = FALSE,
  md = ".all",
  unc = "singly",
  other = FALSE,
  transform = "none",
  ties = "random",
  seed = 0
)

taxa_matrix(
  biom,
  rank = -1,
  taxa = NULL,
  lineage = FALSE,
  sparse = FALSE,
  unc = "singly",
  other = FALSE,
  transform = "none",
  ties = "random",
  seed = 0
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

- taxa:

  Which taxa to display. An integer value will show the top n most
  abundant taxa. A value 0 \<= n \< 1 will show any taxa with that mean
  abundance or greater (e.g. `0.1` implies \>= 10%). A character vector
  of taxa names will show only those named taxa. Default: `6`.

- lineage:

  Include all ranks in the name of the taxa. For instance, setting to
  `TRUE` will produce
  `Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales`.
  Otherwise the taxa name will simply be `Coriobacteriales`. You want to
  set this to TRUE when `unc = "asis"` and you have taxa names (such as
  *Incertae_Sedis*) that map to multiple higher level ranks. Default:
  `FALSE`

- md:

  Dataset field(s) to include in the output data frame, or `'.all'` to
  include all metadata fields. Default: `'.all'`

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

- other:

  Sum all non-itemized taxa into an "Other" taxa. When `FALSE`, only
  returns taxa matched by the `taxa` argument. Specifying `TRUE` adds
  "Other" to the returned set. A string can also be given to imply
  `TRUE`, but with that value as the name to use instead of "Other".
  Default: `FALSE`

- transform:

  Transformation to apply to calculated values. Options are:
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

- sparse:

  If `TRUE`, returns a sparse matrix from the `Matrix` package,
  otherwise returns a normal R matrix object. Default: `FALSE`

## Value

- `taxa_matrix()` - :

  A numeric matrix with taxa as rows, and samples as columns.

- `taxa_table()` - :

  A tibble data frame with column names .sample, .taxa, .abundance, and
  any requested by `md`.

## See also

Other taxa_abundance:
[`sample_sums()`](https://cmmr.github.io/rbiom/reference/sample_sums.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_clusters()`](https://cmmr.github.io/rbiom/reference/taxa_clusters.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md),
[`taxa_sums()`](https://cmmr.github.io/rbiom/reference/taxa_sums.md)

## Examples

``` r
    library(rbiom)
    
    hmp50$ranks
#> [1] ".otu"    "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"  
    
    taxa_matrix(hmp50, 'Phylum')[1:4,1:6]
#>                     HMP01 HMP02 HMP03 HMP04 HMP05 HMP06
#> Actinobacteria         18    60   126   120    30    71
#> Bacteroidetes         276   221   313   218   144   880
#> Cyanobacteria           0     0     0     0     0     0
#> Deinococcus Thermus     0     0     0     0     0     0
    
    taxa_table(hmp50, 'Phylum')
#> # A tibble: 300 × 8
#>    .rank  .sample .taxa          .abundance   Age   BMI `Body Site`   Sex   
#>    <fct>  <chr>   <fct>               <dbl> <dbl> <dbl> <fct>         <fct> 
#>  1 Phylum HMP01   Firmicutes           1208    22    20 Buccal mucosa Female
#>  2 Phylum HMP01   Bacteroidetes         276    22    20 Buccal mucosa Female
#>  3 Phylum HMP01   Actinobacteria         18    22    20 Buccal mucosa Female
#>  4 Phylum HMP01   Proteobacteria        103    22    20 Buccal mucosa Female
#>  5 Phylum HMP01   Fusobacteria           41    22    20 Buccal mucosa Female
#>  6 Phylum HMP01   Tenericutes             0    22    20 Buccal mucosa Female
#>  7 Phylum HMP02   Firmicutes            931    24    23 Buccal mucosa Male  
#>  8 Phylum HMP02   Bacteroidetes         221    24    23 Buccal mucosa Male  
#>  9 Phylum HMP02   Actinobacteria         60    24    23 Buccal mucosa Male  
#> 10 Phylum HMP02   Proteobacteria        112    24    23 Buccal mucosa Male  
#> # ℹ 290 more rows
```
