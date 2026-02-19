# Get summary taxa abundances.

Get summary taxa abundances.

## Usage

``` r
taxa_sums(
  biom,
  rank = -1,
  sort = NULL,
  lineage = FALSE,
  unc = "singly",
  transform = "none"
)

taxa_means(
  biom,
  rank = -1,
  sort = NULL,
  lineage = FALSE,
  unc = "singly",
  transform = "none"
)

taxa_apply(
  biom,
  FUN,
  rank = -1,
  sort = NULL,
  lineage = FALSE,
  unc = "singly",
  transform = "none",
  ...
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

- sort:

  Sort the result. Options: `NULL`, `"asc"`, or `"desc"`, where `NULL`
  will not sort the result. `"asc"` will sort in ascending order
  (smallest to largest), and `"desc"` in descending order (largest to
  smallest). Ignored when the result is not a simple numeric vector.
  Default: `NULL`

- lineage:

  Include all ranks in the name of the taxa. For instance, setting to
  `TRUE` will produce
  `Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales`.
  Otherwise the taxa name will simply be `Coriobacteriales`. You want to
  set this to TRUE when `unc = "asis"` and you have taxa names (such as
  *Incertae_Sedis*) that map to multiple higher level ranks. Default:
  `FALSE`

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

- transform:

  Transformation to apply to calculated values. Options are:
  `c("none", "rank", "log", "log1p", "sqrt", "percent")`. `"rank"` is
  useful for correcting for non-normally distributions before applying
  regression statistics. Default: `"none"`

- FUN:

  The function to apply to each row of the
  [`taxa_matrix()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md).

- ...:

  Optional arguments to `FUN`.

## Value

For `taxa_sums` and `taxa_means`, a named numeric vector. For
`taxa_apply`, a named vector or list with the results of `FUN`. The
names are the taxa IDs.

## See also

Other taxa_abundance:
[`sample_sums()`](https://cmmr.github.io/rbiom/reference/sample_sums.md),
[`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md),
[`taxa_clusters()`](https://cmmr.github.io/rbiom/reference/taxa_clusters.md),
[`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md),
[`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md),
[`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md),
[`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md),
[`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)

## Examples

``` r
    library(rbiom) 
    
    taxa_sums(hmp50) %>% head(4)
#>     Abiotrophia Acidaminococcus   Acinetobacter  Actinobacillus 
#>             135              62              44             242 
    
    taxa_means(hmp50, 'Family') %>% head(5)
#>  Acidaminococcaceae    Actinomycetaceae       Aerococcaceae      Alcaligenaceae 
#>               14.58               16.92                2.86               11.28 
#> Alicyclobacillaceae 
#>                0.12 
    
    taxa_apply(hmp50, max) %>% head(5)
#>     Abiotrophia Acidaminococcus   Acinetobacter  Actinobacillus     Actinomyces 
#>              78              59              42             100             105 
    
    taxa_apply(hmp50, fivenum) %>% head(5)
#>      Abiotrophia Acidaminococcus Acinetobacter Actinobacillus Actinomyces
#> [1,]           0               0             0              0           0
#> [2,]           0               0             0              0           0
#> [3,]           0               0             0              0           1
#> [4,]           1               0             0              2          21
#> [5,]          78              59            42            100         105
#>      Actinotignum Aerococcus Akkermansia Alistipes Alloprevotella Alloscardovia
#> [1,]            0          0           0         0            0.0             0
#> [2,]            0          0           0         0            0.0             0
#> [3,]            0          0           0         0            0.5             0
#> [4,]            0          0           0         3           12.0             0
#> [5,]            1          2         107       583          135.0             1
#>      Anaerococcus Anaeroglobus Anaerostipes Anaerotruncus Atopobium Bacteroides
#> [1,]            0            0            0             0         0           0
#> [2,]            0            0            0             0         0           0
#> [3,]            0            0            0             0         0           2
#> [4,]            0            0            0             0         2           7
#> [5,]          488            1           14            28       107        5786
#>      Barnesiella Bergeyella Bifidobacterium Bilophila Blautia Bosea
#> [1,]           0          0               0         0       0     0
#> [2,]           0          0               0         0       0     0
#> [3,]           0          0               0         0       0     0
#> [4,]           0         10               0         0       0     0
#> [5,]          59        114               3         2      67     2
#>      Bradyrhizobium Brevundimonas Brochothrix Burkholderia Butyricimonas
#> [1,]              0             0           0            0             0
#> [2,]              0             0           0            0             0
#> [3,]              0             0           0            0             0
#> [4,]              0             0           0            0             0
#> [5,]              4            16          11           15            77
#>      Butyrivibrio Butyrivibrio 2 Campylobacter Candidatus Saccharimonas
#> [1,]            0              0             0                        0
#> [2,]            0              0             0                        0
#> [3,]            0              0             0                        0
#> [4,]            0              0             2                        0
#> [5,]          382              5            69                        4
#>      Candidatus Stoquefichus Capnocytophaga Cardiobacterium Catenibacterium
#> [1,]                       0              0               0               0
#> [2,]                       0              0               0               0
#> [3,]                       0              0               0               0
#> [4,]                       0             15               0               0
#> [5,]                       7            235               9               2
#>      Catonella Christensenellaceae R_7 Clostridiales XIII_AD3011
#> [1,]         0                       0                         0
#> [2,]         0                       0                         0
#> [3,]         0                       0                         0
#> [4,]         1                       0                         0
#> [5,]        27                      62                         2
#>      Clostridiales XIII_UCG_001 Clostridium sensu_stricto_1 Collinsella
#> [1,]                          0                           0           0
#> [2,]                          0                           0           0
#> [3,]                          0                           0           0
#> [4,]                          0                           0           0
#> [5,]                          4                           1           8
#>      Comamonas Conchiformibius Coprobacter Coprococcus 1 Coprococcus 2
#> [1,]         0               0           0             0             0
#> [2,]         0               0           0             0             0
#> [3,]         0               0           0             0             0
#> [4,]         0               0           0             0             0
#> [5,]         4               2           1             6            10
#>      Corynebacterium Corynebacterium 1 Defluviitaleaceae UCG_011 Deinococcus
#> [1,]               0                 0                         0           0
#> [2,]               0                 0                         0           0
#> [3,]               0                 0                         0           0
#> [4,]               1                 1                         0           0
#> [5,]              32              4644                        12           3
#>      Dermabacter Desulfovibrio Dialister Dolosigranulum Dorea Eikenella
#> [1,]           0             0         0              0     0         0
#> [2,]           0             0         0              0     0         0
#> [3,]           0             0         0              0     0         0
#> [4,]           0             0         2              0     1         0
#> [5,]           1             3        32            707    59        22
#>      Elizabethkingia Enterococcus Erysipelatoclostridium
#> [1,]               0            0                      0
#> [2,]               0            0                      0
#> [3,]               0            0                      0
#> [4,]               0            0                      0
#> [5,]              12            1                     26
#>      Erysipelotrichaceae UCG_003 Escherichia Shigella Eubacterium brachy
#> [1,]                           0                    0                  0
#> [2,]                           0                    0                  0
#> [3,]                           0                    0                  0
#> [4,]                           0                    0                  0
#> [5,]                          13                   65                  3
#>      Eubacterium coprostanoligenes Eubacterium nodatum Eubacterium ruminantium
#> [1,]                             0                   0                       0
#> [2,]                             0                   0                       0
#> [3,]                             0                   0                       0
#> [4,]                             0                   1                       0
#> [5,]                            55                   8                       3
#>      Eubacterium ventriosum Ezakiella Facklamia Faecalibacterium Faecalitalea
#> [1,]                      0         0         0                0            0
#> [2,]                      0         0         0                0            0
#> [3,]                      0         0         0                0            0
#> [4,]                      0         0         0                2            0
#> [5,]                     21        11         3              933           14
#>      Fastidiosipila Filifactor Finegoldia Flavobacterium Flavonifractor
#> [1,]              0          0          0              0              0
#> [2,]              0          0          0              0              0
#> [3,]              0          0          0              0              0
#> [4,]              0          0          1              0              0
#> [5,]            731          9         79              2             51
#>      Fretibacterium Fusicatenibacter Fusobacterium Gallicola Gardnerella
#> [1,]              0                0             0         0           0
#> [2,]              0                0             0         0           0
#> [3,]              0                0             0         0           0
#> [4,]              0                0             6         0           0
#> [5,]              3               60           482         3           4
#>      Gemella Gordonia Granulicatella Haemophilus Herbaspirillum Holdemania
#> [1,]       0        0              0           0              0          0
#> [2,]       0        0              0           0              0          0
#> [3,]       1        0              2           3              0          0
#> [4,]      25        0             24          92              0          0
#> [5,]     391        3            141         922             38          4
#>      Howardella Hydrogenoanaerobacterium Intestinibacter Intestinimonas
#> [1,]          0                        0               0              0
#> [2,]          0                        0               0              0
#> [3,]          0                        0               0              0
#> [4,]          0                        0               0              0
#> [5,]          1                        1               9              2
#>      Janibacter Jeotgalicoccus Johnsonella Kingella Lachnoanaerobaculum
#> [1,]          0              0           0        0                   0
#> [2,]          0              0           0        0                   0
#> [3,]          0              0           0        0                   0
#> [4,]          0              0           0        0                   5
#> [5,]          2              1          46       17                  45
#>      Lachnoclostridium Lachnospira Lachnospiraceae FCS020
#> [1,]                 0           0                      0
#> [2,]                 0           0                      0
#> [3,]                 0           0                      0
#> [4,]                 0           0                      0
#> [5,]               142          67                      2
#>      Lachnospiraceae NC2004 Lachnospiraceae ND3007 Lachnospiraceae NK4A136
#> [1,]                      0                      0                       0
#> [2,]                      0                      0                       0
#> [3,]                      0                      0                       0
#> [4,]                      0                      0                       0
#> [5,]                      5                      9                      44
#>      Lachnospiraceae NK4B4 Lachnospiraceae UCG_003 Lachnospiraceae UCG_004
#> [1,]                     0                       0                       0
#> [2,]                     0                       0                       0
#> [3,]                     0                       0                       0
#> [4,]                     0                       0                       0
#> [5,]                     2                     105                       6
#>      Lachnospiraceae UCG_008 Lactobacillus Lactococcus Lautropia Leptotrichia
#> [1,]                       0             0           0         0            0
#> [2,]                       0             0           0         0            0
#> [3,]                       0             0           0         0            0
#> [4,]                       1             1           0         9           13
#> [5,]                      50         16220          19       171           72
#>      Luteibacter Megasphaera Methylobacterium Methyloversatilis Micrococcus
#> [1,]           0           0                0                 0           0
#> [2,]           0           0                0                 0           0
#> [3,]           0           0                0                 0           0
#> [4,]           0           2                0                 0           0
#> [5,]           2         217                3                 2           8
#>      Mitsuokella Mobiluncus Mogibacterium Moraxella Moryella Mycobacterium
#> [1,]           0          0             0         0        0             0
#> [2,]           0          0             0         0        0             0
#> [3,]           0          0             0         0        0             0
#> [4,]           0          0             0         0        0             0
#> [5,]           2          7             6         4        2             6
#>      Mycoplasma Negativicoccus Neisseria Nevskia Odoribacter Oribacterium
#> [1,]          0              0       0.0       0           0            0
#> [2,]          0              0       0.0       0           0            0
#> [3,]          0              0       1.5       0           0            0
#> [4,]          0              0      57.0       0           0            9
#> [5,]          9              6     623.0       1         114          258
#>      Oscillospira Oxalobacter Parabacteroides Paracocccus Paraprevotella
#> [1,]            0           0               0           0              0
#> [2,]            0           0               0           0              0
#> [3,]            0           0               0           0              0
#> [4,]            0           0               2           0              0
#> [5,]            3           5             403           9              7
#>      Parasutterella Parvimonas Pelomonas Peptoclostridium Peptococcus
#> [1,]              0          0         0                0           0
#> [2,]              0          0         0                0           0
#> [3,]              0          0         0                0           0
#> [4,]              0          0         0                0           0
#> [5,]            360         52         2               43           6
#>      Peptoniphilus Peptostreptococcus Phascolarctobacterium Porphyromonas
#> [1,]             0                  0                     0             0
#> [2,]             0                  0                     0             0
#> [3,]             0                  0                     0             0
#> [4,]             1                  1                     0             6
#> [5,]           165                 23                   290            63
#>      Prevotella Prevotella 1 Prevotella 2 Prevotella 6 Prevotella 7
#> [1,]        0.0            0            0            0            0
#> [2,]        0.0            0            0            0            0
#> [3,]        3.5            0            0            0            0
#> [4,]       37.0            0            1            3           42
#> [5,]      290.0           13           28           36          393
#>      Prevotella 9 Prevotellaceae NK3B31 Propionibacterium Propionivibrio
#> [1,]            0                     0                 0              0
#> [2,]            0                     0                 0              0
#> [3,]            0                     0                 0              0
#> [4,]            0                     0                 2              0
#> [5,]            2                    80              4468              1
#>      Providencia Pseudobutyrivibrio Pseudomonas Ralstonia Rhodobacter
#> [1,]           0                  0           0         0           0
#> [2,]           0                  0           0         0           0
#> [3,]           0                  0           0         0           0
#> [4,]           0                  1           0         0           0
#> [5,]           1                346          18         4           1
#>      Rikenellaceae RC9_gut Roseburia Rothia Ruminiclostridium
#> [1,]                     0         0      0                 0
#> [2,]                     0         0      0                 0
#> [3,]                     0         0      0                 0
#> [4,]                     0         0      6                 0
#> [5,]                    49         2    107                 3
#>      Ruminiclostridium 1 Ruminiclostridium 5 Ruminiclostridium 6
#> [1,]                   0                   0                   0
#> [2,]                   0                   0                   0
#> [3,]                   0                   0                   0
#> [4,]                   0                   0                   0
#> [5,]                   1                  37                  13
#>      Ruminiclostridium 9 Ruminococcaceae NK4A214 Ruminococcaceae UCG_002
#> [1,]                   0                       0                       0
#> [2,]                   0                       0                       0
#> [3,]                   0                       0                       0
#> [4,]                   0                       0                       0
#> [5,]                  30                       4                     132
#>      Ruminococcaceae UCG_003 Ruminococcaceae UCG_005 Ruminococcaceae UCG_010
#> [1,]                       0                       0                       0
#> [2,]                       0                       0                       0
#> [3,]                       0                       0                       0
#> [4,]                       0                       0                       0
#> [5,]                      69                       9                      35
#>      Ruminococcaceae UCG_013 Ruminococcaceae UCG_014 Ruminococcus 1
#> [1,]                       0                       0              0
#> [2,]                       0                       0              0
#> [3,]                       0                       0              0
#> [4,]                       0                       5              0
#> [5,]                      34                     193            186
#>      Ruminococcus 2 Scardovia Selenomonas Selenomonas 3 Selenomonas 4
#> [1,]              0         0           0             0             0
#> [2,]              0         0           0             0             0
#> [3,]              0         0           0             0             0
#> [4,]              0         0           0             2             0
#> [5,]             47         2           3           123             6
#>      Senegalimassilia Shuttleworthia Sneathia Solobacterium Sphingobium
#> [1,]                0              0        0             0           0
#> [2,]                0              0        0             0           0
#> [3,]                0              0        0             0           0
#> [4,]                0              0        0             0           0
#> [5,]              100           1654       90             2           3
#>      Sphingomonas Sphingopyxis Sporolactobacillus Staphylococcus
#> [1,]            0            0                  0              0
#> [2,]            0            0                  0              0
#> [3,]            0            0                  0              0
#> [4,]            0            0                  0              1
#> [5,]            4            7                  3           9581
#>      Stenotrophomonas Stomatobaculum Streptobacillus Streptococcus
#> [1,]                0              0               0             0
#> [2,]                0              0               0             0
#> [3,]                0              0               0            27
#> [4,]                0              1               0           411
#> [5,]                2             54               2          6308
#>      Subdoligranulum Sutterella Tannerella Thalassospira Treponema 2
#> [1,]               0          0          0             0           0
#> [2,]               0          0          0             0           0
#> [3,]               0          0          0             0           0
#> [4,]               1          0          0             0           0
#> [5,]              40         27         18           170          49
#>      Tumebacillus Turicibacter Tyzzerella Unc. BactOral Unc. FctCyli5
#> [1,]            0            0          0             0             0
#> [2,]            0            0          0             0             0
#> [3,]            0            0          0             0             0
#> [4,]            0            0          0             0             0
#> [5,]            3            2          4             6             3
#>      Unc. FirmOra2 Unc. GcbBacte Unc. Hu4Lup30 Unc. HunGu114 Unc. PinJeffr
#> [1,]             0             0             0             0             0
#> [2,]             0             0             0             0             0
#> [3,]             0             0             0             0             0
#> [4,]             0             0             0             0             0
#> [5,]             6            14            93             5            11
#>      Unc. Unc006vd Unc. Unc00a1s Unc. Unc00g7e Unc. Unc00lvf Unc. Unc00zke
#> [1,]             0             0             0             0             0
#> [2,]             0             0             0             0             0
#> [3,]             0             0             0             0             0
#> [4,]             0             0             0             0             0
#> [5,]            24             4             1           898             1
#>      Unc. Unc0168y Unc. Unc01skg Unc. Unc021pc Unc. Unc0266w Unc. Unc029kz
#> [1,]             0             0             0             0             0
#> [2,]             0             0             0             0             0
#> [3,]             0             0             0             0             0
#> [4,]             0             0             0             0             0
#> [5,]             2             1             3          2992             7
#>      Unc. Unc02cxe Unc. Unc02jgy Unc. Unc02oth Unc. Unc02vkz Unc. Unc037e3
#> [1,]             0             0             0             0             0
#> [2,]             0             0             0             0             0
#> [3,]             0             0             0             0             0
#> [4,]             0             0             0             0             0
#> [5,]             2             1            15             9             2
#>      Unc. Unc038h1 Unc. Unc038t2 Unc. Unc0397y Unc. Unc039sp Unc. Unc03dfb
#> [1,]             0             0             0             0             0
#> [2,]             0             0             0             0             0
#> [3,]             0             0             0             0             0
#> [4,]             0             0             0             0             0
#> [5,]            11            16            93           120             1
#>      Unc. Unc03f4o Unc. Unc03iw8 Unc. Unc03j4k Unc. Unc03m1r Unc. Unc03php
#> [1,]             0             0             0             0             0
#> [2,]             0             0             0             0             0
#> [3,]             0             0             0             0             0
#> [4,]             0             0             0             0             0
#> [5,]             3             5             5             4            15
#>      Unc. Unc03tng Unc. Unc03to4 Unc. Unc0406b Unc. Unc0434x Unc. Unc04cb2
#> [1,]             0             0             0             0             0
#> [2,]             0             0             0             0             0
#> [3,]             0             0             0             0             0
#> [4,]             0             0             0             0             0
#> [5,]             3             1             4            26             9
#>      Unc. Unc13404 Unc. Unc23927 Unc. Unc28979 Unc. Unc29437 Unc. Unc48787
#> [1,]             0             0             0             0             0
#> [2,]             0             0             0             0             0
#> [3,]             0             0             0             0             0
#> [4,]             0             0             0             0             0
#> [5,]             2             2            34             1             1
#>      Unc. Unc58411 Unc. Unc61746 Unc. Unc66966 Unc. Unc71514 Unc. Unc72229
#> [1,]             0             0             0             0             0
#> [2,]             0             0             0             0             0
#> [3,]             0             0             0             0             0
#> [4,]             0             0             0             0             0
#> [5,]             2           923            12             4             1
#>      Unc. Unc75035 Unc. Unc79323 Unc. Unc82785 Unc. Unc89065 Unc. Unc90828
#> [1,]             0             0             0             0             0
#> [2,]             0             0             0             0             0
#> [3,]             0             0             0             0             0
#> [4,]             0             0             0             0             0
#> [5,]             2             1             2             3            71
#>      Unc. UncO2854 Veillonella
#> [1,]             0         0.0
#> [2,]             0         0.0
#> [3,]             0         1.5
#> [4,]             0        90.0
#> [5,]            66      1223.0
```
