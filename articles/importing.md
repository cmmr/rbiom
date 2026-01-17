# Importing Datasets

## Create an rbiom object

The fastest way to make an rbiom object is with
[`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md),
which accepts:

- A filepath or URL to a BIOM file.
- An abundance matrix with OTUs in rows and samples in columns.
- A `phyloseq`-class object, from the phyloseq Bioconductor R package.
- A list with `counts` and optionally `metadata`, `taxonomy`, `tree`,
  etc (see
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md)).

``` r
library(rbiom)

# create a simple matrix ------------------------
mtx <- matrix(
  data     = floor(runif(24) * 1000), 
  nrow     = 6, 
  dimnames = list(paste0("OTU", 1:6), paste0("Sample", 1:4)) )
mtx
#>      Sample1 Sample2 Sample3 Sample4
#> OTU1      80     497      34     388
#> OTU2     834     289     320     975
#> OTU3     600     732     402     289
#> OTU4     157     772     195     678
#> OTU5       7     874     403     735
#> OTU6     466     174      63     195

# convert matrix to rbiom -----------------------
biom <- as_rbiom(biom = mtx)
biom
#> 
#> ══ Untitled Dataset ════════════════════════════════════════
#> 
#>       4 Samples: Sample1, Sample2, Sample3, and Sample4
#>       6 OTUs:    OTU1, OTU2, OTU3, OTU4, OTU5, and OTU6
#>       1 Ranks:   .otu
#>       1 Fields:  .sample
#>         Tree:    <absent>
#> 
#> ── 1.4k - 3.3k reads/sample ────────────────── 2026-01-17 ──
#> 

# convert from phyloseq to rbiom ----------------
file <- system.file("extdata", "rich_sparse_otu_table.biom", package="phyloseq")
phy  <- phyloseq::import_biom(file)
phy
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 5 taxa and 6 samples ]
#> sample_data() Sample Data:       [ 6 samples by 4 sample variables ]
#> tax_table()   Taxonomy Table:    [ 5 taxa by 7 taxonomic ranks ]

biom <- as_rbiom(biom = phy)
biom
#> 
#> ══ Imported phyloseq Data ══════════════════════════════════
#> 
#>       6 Samples: Sample1, Sample2, ..., and Sample6
#>       5 OTUs:    GG_OTU_1, GG_OTU_2, GG_OTU_3, ...
#>       8 Ranks:   .otu, Rank1, Rank2, ..., and Rank7
#>       5 Fields:  .sample, BarcodeSequence, ...
#>         Tree:    <absent>
#> 
#> ── 3 - 7 reads/sample ──────────────────────── 2026-01-17 ──
#> 
```

Now we have `biom`, an rbiom-class object that can be used with this
package’s functions. If you loaded your data from a BIOM file or
phyloseq object, it might already include metadata, ranks, and a tree.
These attributes are technically optional. However, more analyses are
possible when extra information about samples and OTUs are present.

See the \[QIIME 2, mothur, and BioConductor\]\[convert\] page for
instructions on how to import data from those packages.

## Attach metadata

`$metadata` lets you set arbitrary data for each sample.

A few quick rules:

- `.sample` should be the first column.
- Other column names cannot start with a dot (`.`).
- Sample names need to match `biom$samples`.

``` r
# create example metadata -----------------------
md <- data.frame(
  .sample   = paste0("Sample", 1:4),
  state     = c("TX", "TX", "WA", "WA"),
  age       = c(32, 19, 36, 40),
  treatment = c(1, 2, 1, 2) )
md
#>   .sample state age treatment
#> 1 Sample1    TX  32         1
#> 2 Sample2    TX  19         2
#> 3 Sample3    WA  36         1
#> 4 Sample4    WA  40         2

# add metadata to rbiom object ------------------
biom <- as_rbiom(biom = mtx)
biom$metadata <- md
biom
#> 
#> ══ Untitled Dataset ════════════════════════════════════════
#> 
#>       4 Samples: Sample1, Sample2, Sample3, and Sample4
#>       6 OTUs:    OTU1, OTU2, OTU3, OTU4, OTU5, and OTU6
#>       1 Ranks:   .otu
#>       4 Fields:  .sample, state, age, and treatment
#>         Tree:    <absent>
#> 
#> ── 1.4k - 3.3k reads/sample ────────────────── 2026-01-17 ──
#> 

# or in a single step ---------------------------
biom <- as_rbiom(biom = list(counts = mtx, metadata = md))
biom
#> 
#> ══ Untitled Dataset ════════════════════════════════════════
#> 
#>       4 Samples: Sample1, Sample2, Sample3, and Sample4
#>       6 OTUs:    OTU1, OTU2, OTU3, OTU4, OTU5, and OTU6
#>       1 Ranks:   .otu
#>       4 Fields:  .sample, state, age, and treatment
#>         Tree:    <absent>
#> 
#> ── 1.4k - 3.3k reads/sample ────────────────── 2026-01-17 ──
#> 
```

### Setting categorical variables

Any categorical metadata variable that looks numerical, such as
“treatment” in the above example, will need to be manually changed to a
categorical variable.

``` r
class(pull(biom, 'treatment'))
#> [1] "numeric"

biom$metadata$treatment %<>% as.factor()

class(pull(biom, 'treatment'))
#> [1] "factor"
pull(biom, 'treatment')
#> Sample1 Sample2 Sample3 Sample4 
#>       1       2       1       2 
#> Levels: 1 2
```

## Attach a tree

Use `$tree` to set the tree. You can specify a phylo object directly, or
a newick file/string.

``` r
# define a random tree --------------------------
biom$tree <- "(((OTU6,(OTU5,OTU4)),OTU3),(OTU2,OTU1));"
biom
#> 
#> ══ Untitled Dataset ════════════════════════════════════════
#> 
#>       4 Samples: Sample1, Sample2, Sample3, and Sample4
#>       6 OTUs:    OTU1, OTU2, OTU3, OTU4, OTU5, and OTU6
#>       1 Ranks:   .otu
#>       4 Fields:  .sample, state, age, and treatment
#>         Tree:    <present>
#> 
#> ── 1.4k - 3.3k reads/sample ────────────────── 2026-01-17 ──
#> 
```

## Attach taxonomy

Use `$taxonomy` to define taxonomic clades for each OTU.

``` r
# .otu must match otu_names(biom) ---------------
map <- data.frame(
  .otu   = paste0("OTU", 1:6),
  Phylum = c("Bacteroidetes", "Firmicutes", "Firmicutes"),
  Order  = c("Bacteroidia", "Clostridiales", "Bacillales") )
map
#>   .otu        Phylum         Order
#> 1 OTU1 Bacteroidetes   Bacteroidia
#> 2 OTU2    Firmicutes Clostridiales
#> 3 OTU3    Firmicutes    Bacillales
#> 4 OTU4 Bacteroidetes   Bacteroidia
#> 5 OTU5    Firmicutes Clostridiales
#> 6 OTU6    Firmicutes    Bacillales

biom$taxonomy <- map
biom
#> 
#> ══ Untitled Dataset ════════════════════════════════════════
#> 
#>       4 Samples: Sample1, Sample2, Sample3, and Sample4
#>       6 OTUs:    OTU1, OTU2, OTU3, OTU4, OTU5, and OTU6
#>       3 Ranks:   .otu, Phylum, and Order
#>       4 Fields:  .sample, state, age, and treatment
#>         Tree:    <present>
#> 
#> ── 1.4k - 3.3k reads/sample ────────────────── 2026-01-17 ──
#> 
```
