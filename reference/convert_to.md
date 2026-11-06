# Convert biom data to an external package class

Converts your `rbiom` object into other common Bioconductor data
structures. Each function requires the corresponding target package to
be installed.

## Usage

``` r
convert_to_animalcules(biom, ...)

convert_to_biomformat(biom, ...)

convert_to_phyloseq(biom, ...)

convert_to_SE(biom, ...)

convert_to_TSE(biom, ...)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- ...:

  Not Used.

## Value

An `animalcules` (`MultiAssayExperiment` class), `biomformat` (`biom`
class), `phyloseq`, `SummarizedExperiment`, or
`TreeSummarizedExperiment` object.

## Details

- **`convert_to_animalcules()`**: Converts to a `MultiAssayExperiment`
  object tailored for the
  [animalcules](https://bioconductor.org/packages/animalcules/)
  interactive microbiome analysis toolkit. *Includes: counts, metadata,
  and taxonomy.*

- **`convert_to_biomformat()`**: Converts to a `biom` object used by the
  [biomformat](https://bioconductor.org/packages/biomformat/) package,
  the standard Bioconductor class for reading and writing BIOM data.
  *Includes: counts, metadata, and taxonomy.*

- **`convert_to_phyloseq()`**: Converts to a `phyloseq` object for use
  with the comprehensive
  [phyloseq](https://bioconductor.org/packages/phyloseq/) ecosystem.
  Rank names, metadata field names, and metadata column types are
  preserved. Sequences are included only when every OTU has one
  (requires the Biostrings package). *Includes: counts, metadata,
  taxonomy, phylogenetic tree, and sequences.*

- **`convert_to_SE()`**: Converts to a `SummarizedExperiment` object, a
  core
  [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/)
  Bioconductor container for matrix-like data and annotations.
  *Includes: counts, metadata, and taxonomy.*

- **`convert_to_TSE()`**: Converts to a `TreeSummarizedExperiment`
  object. This extends the SE class to natively support hierarchical
  [TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment/)
  relationships. The tree and sequences are included when present.
  *Includes: counts, metadata, taxonomy, phylogenetic tree, and
  sequences.*

## Examples

``` r
    library(rbiom)

    print(hmp50)
#> 
#> ══ Human Microbiome Project - 50 Sample Demo ═══════════════
#> 
#> Oral, nasal, vaginal, and fecal samples from a diverse set
#> of healthy volunteers. Source: Human Microbiome Project
#> (<https://hmpdacc.org>).
#> 
#>      50 Samples: HMP01, HMP02, HMP03, ..., and HMP50 
#>     490 OTUs:    Unc01yki, Unc53100, LtbAci52, ... 
#>       7 Ranks:   .otu, Kingdom, Phylum, ..., and Genus 
#>       5 Fields:  .sample, Age, BMI, Body Site, and Sex 
#>         Tree:    <present>
#> 
#> ── 182 - 22k reads/sample ──────────────────── 2023-09-22 ──
#> 

    # \donttest{
    # Requires 'animalcules' dependencies, Bioconductor R packages
    if (requireNamespace("MultiAssayExperiment", quietly = TRUE) &&
        requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      ani <- convert_to_animalcules(hmp50)
      print(ani)
    }
#> A MultiAssayExperiment object of 1 listed
#>  experiment with a user-defined name and respective class.
#>  Containing an ExperimentList class object of length 1:
#>  [1] MicrobeGenetics: SummarizedExperiment with 490 rows and 50 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files

    # Requires 'biomformat', a Bioconductor R package
    if (requireNamespace("biomformat", quietly = TRUE)) {
      bio <- convert_to_biomformat(hmp50)
      print(bio)
    }
#> biom object. 
#> type: OTU table 
#> matrix_type: dense 
#> 490 rows and 50 columns 

    # Requires 'phyloseq', a Bioconductor R package
    if (requireNamespace("phyloseq", quietly = TRUE)) {
      phy <- convert_to_phyloseq(hmp50)
      print(phy)
    }
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 490 taxa and 50 samples ]
#> sample_data() Sample Data:       [ 50 samples by 4 sample variables ]
#> tax_table()   Taxonomy Table:    [ 490 taxa by 6 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 490 tips and 489 internal nodes ]
#> refseq()      DNAStringSet:      [ 490 reference sequences ]

    # Requires 'SummarizedExperiment', a Bioconductor R package
    if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      se <- convert_to_SE(hmp50)
      print(se)
    }
#> class: SummarizedExperiment 
#> dim: 490 50 
#> metadata(0):
#> assays(1): OTU table
#> rownames(490): Unc01yki Unc53100 ... UncOr431 UncTr598
#> rowData names(6): Kingdom Phylum ... Family Genus
#> colnames(50): HMP01 HMP02 ... HMP49 HMP50
#> colData names(4): Age BMI Body Site Sex

    # Requires 'TreeSummarizedExperiment', a Bioconductor R package
    if (requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
      tse <- convert_to_TSE(hmp50)
      print(tse)
    }
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> class: TreeSummarizedExperiment 
#> dim: 490 50 
#> metadata(0):
#> assays(1): OTU table
#> rownames(490): Unc01yki Unc53100 ... UncOr431 UncTr598
#> rowData names(6): Kingdom Phylum ... Family Genus
#> colnames(50): HMP01 HMP02 ... HMP49 HMP50
#> colData names(4): Age BMI Body Site Sex
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> reducedDimNames(0):
#> mainExpName: NULL
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> altExpNames(0):
#> rowLinks: a LinkDataFrame (490 rows)
#> rowTree: 1 phylo tree(s) (490 leaves)
#> colLinks: NULL
#> colTree: NULL
#> referenceSeq: a DNAStringSet (490 sequences)
    # }
```
