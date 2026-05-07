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
  *Includes: counts, metadata, taxonomy, phylogenetic tree, and
  sequences.*

- **`convert_to_SE()`**: Converts to a `SummarizedExperiment` object, a
  core
  [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/)
  Bioconductor container for matrix-like data and annotations.
  *Includes: counts, metadata, and taxonomy.*

- **`convert_to_TSE()`**: Converts to a `TreeSummarizedExperiment`
  object. This extends the SE class to natively support hierarchical
  [TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment/)
  relationships. *Includes: counts, metadata, taxonomy, phylogenetic
  tree, and sequences.*

## Examples

``` r
if (FALSE) { # \dontrun{
    library(rbiom)

    print(hmp50)

    # Requires 'animalcules', a Bioconductor R package
    if (nzchar(system.file(package = "animalcules"))) {
      ani <- convert_to_animalcules(hmp50)
      print(ani)
    }

    # Requires 'biomformat', a Bioconductor R package
    if (nzchar(system.file(package = "biomformat"))) {
      bio <- convert_to_biomformat(hmp50)
      print(bio)
    }

    # Requires 'phyloseq', a Bioconductor R package
    if (nzchar(system.file(package = "phyloseq"))) {
      phy <- convert_to_phyloseq(hmp50)
      print(phy)
    }

    # Requires 'SummarizedExperiment', a Bioconductor R package
    if (nzchar(system.file(package = "SummarizedExperiment"))) {
      se <- convert_to_SE(hmp50)
      print(se)
    }

    # Requires 'TreeSummarizedExperiment', a Bioconductor R package
    if (nzchar(system.file(package = "TreeSummarizedExperiment"))) {
      tse <- convert_to_TSE(hmp50)
      print(tse)
    }
} # }
```
