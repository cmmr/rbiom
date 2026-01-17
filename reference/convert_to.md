# Convert biom data to an external package class.

Requires the relevant Bioconductor R package to be installed:

- `convert_to_phyloseq` - :

  [phyloseq](https://bioconductor.org/packages/phyloseq/)

- `convert_to_animalcules` - :

  [animalcules](https://bioconductor.org/packages/animalcules/)

- `convert_to_SE` - :

  [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/)

- `convert_to_TSE` - :

  [TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment/)

## Usage

``` r
convert_to_SE(biom, ...)

convert_to_TSE(biom, ...)

convert_to_phyloseq(biom, ...)

convert_to_animalcules(biom, ...)
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

A phyloseq, animalcules (MultiAssayExperiment), SummarizedExperiment, or
TreeSummarizedExperiment object.

## Details

animalcules and SummarizedExperiment objects includes counts, metadata,
and taxonomy.

phyloseq and TreeSummarizedExperiment additionally includes the tree and
sequences.

## Examples

``` r
if (FALSE) { # \dontrun{
    library(rbiom)

    print(hmp50)

    # Requires 'phyloseq', a Bioconductor R package
    if (nzchar(system.file(package = "phyloseq"))) {
      physeq <- convert_to_phyloseq(hmp50)
      print(physeq)
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
