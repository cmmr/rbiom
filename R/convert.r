
#' Convert biom data to an external package class
#'
#' Converts your `rbiom` object into other common Bioconductor data structures.
#' Each function requires the corresponding target package to be installed.
#'
#' @details
#' * **`convert_to_animalcules()`**: Converts to a `MultiAssayExperiment` object tailored for the
#'   [animalcules](https://bioconductor.org/packages/animalcules/) interactive microbiome analysis toolkit.
#'   *Includes: counts, metadata, and taxonomy.*
#'
#' * **`convert_to_biomformat()`**: Converts to a `biom` object used by the
#'   [biomformat](https://bioconductor.org/packages/biomformat/) package, the standard Bioconductor
#'   class for reading and writing BIOM data.
#'   *Includes: counts, metadata, and taxonomy.*
#'
#' * **`convert_to_phyloseq()`**: Converts to a `phyloseq` object for use with the comprehensive
#'   [phyloseq](https://bioconductor.org/packages/phyloseq/) ecosystem. Rank names, metadata
#'   field names, and metadata column types are preserved. Sequences are included only
#'   when every OTU has one (requires the Biostrings package).
#'   *Includes: counts, metadata, taxonomy, phylogenetic tree, and sequences.*
#'
#' * **`convert_to_SE()`**: Converts to a `SummarizedExperiment` object, a core
#'   [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/) Bioconductor
#'   container for matrix-like data and annotations.
#'   *Includes: counts, metadata, and taxonomy.*
#'
#' * **`convert_to_TSE()`**: Converts to a `TreeSummarizedExperiment` object. This extends the SE class
#'   to natively support hierarchical [TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment/)
#'   relationships. The tree and sequences are included when present.
#'   *Includes: counts, metadata, taxonomy, phylogenetic tree, and sequences.*
#'
#' @name convert_to
#' @inherit documentation_default
#'
#' @param ...  Not Used.
#'
#' @return An `animalcules` (`MultiAssayExperiment` class), `biomformat` (`biom` class), `phyloseq`,
#' `SummarizedExperiment`, or `TreeSummarizedExperiment` object.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     print(hmp50)
#'
#'     \donttest{
#'     # Requires 'animalcules' dependencies, Bioconductor R packages
#'     if (requireNamespace("MultiAssayExperiment", quietly = TRUE) &&
#'         requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'       ani <- convert_to_animalcules(hmp50)
#'       print(ani)
#'     }
#'
#'     # Requires 'biomformat', a Bioconductor R package
#'     if (requireNamespace("biomformat", quietly = TRUE)) {
#'       bio <- convert_to_biomformat(hmp50)
#'       print(bio)
#'     }
#'
#'     # Requires 'phyloseq', a Bioconductor R package
#'     if (requireNamespace("phyloseq", quietly = TRUE)) {
#'       phy <- convert_to_phyloseq(hmp50)
#'       print(phy)
#'     }
#'
#'     # Requires 'SummarizedExperiment', a Bioconductor R package
#'     if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'       se <- convert_to_SE(hmp50)
#'       print(se)
#'     }
#'
#'     # Requires 'TreeSummarizedExperiment', a Bioconductor R package
#'     if (requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
#'       tse <- convert_to_TSE(hmp50)
#'       print(tse)
#'     }
#'     }

convert_to_animalcules <- function (biom, ...) {

  require_package('S4Vectors',            'to use convert_to_animalcules()')
  require_package('SummarizedExperiment', 'to use convert_to_animalcules()')
  require_package('MultiAssayExperiment', 'to use convert_to_animalcules()')

  biom <- as_rbiom(biom)

  assays <- S4Vectors::SimpleList('MGX' = as.matrix(biom$counts))

  colData <- metadata_df(biom)
  for (i in which(vapply(colData, is.factor, TRUE)))
    colData[[i]] <- as.character(colData[[i]])
  colData <- S4Vectors::DataFrame(colData, check.names = FALSE)

  rowData <- S4Vectors::DataFrame(taxonomy_df(biom), check.names = FALSE)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = assays,
    colData = colData,
    rowData = rowData )

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList('MicrobeGenetics' = se),
    colData     = colData )
}



#' @rdname convert_to
#' @export

convert_to_biomformat <- function (biom, ...) {

  require_package('biomformat', 'to use convert_to_biomformat()')

  biom <- as_rbiom(biom)

  biomformat::make_biom(
    data                 = as.matrix(biom$counts),
    sample_metadata      = if (biom$n_fields > 1) metadata_df(biom),
    observation_metadata = if (biom$n_ranks  > 1) taxonomy_df(biom),
    id                   = biom$id )
}



#' @rdname convert_to
#' @export

convert_to_phyloseq <- function (biom, ...) {

  require_package('phyloseq', 'to use convert_to_phyloseq()')

  biom <- as_rbiom(biom)

  parts <- list(phyloseq::otu_table(as.matrix(biom$counts), taxa_are_rows = TRUE))

  if (biom$n_fields > 1)
    parts <- c(parts, list(phyloseq::sample_data(metadata_df(biom))))

  if (biom$n_ranks > 1)
    parts <- c(parts, list(phyloseq::tax_table(as.matrix(taxonomy_df(biom)))))

  if (!is.null(biom$tree))
    parts <- c(parts, list(phyloseq::phy_tree(biom$tree)))

  if (!is.null(seqs <- biom$sequences)) {
    if (anyNA(seqs)) {
      cli_warn("Sequences are omitted from the phyloseq object because some OTUs have no sequence.")
    } else {
      require_package('Biostrings', 'to include sequences in a phyloseq object')
      parts <- c(parts, list(phyloseq::refseq(Biostrings::DNAStringSet(seqs))))
    }
  }

  # phyloseq() returns a lone component as-is rather than wrapping it.
  if (length(parts) == 1)
    return (methods::new("phyloseq", otu_table = parts[[1]]))

  do.call(phyloseq::phyloseq, parts)
}



#' @rdname convert_to
#' @export

convert_to_SE <- function (biom, ...) {

  require_package('S4Vectors',            'to use convert_to_SE()')
  require_package('SummarizedExperiment', 'to use convert_to_SE()')

  biom <- as_rbiom(biom)

  SummarizedExperiment::SummarizedExperiment(
    assays  = list('OTU table' = as.matrix(biom$counts)),
    rowData = S4Vectors::DataFrame(taxonomy_df(biom), check.names = FALSE),
    colData = S4Vectors::DataFrame(metadata_df(biom), check.names = FALSE) )
}



#' @rdname convert_to
#' @export

convert_to_TSE <- function (biom, ...) {

  require_package('S4Vectors',                'to use convert_to_TSE()')
  require_package('TreeSummarizedExperiment', 'to use convert_to_TSE()')

  biom <- as_rbiom(biom)

  refseq <- NULL
  if (!is.null(seqs <- biom$sequences) && !anyNA(seqs)) {
    require_package('Biostrings', 'to include sequences in a TreeSummarizedExperiment object')
    refseq <- Biostrings::DNAStringSet(seqs)
  }

  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays       = list('OTU table' = as.matrix(biom$counts)),
    rowData      = S4Vectors::DataFrame(taxonomy_df(biom), check.names = FALSE),
    colData      = S4Vectors::DataFrame(metadata_df(biom), check.names = FALSE),
    rowTree      = biom$tree,
    referenceSeq = refseq )
}



#________________________________________________________
# Base data.frames with sample / OTU names as row names.
#________________________________________________________
metadata_df <- function (biom) {
  df <- as.data.frame(biom$metadata[,-1,drop=FALSE], check.names = FALSE, stringsAsFactors = FALSE)
  rownames(df) <- biom$samples
  return (df)
}

taxonomy_df <- function (biom) {
  df <- as.data.frame(biom$taxonomy[,-1,drop=FALSE], check.names = FALSE, stringsAsFactors = FALSE)
  for (i in seq_along(df)) df[[i]] <- as.character(df[[i]])
  rownames(df) <- biom$otus
  return (df)
}
