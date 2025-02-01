
#' Convert biom data to an external package class.
#'
#' Requires the relevant Bioconductor R package to be installed:
#' \describe{
#'   \item{`convert_to_phyloseq` - }{ [phyloseq](https://bioconductor.org/packages/phyloseq/) }
#'   \item{`convert_to_SE` - }{ [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/) }
#'   \item{`convert_to_TSE` - }{ [TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment/) }
#' }
#'
#' A SummarizedExperiment object includes counts, metadata, and taxonomy.
#'
#' phyloseq and TreeSummarizedExperiment additionally includes the tree and sequences.
#'
#' @name convert_to
#' @inherit documentation_default
#'
#' @param ...  Not Used.
#'
#' @return A phyloseq, SummarizedExperiment, or TreeSummarizedExperiment object.
#'
#' @export
#' @examples
#' \dontrun{
#'     library(rbiom)
#'
#'     print(hmp50)
#'
#'     # Requires 'phyloseq', a Bioconductor R package
#'     if (nzchar(system.file(package = "phyloseq"))) {
#'       physeq <- convert_to_phyloseq(hmp50)
#'       print(physeq)
#'     }
#'
#'     # Requires 'SummarizedExperiment', a Bioconductor R package
#'     if (nzchar(system.file(package = "SummarizedExperiment"))) {
#'       se <- convert_to_SE(hmp50)
#'       print(se)
#'     }
#'
#'     # Requires 'TreeSummarizedExperiment', a Bioconductor R package
#'     if (nzchar(system.file(package = "TreeSummarizedExperiment"))) {
#'       tse <- convert_to_TSE(hmp50)
#'       print(tse)
#'     }
#' }

convert_to_SE <- function (biom, ...) {

  require_package('SummarizedExperiment', 'to use convert_to_SE()')
 
  dots <- list(...)
  biom <- as_rbiom(biom)
 
  SummarizedExperiment <- getFromNamespace('SummarizedExperiment', 'SummarizedExperiment')
  SummarizedExperiment(
    assays  = list('OTU table' = as.matrix(biom$counts)),
    rowData = biom$taxonomy %>% tibble::column_to_rownames(".otu") %>% as.matrix(),
    colData = biom$metadata %>% tibble::column_to_rownames(".sample") )
  
}



#' @rdname convert_to
#' @export

convert_to_TSE <- function (biom, ...) {

  require_package('TreeSummarizedExperiment', 'to use convert_to_TSE()')

  dots <- list(...)
  biom <- as_rbiom(biom)
  
  TreeSummarizedExperiment <- getFromNamespace('TreeSummarizedExperiment', 'TreeSummarizedExperiment')
  DNAStringSet             <- getFromNamespace('DNAStringSet', 'Biostrings')
  
  TreeSummarizedExperiment(
    assays       = list('OTU table' = as.matrix(biom$counts)),
    rowData      = biom$taxonomy %>% tibble::column_to_rownames(".otu") %>% as.matrix(),
    colData      = biom$metadata %>% tibble::column_to_rownames(".sample"),
    rowTree      = biom$tree,
    referenceSeq = DNAStringSet(biom$sequences) )
  
}



#' @rdname convert_to
#' @export

convert_to_phyloseq <- function (biom, ...) {
  
  require_package('phyloseq', 'to use convert_to_phyloseq()')
  
  dots <- list(...)
  biom <- as_rbiom(biom)
  
  import_biom <- getFromNamespace('import_biom', 'phyloseq')
  
  BIOMfilename   <- write_biom(biom, file = tempfile())
  treefilename   <- if (!is.null(biom$tree))      write_tree(biom,  file = tempfile())
  refseqfilename <- if (!is.null(biom$sequences)) write_fasta(biom, file = tempfile())
  
  on.exit(unlink(c(BIOMfilename, treefilename, refseqfilename)), add = TRUE)
  
  import_biom(
    BIOMfilename   = BIOMfilename, 
    treefilename   = treefilename, 
    refseqfilename = refseqfilename )
}


