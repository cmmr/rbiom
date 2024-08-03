
#' Convert an rbiom object to a SummarizedExperiment object.
#' 
#' Requires the relevant Bioconductor R package to be installed:
#' \describe{
#'   \item{`convert_to_SE` - }{ [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/) }
#'   \item{`convert_to_TSE` - }{ [TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment/) }
#' }
#' 
#' A SummarizedExperiment object includes counts, metadata, and taxonomy.
#' 
#' TreeSummarizedExperiment additionally includes the tree and sequences.
#' 
#' @inherit documentation_default
#' 
#' @return A SummarizedExperiment or TreeSummarizedExperiment object.
#' 
#' @export
#' @examples
#'     
#'     library(rbiom) 
#'     
#'     print(hmp50)
#'       
#'     # Requires 'SummarizedExperiment' Bioconductor R package
#'     if (nzchar(system.file(package = "SummarizedExperiment"))) {
#'       se <- convert_to_SE(hmp50)
#'       print(se)
#'     }
#'       
#'     # Requires 'TreeSummarizedExperiment' Bioconductor R package
#'     if (nzchar(system.file(package = "TreeSummarizedExperiment"))) {
#'       tse <- convert_to_TSE(hmp50)
#'       print(tse)
#'     }

convert_to_SE <- function (biom) {
  
  biom <- as_rbiom(biom)
  
  if (nchar(system.file(package = "SummarizedExperiment")) == 0)
    stop("Bioconductor R package 'SummarizedExperiment' must be installed to use convert_to_SE().")
  
  SummarizedExperiment::SummarizedExperiment(
    assays  = list('OTU table' = as.matrix(biom$counts)),
    rowData = biom$taxonomy %>% tibble::column_to_rownames(".otu") %>% as.matrix(),
    colData = biom$metadata %>% tibble::column_to_rownames(".sample") )
}



#' @rdname convert_to_SE
#' @export

convert_to_TSE <- function (biom) {
  
  biom <- as_rbiom(biom)
  
  if (nchar(system.file(package = "TreeSummarizedExperiment")) == 0)
    stop("Bioconductor R package 'TreeSummarizedExperiment' must be installed to use convert_to_TSE().")
  
  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays       = list('OTU table' = as.matrix(biom$counts)),
    rowData      = biom$taxonomy %>% tibble::column_to_rownames(".otu") %>% as.matrix(),
    colData      = biom$metadata %>% tibble::column_to_rownames(".sample"),
    rowTree      = biom$tree, 
    referenceSeq = Biostrings::DNAStringSet(biom$sequences) )
}

