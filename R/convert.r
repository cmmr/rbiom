
#  
#  |NOTE|  Commented out 1/7/2025 due to build failure of matrixStats.
#  |NOTE|  TreeSummarizedExperiment -> SummarizedExperiment -> MatrixGenerics -> matrixStats
#  |NOTE|  
#  |NOTE|  Also, this little bit of code requires a large dependency base:
#  |NOTE|  Biobase, BiocGenerics, BiocParallel, Biostrings, DelayedArray, GenomeInfoDb, 
#  |NOTE|  GenomeInfoDbData, GenomicRanges, IRanges, Matrix, MatrixGenerics, S4Arrays, S4Vectors, 
#  |NOTE|  SingleCellExperiment, SparseArray, SummarizedExperiment, TreeSummarizedExperiment, 
#  |NOTE|  UCSC.utils, XVector, matrixStats, tidytree
#  |NOTE|  
#  |NOTE|  If re-enabling, add Biostrings, SummarizedExperiment, and TreeSummarizedExperiment
#  |NOTE|  to DESCRIPTION->Suggests, export(`convert_to_TSE`); export(`convert_to_SE`) to NAMESPACE,
#  |NOTE|  and uncomment tests/testthat/test-convert.r
#  
#  
#  #' Convert an rbiom object to a SummarizedExperiment object.
#  #' 
#  #' Requires the relevant Bioconductor R package to be installed:
#  #' \describe{
#  #'   \item{`convert_to_SE` - }{ [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/) }
#  #'   \item{`convert_to_TSE` - }{ [TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment/) }
#  #' }
#  #' 
#  #' A SummarizedExperiment object includes counts, metadata, and taxonomy.
#  #' 
#  #' TreeSummarizedExperiment additionally includes the tree and sequences.
#  #' 
#  #' @inherit documentation_default
#  #' 
#  #' @param ...  Not Used.
#  #' 
#  #' @return A SummarizedExperiment or TreeSummarizedExperiment object.
#  #' 
#  #' @export
#  #' @examples
#  #' \dontrun{
#  #'     library(rbiom)
#  #'     
#  #'     print(hmp50)
#  #'     
#  #'     # Requires 'SummarizedExperiment' Bioconductor R package
#  #'     if (nzchar(system.file(package = "SummarizedExperiment"))) {
#  #'       se <- convert_to_SE(hmp50)
#  #'       print(se)
#  #'     }
#  #'     
#  #'     # Requires 'TreeSummarizedExperiment' Bioconductor R package
#  #'     if (nzchar(system.file(package = "TreeSummarizedExperiment"))) {
#  #'       tse <- convert_to_TSE(hmp50)
#  #'       print(tse)
#  #'     }
#  #' }
#  
#  convert_to_SE <- function (biom, ...) {
#    
#    require_package('SummarizedExperiment', 'to use convert_to_SE()')
#    
#    dots <- list(...)
#    biom <- as_rbiom(biom)
#    
#    SummarizedExperiment::SummarizedExperiment(
#      assays  = list('OTU table' = as.matrix(biom$counts)),
#      rowData = biom$taxonomy %>% tibble::column_to_rownames(".otu") %>% as.matrix(),
#      colData = biom$metadata %>% tibble::column_to_rownames(".sample") )
#  }
#  
#  
#  
#  #' @rdname convert_to_SE
#  #' @export
#  
#  convert_to_TSE <- function (biom, ...) {
#    
#    require_package('TreeSummarizedExperiment', 'to use convert_to_TSE()')
#    require_package('Biostrings',               'to use convert_to_TSE()')
#    
#    dots <- list(...)
#    biom <- as_rbiom(biom)
#    
#    TreeSummarizedExperiment::TreeSummarizedExperiment(
#      assays       = list('OTU table' = as.matrix(biom$counts)),
#      rowData      = biom$taxonomy %>% tibble::column_to_rownames(".otu") %>% as.matrix(),
#      colData      = biom$metadata %>% tibble::column_to_rownames(".sample"),
#      rowTree      = biom$tree, 
#      referenceSeq = Biostrings::DNAStringSet(biom$sequences) )
#  }
#  
#  
