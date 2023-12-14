#' rbiom: Read/Write, Transform, and Summarize BIOM Data
#' 
#' A toolkit for working with Biological Observation Matrix (BIOM) files.
#' Features include reading/writing all BIOM formats, rarefaction, alpha
#' diversity, beta diversity (including UniFrac), summarizing counts by 
#' taxonomic level, and sample subsetting. Standalone functions for 
#' reading, writing, and subsetting phylogenetic trees are also provided. 
#' All CPU intensive operations are encoded in C with multi-thread support.
#' 
#' @section Multithreading:
#' Many rbiom functions support multithreading:
#' 
#' The default behavior of these function is to run on as many cores as are
#' available in the local compute environment. If you wish to limit the number
#' of simultaneous threads, set \code{RcppParallel}'s \code{numThreads} option.
#' For instance:
#' 
#' \preformatted{
#'     RcppParallel::setThreadOptions(numThreads = 4)
#' }
#' 
#'
#' @docType package
#' @keywords internal
#' @aliases rbiom-package
"_PACKAGE"

NULL
