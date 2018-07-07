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
#' Many \code{rbiom} functions support multithreading:
#' 
#' The default behavior of these function is to run on as many cores as are
#' available in the local compute environment. If you wish to limit the number
#' of simultaneous threads, set the \code{rbiom.max.threads} option. For
#' instance:
#' 
#' \preformatted{
#'     options('rbiom.max.threads' = 6)
#' }
#' 
#' Alternatively, you can register a custom parallel backend, which will be
#' detected and used for multithreaded functions. Please note that only the
#' \code{doSNOW} backend currently supports progress tracking. With other
#' backends, you will only be notified of the beginning and ending of computations.
#' 
#' \preformatted{
#'     cl <- parallel::makeCluster(ncores)
#'     doSNOW::registerDoSNOW(cl)
#' }
#'
#' @docType package
#' @name rbiom
NULL
