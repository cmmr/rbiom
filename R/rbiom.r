#' rbiom: A package for computating the notorious bar statistic.
#'
#' The foo package provides three categories of important functions:
#' foo, bar and baz.
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
#'     options('rbiom.max.threads') <- 6
#' }
#' 
#' Alternatively, you can register a custom parallel backend, which will be
#' detected and used for multithreaded functions. Please note that only the
#' \code{doSNOW} backend currently supports progress tracking. With other
#' backends, you will only be notified of beginning and end of computations.
#' 
#' \preformatted{
#'     cl <- parallel::makeCluster(ncores)
#'     doSNOW::registerDoSNOW(cl)
#' }
#'
#' @docType package
#' @name rbiom
NULL
