#' Subset samples using the BIOM object's metadata
#' 
#' @name subset
#' @param x  A BIOM object, as returned from \link{read.biom}.
#' @param ...  Test to run on the metadata to identify samples to retain.
#' @return A \code{BIOM} object.
#' @export
#' @seealso \code{\link{select}}
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     ex1 <- subset(biom, Age > 30)
#'     ex2 <- subset(biom, `Body Site` %in% c("Saliva", "Stool"))
#'     ex3 <- subset(biom, Age < 25 & BMI > 22)
#'
subset.BIOM <- function (x, ...) {
  
  stopifnot(is(x, 'BIOM'))
  
  res <- try(eval(expr = substitute(...), envir = x$metadata), silent=TRUE)
  if (is(res, "try-error") || is(res, 'error'))
    stop(simpleError(sprintf("Subset failed: %s", res)))
  
  select(x, res)
}

