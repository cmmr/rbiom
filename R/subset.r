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
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     ex1 <- subset(biom, Age > 30)
#'     ex2 <- subset(biom, `Body Site` %in% c("Saliva", "Stool"))
#'     ex3 <- subset(biom, Age < 25 & BMI > 22)
#'
subset.BIOM <- function (x, ...) {
  
  res <- try(eval(substitute(...), x$metadata), silent=TRUE)
  if (is(res, "try-error"))
    stop(simpleError(sprintf("Subset failed: %s", res)))
  
  select(x, res)
}

