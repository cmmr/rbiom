#' Subset samples using the BIOM object's metadata
#' 
#' @name subset
#' @param x         A BIOM object, as returned from \link{read.biom}.
#' @param ...       Test to run on the metadata to identify samples to retain.
#' @param drop.na   When \code{...} is e.g. \code{Age > 30}, should \code{!is.na(Age)} be automatically applied too?
#' @param refactor  When \code{...} is e.g. \code{`Body Site` \%in\% c("Stool", "Saliva")}, should
#'                  \code{`Body Site`} be redefined as \code{factor(`Body Site`, levels=c('Stool', 'Saliva'))}? 
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
subset.BIOM <- function (x, ..., drop.na = TRUE, refactor = TRUE) {
  
  stopifnot(is(x, 'BIOM'))
  
  expr = substitute(...)
  
  res <- try(eval(expr = expr, envir = x$metadata), silent=TRUE)
  if (is(res, "try-error") || is(res, 'error'))
    stop(simpleError(sprintf("Subset failed: %s", res)))
  
  x <- select(x, res)
  
  # Remove NAs in subsetted columns
  if (isTRUE(drop.na))
    for (i in all.vars(expr))
      x <- select(x, !is.na(x$metadata[[i]]))
  
  # Convert/update factors according to subset specs
  if (isTRUE(refactor)) {
    specs <- in_lists(expr)
    for (key in names(specs)) {
      val <- specs[[key]]
      col <- x$metadata[[key]]
      if (all(unique(col) %in% val))
        x$metadata[[key]] <- factor(as.character(col), levels = val)
    }
  }
  
  return (x)
}



# Recursively traverse subsetting expression to find all instances 
# of, e.g., `Body Site` %in% c('Stool', 'Saliva')

in_lists <- function (expr) {
  result <- list()
    
  if (isTRUE(expr[[1]] == "%in%")) {
    
    key <- as.character(expr[[2]])
    val <- expr[[3]]
    
    if (is(val, "call") && val[[1]] == "c")
        val <- sapply(val[-1], unlist)
    
    result[[key]] <- as.character(val)
    
  } else {
    
    for (i in seq_len(length(expr))) {
      if (length(expr[[i]]) > 1) {
        result <- c(result, in_lists(expr[[i]]))
      }
    }
  }
  
  result <- result[!duplicated(names(result))]
  
  return (result)
}



