#' Subset samples using the BIOM object's metadata
#' 
#' @name subset
#' @param x   A BIOM object, as returned from \link{read.biom}.
#' 
#' @param expr   Logical expression to run on the metadata to identify samples
#'        to retain.
#'        
#' @param env   The environment to search for variables used in \code{expr}.
#'        (Default: \code{parent.frame()})
#'        
#' @param drop.na   When \code{expr} is e.g. \code{Age > 30}, should 
#'        \code{!is.na(Age)} be automatically applied too?
#'        
#' @param refactor  When \code{expr} is e.g. 
#'        \code{`Body Site` \%in\% c("Stool", "Saliva")}, should
#'        \code{`Body Site`} be redefined as 
#'        \code{factor(`Body Site`, levels=c('Stool', 'Saliva'))}? 
#'        
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
#'     ex3 <- subset(biom, a == b, list(a = as.name("Body Site"), b ="Saliva"))
#'     ex4 <- subset(biom, Age < 25 & BMI > 22)
#'
subset.BIOM <- function (x, expr, env = parent.frame(), drop.na = TRUE, refactor = TRUE) {
  
  stopifnot(is(x, 'BIOM'))
  
  #--------------------------------------------------------------
  # Don't record sub-calls in the biom's history
  #--------------------------------------------------------------
  biom <- x
  hist <- attr(biom, 'history')
  
  
  #--------------------------------------------------------------
  # Evaluate expr to determine which samples to keep
  #--------------------------------------------------------------
  expr <- do.call(substitute, list(substitute(expr), as.list(env)))
  keep <- try(eval(expr, envir = x$metadata), silent=TRUE)
  if (is(keep, "try-error") || is(keep, 'error'))
    stop(simpleError(sprintf("Subset failed: %s", keep)))
  
  
  #--------------------------------------------------------------
  # Keep only the indicated samples
  #--------------------------------------------------------------
  biom <- select(biom, keep)
  
  
  #--------------------------------------------------------------
  # Remove NAs in subsetted columns
  #--------------------------------------------------------------
  if (isTRUE(drop.na))
    for (i in all.vars(expr))
      if (i %in% colnames(biom$metadata))
        biom <- select(biom, !is.na(biom$metadata[[i]]))
  
  
  
  #--------------------------------------------------------------
  # Recursively traverse subsetting expression to find all
  # instances of, e.g., `Body Site` %in% c('Stool', 'Saliva')
  #--------------------------------------------------------------
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
  
  
  #--------------------------------------------------------------
  # Convert/update factors according to subset specs
  #--------------------------------------------------------------
  if (isTRUE(refactor)) {
    specs <- in_lists(expr)
    for (key in names(specs)) {
      val <- specs[[key]]
      col <- biom$metadata[[key]]
      if (all(unique(col) %in% val))
        biom$metadata[[key]] <- factor(as.character(col), levels = val)
    }
  }
  
  
  #--------------------------------------------------------------
  # Attach subset() call to provenance tracking
  #--------------------------------------------------------------
  
  cl <- match.call()
  cl <- cl[names(cl) != "env"]
  cl[[1]] <- as.name("subset")
  cl[[2]] <- as.name("biom")
  cl[[3]] <- expr
  for (i in seq_along(names(cl))[-(1:3)]) {
    cl[[i]] <- eval.parent(cl[[i]])
  }
  names(cl)[[2]] <- ""
  names(cl)[[3]] <- ""
  attr(biom, 'history') <- c(hist, paste("biom <-", deparse1(cl)))
  
  return (biom)
}



