#' Subset samples using the BIOM's metadata or taxonomy.
#' 
#' @name subset
#' @param x   A BIOM object, as returned from \link{read.biom}.
#' 
#' @param expr   Logical expression to run on the metadata or taxonomy (not both)
#'        to identify samples or taxa to retain.
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
#'        Applies to non-numeric metadata only.
#'
#' @param fast  Should subsetting the phylogenetic tree and sequences be 
#'        skipped? These slow steps are often not necessary. (Default: FALSE)
#'        
#' @return A \code{BIOM} object.
#' @export
#' @seealso \code{\link{select}}
#' @examples
#'     library(rbiom)
#'     
#'     ex1 <- subset(hmp50, Age > 30)
#'     ex2 <- subset(hmp50, Phylum %in% c("Firmicutes", "Actinobacteria"))
#'     ex3 <- subset(hmp50, `Body Site` %in% c("Saliva", "Stool"))
#'     ex4 <- subset(hmp50, a == b, list(a = as.name("Body Site"), b ="Saliva"))
#'     ex5 <- subset(hmp50, Age < 25 & BMI > 22)
#'
subset.BIOM <- function (x, expr, env = parent.frame(), drop.na = TRUE, refactor = TRUE, fast = FALSE) {
  
  stopifnot(is(x, 'BIOM'))
  
  #--------------------------------------------------------------
  # Don't record sub-calls in the biom's history
  #--------------------------------------------------------------
  biom <- x
  hist <- attr(biom, 'history')
  
  
  #--------------------------------------------------------------
  # Replace variable names with their values
  #--------------------------------------------------------------
  expr <- do.call(substitute, list(substitute(expr), as.list(env)))
  
  
  #--------------------------------------------------------------
  # Metadata or Taxonomy based subsetting?
  #--------------------------------------------------------------
  vars   <- all.vars(expr)
  ranks  <- taxa.ranks(biom)
  mdcols <- colnames(metadata(biom))
  if (length(missing <- setdiff(vars, c(ranks, mdcols))) > 0)
    stop("Unknown object(s): ", paste(missing, collapse = ", "))
  if (!all(vars %in% ranks) && !all(vars %in% mdcols))
    stop("subset expression must be either all metadata or all taxonomy.")
  mode  <- ifelse(all(vars %in% mdcols), "metadata", "taxonomy")
  envir <- if (mode == "metadata") metadata(biom) else as.data.frame(taxonomy(biom))
  remove("vars", "ranks", "mdcols")
  
  
  #--------------------------------------------------------------
  # Evaluate expr to determine which samples to keep
  #--------------------------------------------------------------
  keep <- try(eval(expr, envir = envir), silent=TRUE)
  if (is(keep, "try-error") || is(keep, 'error'))
    stop(simpleError(sprintf("Subset failed: %s", keep)))
  
  
  #--------------------------------------------------------------
  # Taxonomy vs metadata subsetting are done differently
  #--------------------------------------------------------------
  if (mode == "taxonomy") {
    
    biom$taxonomy <- biom$taxonomy[keep,,drop=F]
    biom <- repair(biom, fast=fast)
    
  } else {
    
    #--------------------------------------------------------------
    # Keep only the indicated samples
    #--------------------------------------------------------------
    biom <- select(biom, keep, fast=fast)
    
    
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



