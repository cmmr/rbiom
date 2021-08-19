#' Subset samples using the BIOM's metadata or taxonomy.
#' 
#' @md
#' @section Taxonomic abundance filtering:
#' 
#' For taxonomic subsetting, several functions are added or overridden to
#' behave as expected within the subsetting expression. They are:
#' 
#' \code{mean()}, \code{median()}, \code{min()}, \code{max()}, \code{n()}, 
#' \code{count()}, and \code{pct()}.
#' 
#' Therefore you can write 
#' \code{subset(hmp50, mean(Genus) >= 0.1)} and the returned BIOM object will 
#' contain only the genera that average at least 10% relative abundance across
#' all the samples.
#' 
#' If you want only orders that are present in three or more samples, you can 
#' do: \code{subset(hmp50, count(Order) >= 3)}. To require presence in 25% of 
#' samples, you'd use: \code{subset(hmp50, pct(Order) >= 0.25)}.
#' 
#' Both \code{count()} and \code{pct()} have default arguments of 
#' \code{gt=0, le=1, ge=NULL, lt=NULL}, which can be overridden to find, e.g., 
#' which genera comprise at least 2% of the community in 10% or more of the 
#' samples: \code{subset(hmp50, pct(Genus, ge=0.02) >= 0.10)}. 
#' 
#' \emph{gt = greater than, ge = greater than or equal to. lt/le similarly with 
#' 'less than'.}
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
#'   \dontrun{
#'     library(rbiom)
#'     
#'     subset(hmp50, `Body Site` %in% c("Saliva", "Stool"))
#'     subset(hmp50, Age < 25 & BMI > 22)
#'     subset(hmp50, Phylum %in% c("Firmicutes", "Actinobacteria"))
#'     subset(hmp50, mean(Genus) > 0.1)
#'     subset(hmp50, a == b, list(a = as.name("Body Site"), b ="Saliva"))
#'  }
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
  remove("vars", "ranks", "mdcols")
  
  
  #--------------------------------------------------------------
  # Taxonomy vs metadata subsetting are done differently
  #--------------------------------------------------------------
  if (mode == "taxonomy") {
    envir <- c(
      as.list(as.data.frame(taxonomy(biom))),
      list(OTU = taxa.names(biom)) )
    
    
    #--------------------------------------------------------------
    # Formula functions that operate on relative abundances
    #--------------------------------------------------------------
    RelAbFunc <- function (fun, x, ...) {
      assign(as.character(substitute(fun)),fun)
      rank <- as.character(substitute(x))
      raw  <- taxa.rollup(biom, rank)
      rel  <- raw / rowSums(raw)
      res  <- apply(rel, 2L, fun, ...)
      res[x]
    }
    envir[['mean']]     <- function (...) RelAbFunc(mean,   ...)
    envir[['median']]   <- function (...) RelAbFunc(median, ...)
    envir[['min']]      <- function (...) RelAbFunc(min,    ...)
    envir[['max']]      <- function (...) RelAbFunc(max,    ...)
    envir[['n']]        <- function (...) RelAbFunc(length, ...)
    
    count <- function (x, gt=0, le=1, ge=NULL, lt=NULL, pct=FALSE) {
      lower_pass <- if (is.null(ge)) x >  gt else x >= ge
      upper_pass <- if (is.null(lt)) x <= le else x >  lt
      sum(lower_pass & upper_pass) / ifelse(isTRUE(pct), length(x), 1)
    }
    envir[['count']] <- function (...) RelAbFunc(count, ...)
    envir[['pct']]   <- function (...) RelAbFunc(count, ..., pct=TRUE)
    
    
    #--------------------------------------------------------------
    # Evaluate expr to determine which taxa to keep
    #--------------------------------------------------------------
    keep <- try(eval(expr, envir = envir), silent=TRUE)
    if (is(keep, "try-error") || is(keep, 'error'))
      stop(simpleError(sprintf("Subset failed: %s", keep)))
    
    
    biom$taxonomy <- biom$taxonomy[keep,,drop=F]
    biom <- repair(biom, fast=fast)
    
  } else {
    envir <- metadata(biom)
    
    #--------------------------------------------------------------
    # Evaluate expr to determine which samples to keep
    #--------------------------------------------------------------
    keep <- try(eval(expr, envir = envir), silent=TRUE)
    if (is(keep, "try-error") || is(keep, 'error'))
      stop(simpleError(sprintf("Subset failed: %s", keep)))
    
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



