#' Subset counts so that all samples have the same number of observations.
#'
#' @param biom   A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'        object, as returned from [read_biom()]. For matrices, the rows and
#'        columns are assumed to be the taxa and samples, respectively.
#'     
#' @param depth   The number of observations to keep, per sample. If set to
#'        \code{NULL}, a depth will be automatically selected. Samples that 
#'        have fewer than this number of observations will be dropped. If 
#'        called on data with non-integer abundances, values will be re-scaled 
#'        to integers between 1 and \code{depth} such that they sum to 
#'        \code{depth}.
#'     
#' @param seed   An integer to use for seeding the random number generator. If
#'        you need to create different random rarefactions of the same 
#'        \code{BIOM} object, set this seed value to a different number each 
#'        time.
#'     
#' @return A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'         object, depending on the input object type. The type of object 
#'         provided is the same type that is returned. The retained 
#'         observations are randomly selected, based on a seed value derived 
#'         from the \code{BIOM} object. Therefore, rarefying the same biom to 
#'         the same depth will always produce the same resultant rarefaction.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'
#'     range(slam::col_sums(hmp50$counts))
#'
#'     biom <- sample_rarefy(hmp50, depth=1000)
#'     range(slam::col_sums(biom$counts))
#'

sample_rarefy <- function (biom, depth=NULL, seed=0) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("sample_rarefy", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Get the input into a simple_triplet_matrix
  #________________________________________________________
  if (is(biom, "simple_triplet_matrix")) { counts <- biom
  } else if (is(biom, "BIOM"))           { counts <- biom$counts
  } else if (is(biom, "matrix"))         { counts <- slam::as.simple_triplet_matrix(biom)
  } else {
    stop("biom must be a matrix, simple_triplet_matrix, or BIOM object.")
  }
  
  
  #________________________________________________________
  # Remove rows/cols with zero observations
  #________________________________________________________
  counts <- counts[slam::row_sums(counts) > 0, slam::col_sums(counts) > 0]
  
  
  #________________________________________________________
  # Choose the rarefaction depth to sample at
  #________________________________________________________

  if (is_null(depth))
    depth <- rare_suggest(counts)
  
  stopifnot(is_scalar_integerish(depth) && !is_na(depth))
  stopifnot(depth > 0)
  
  
  
  
  #________________________________________________________
  # Integer data. Randomly select observations to keep.
  #________________________________________________________
  
  if (all(counts$v %% 1 == 0)) {
    
    counts <- rcpp_rarefy(counts, depth, seed)
  
    
  #________________________________________________________
  # Fractional data. Re-scale between 1 and depth. Total = depth.
  #________________________________________________________
  
  } else {
    
    counts <- as.matrix(counts)
    
    # Does this, but then nudges up and down until == depth:
    # round(t((t(counts) / colSums(counts))) * depth)
    
    counts <- apply(counts, 2L, function (x) { 
      x        <- x / sum(x)
      y        <- round(x * depth)
      maxIters <- length(x)
      
      while (sum(y) < depth && maxIters > 0) {
        i        <- which.min( ((y + 1) / depth) - x )
        y[i]     <- y[i] + 1
        maxIters <- maxIters - 1
      }
      while (sum(y) > depth && maxIters > 0) {
        i        <- which.min( x - ((y - 1) / depth) )
        y[i]     <- y[i] - 1
        maxIters <- maxIters - 1
      }
      
      return (y)
    })
    
    counts    <- slam::as.simple_triplet_matrix(counts)
    
  }
  
  TaxaIDs   <- rownames(counts)
  SampleIDs <- colnames(counts)
  
  
  #________________________________________________________
  # If biom isn't a BIOM object, return now.
  #________________________________________________________
  if (is(biom, "simple_triplet_matrix")) return (counts)
  if (is(biom, "matrix")) return (as.matrix(counts))
  


  #________________________________________________________
  # Drop unobserved taxa and excluded samples from BIOM object
  #________________________________________________________
  
  biom$counts   <- counts
  biom$taxonomy <- biom$taxonomy[TaxaIDs,,drop=FALSE]
  biom$metadata <- biom$metadata[SampleIDs,,drop=FALSE]

  if (!is_null(biom$phylogeny))
    biom$phylogeny <- tree_subset(biom$phylogeny, TaxaIDs)
  
  if (!is_null(biom$sequences))
    biom$sequences <- biom$sequences[TaxaIDs]
  
  biom$info$shape <- c(length(TaxaIDs), length(SampleIDs))
  
  
  #________________________________________________________
  # Drop missing factor levels from metadata
  #________________________________________________________
  for (i in colnames(biom$metadata))
    if (is.factor(biom$metadata[[i]]))
      biom$metadata[[i]] <- factor(
        x      = as.character(biom$metadata[[i]]), 
        levels = intersect(
          levels(biom$metadata[[i]]),
          unique(as.character(biom$metadata[[i]])) ))
  
  
  #________________________________________________________
  # Attach sample_rarefy() call to provenance tracking
  #________________________________________________________
  cl      <- match.call()
  cl[[1]] <- as.name("sample_rarefy")
  cl[[2]] <- as.name("biom")
  for (i in seq_along(cl)[-(1:2)]) {
    cl[i] <- list(eval.parent(cl[[i]]))
  }
  names(cl)[[2]] <- ""
  
  attr(biom, 'history') <- paste0(collapse = "\n", c(
    attr(biom, 'history', exact = TRUE),
    sprintf("biom <- %s", deparse1(cl)) ))
  
  
  
  #________________________________________________________
  # Record rarefaction level
  #________________________________________________________
  attr(biom, 'rarefaction') <- depth

  
  set_cache_value(cache_file, biom)
  return (biom)
}


#' Suggest a 'good' rarefaction depth.
#'
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'        object, as returned from [read_biom()]. For matrices, the rows and
#'        columns are assumed to be the taxa and samples, respectively.
#'     
#' @return An integer.
#' 
#' @export

rare_suggest <- function (biom) {
  
  
  #________________________________________________________
  # Get the input into a simple_triplet_matrix
  #________________________________________________________
  if (is(biom, "simple_triplet_matrix")) { counts <- biom
  } else if (is(biom, "BIOM"))           { counts <- biom$counts
  } else if (is(biom, "matrix"))         { counts <- slam::as.simple_triplet_matrix(biom)
  } else {
    stop("biom must be a matrix, simple_triplet_matrix, or BIOM object.")
  }
  
  
  if (all(counts$v %% 1 == 0)) {
    
    sample_sums <- slam::col_sums(counts)
    depth       <- (sum(sample_sums) * .1) / length(sample_sums)
    depth       <- min(sample_sums[sample_sums >= depth])
    remove("sample_sums")
    
  } else {
    
    depth <- 10000
    
  }
  
  return (as.integer(depth))
}

