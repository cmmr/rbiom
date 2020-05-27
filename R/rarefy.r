#' Subset counts so that all samples have the same number of observations.
#'
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'     object, as returned from \link{read.biom}. For matrices, the rows and 
#'     columns are assumed to be the taxa and samples, respectively.
#' @param depth  The number of observations to keep, per sample. If set to
#'     \code{NULL}, a depth will be automatically selected. Samples that have
#'     fewer than this number of observations will be dropped. If called on
#'     data with non-integer abundances, values will be re-scaled to integers 
#'     between 1 and \code{depth} such that they sum to \code{depth}.
#' @param seed  An integer to use for seeding the random number generator. If
#'     you need to create different random rarefactions of the same \code{BIOM}
#'     object, set this seed value to a different number each time.
#' @return A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'     object, depending on the input object type. The type of object provided
#'     is the same type that is returned. The retained observations are randomly
#'     selected, based on a seed value derived from the \code{BIOM} object. 
#'     Therefore, rarefying the same biom to the same depth will always produce
#'     the same resultant rarification.
#' @export
#' @examples
#'     library(rbiom)
#'
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     range(slam::col_sums(biom$counts))
#'
#'     biom <- rarefy(biom, depth=1000)
#'     range(slam::col_sums(biom$counts))
#'


rarefy <- function (biom, depth=NULL, seed=0) {
  
  
  #--------------------------------------------------------------
  # Get the input into a simple_triplet_matrix
  #--------------------------------------------------------------
  if (is(biom, "simple_triplet_matrix")) { counts <- biom
  } else if (is(biom, "BIOM"))           { counts <- biom$counts
  } else if (is(biom, "matrix"))         { counts <- slam::as.simple_triplet_matrix(biom)
  } else {
    stop(simpleError("biom must be a matrix, simple_triplet_matrix, or BIOM object."))
  }
  
  
  #--------------------------------------------------------------
  # Choose the rarefaction depth to sample at
  #--------------------------------------------------------------

  if (is.null(depth)) {
    
    if (all(counts$v %% 1 == 0)) {
      
      sample_sums <- slam::col_sums(counts)
      depth       <- (sum(sample_sums) * .1) / length(sample_sums)
      depth       <- min(sample_sums[sample_sums >= depth])
      remove("sample_sums")
      
    } else {
      
      depth <- 1000
      
    }

  } else {

    if (!is(depth, "numeric"))
      stop(simpleError("In rarefy(), depth must be an integer."))

    if (depth %% 1 != 0)
      stop(simpleError("In rarefy(), depth must be an integer."))
    
    if (depth <= 0)
      stop(simpleError("In rarefy(), depth must be positive."))
  }
  
  
  
  
  #--------------------------------------------------------------
  # Integer data. Randomly select observations to keep.
  #--------------------------------------------------------------
  
  if (all(counts$v %% 1 == 0)) {
    
    counts <- rcpp_rarefy(counts, depth, seed)
  
    
  #--------------------------------------------------------------
  # Fractional data. Re-scale between 1 and depth. Total = depth.
  #--------------------------------------------------------------
  
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
  
  
  #--------------------------------------------------------------
  # If biom isn't a BIOM object, return now.
  #--------------------------------------------------------------
  if (is(biom, "simple_triplet_matrix")) return (counts)
  if (is(biom, "matrix")) return (as.matrix(counts))
  


  #--------------------------------------------------------------
  # Drop unobserved taxa and excluded samples from BIOM object
  #--------------------------------------------------------------
  
  biom$counts   <- counts
  biom$taxonomy <- biom$taxonomy[TaxaIDs,,drop=FALSE]
  biom$metadata <- biom$metadata[SampleIDs,,drop=FALSE]

  if (!is.null(biom$phylogeny))
    biom$phylogeny <- rbiom::subtree(biom$phylogeny, TaxaIDs)
  
  if (!is.null(biom$sequences))
    biom$sequences <- biom$sequences[TaxaIDs]
  
  biom$info$shape <- c(length(TaxaIDs), length(SampleIDs))

  
  return (biom)

}
