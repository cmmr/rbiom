#' Subset counts so that all samples have the same number of observations.
#'
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'     object, as returned from \link{read.biom}. For matrices, the rows and 
#'     columns are assumed to be the taxa and samples, respectively.
#' @param depth  The number of observations to keep, per sample. If set to
#'     \code{NULL}, a depth will be automatically selected. Samples that have
#'     fewer than this number of observations will be dropped.
#' @param seed  An integer to use for seeding the random number generator. If
#'     you need to create different random rarefactions of the same \code{BIOM}
#'     object, set this seed value to a different number each time.
#' @param progressbar  An object of class \code{Progress}.
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
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     range(slam::col_sums(biom$counts))
#'
#'     biom <- rarefy(biom, depth=1000)
#'     range(slam::col_sums(biom$counts))
#'


rarefy <- function (biom, depth=NULL, seed=0, progressbar=NULL) {
  
  
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
  # Choose the rarefaction depths to sample at
  #--------------------------------------------------------------

  if (is.null(depth)) {

    sample_sums <- slam::col_sums(counts)
    depth       <- (sum(sample_sums) * .1) / length(sample_sums)
    depth       <- min(sample_sums[sample_sums >= depth])
    remove("sample_sums")

  } else {

    if (!is(depth, "numeric"))
      stop(simpleError("In rarefy(), depth must be an integer."))

    if (depth %% 1 != 0)
      stop(simpleError("In rarefy(), depth must be an integer."))
  }
  
  
  pb <- progressBar(progressbar, paste("Rarefying to", depth))
  cl <- configCluster(nTasks=counts$ncol, pb)
  
  res <- {
    
    set <- idx <- NULL
    
    foreach (set=cl$sets, .combine='rbind', .options.snow=cl$opts, .packages='foreach') %dopar% {
      foreach (idx=set, .combine='rbind') %do% {
        
        x <- ceiling(counts$v[counts$j == idx])
        nReads <- sum(x)
        
        if (nReads < depth)
          return (NULL)
        
        sampleID <- ifelse(is.null(colnames(counts)), paste0("Sample", idx), colnames(counts)[[idx]])
        set.seed(sum(utf8ToInt(sampleID) + seed) %% 2147483647)
        
        breaks <- c(0, cumsum(x))
        labels <- counts$i[counts$j == idx]
        retain <- sample.int(nReads, depth)
        retain <- table(cut(retain, breaks=breaks, labels=labels))
        retain <- retain[retain > 0]
        
        matrix(c(as.numeric(names(retain)), rep(idx, length(retain)), unname(retain)), ncol=3)
      }
    }
  }



  #--------------------------------------------------------------
  # Construct a new simple triplet matrix
  #--------------------------------------------------------------

  TaxaIDs   <- counts$dimnames[[1]][sort(unique(res[,1]))]
  SampleIDs <- counts$dimnames[[2]][sort(unique(res[,2]))]

  counts <- slam::simple_triplet_matrix(
    i = as.numeric(factor(res[,1])),
    j = as.numeric(factor(res[,2])),
    v = res[,3],
    dimnames = list(TaxaIDs, SampleIDs)
  )
  
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

  return (biom)

}
