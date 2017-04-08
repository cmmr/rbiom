#' Subset counts so that all samples have the same number of observations.
#'
#' @param biom  A BIOM object, as returned from \link{read.biom}.
#' @param depth  The number of observations to keep, per sample. If set to
#'     \code{NULL}, a depth will be automatically selected. Samples that have
#'     fewer than this number of observations will be dropped.
#' @param seed  An integer to use for seeding the random number generator. If
#'     you need to create different random rarefactions of the same \code{BIOM}
#'     object, set this seed value to a different number each time.
#' @param progressbar  Whether to display a progress bar and status messages
#'     (logical). Will automatically tie in with \pkg{shiny} if run within a
#'     \pkg{shiny} session.
#' @return A \code{BIOM} class object. See \link{read.biom} return value.
#'     The retained observations are randomly selected, based on a seed value
#'     derived from the \code{BIOM} object. Therefore, rarefying the same biom
#'     to the same depth will always produce the same resultant rarified biom.
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


rarefy <- function (biom, depth=NULL, seed=0, progressbar=FALSE) {


  #--------------------------------------------------------------
  # Choose the rarefaction depths to sample at
  #--------------------------------------------------------------

  if (is.null(depth)) {

    sample_sums <- slam::col_sums(biom[['counts']])
    depth       <- (sum(sample_sums) * .1) / length(sample_sums)
    depth       <- min(sample_sums[sample_sums >= depth])

  } else {

    if (!is(depth, "numeric"))
      stop(simpleError("In rarefy(), depth must be an integer."))

    if (depth %% 1 != 0)
      stop(simpleError("In rarefy(), depth must be an integer."))
  }
  
  
  pb <- progressBar(progressbar)
  cl <- configCluster(nTasks=biom$counts$ncol, pb, "Rarefying")
  
  res <- {
    
    set <- idx <- NULL
    
    foreach (set=cl$sets, .combine='rbind', .options.snow=cl$opts) %dopar% {
      foreach (idx=set, .combine='rbind') %do% {
        
        counts <- floor(biom$counts$v[biom$counts$j == idx])
        nReads <- sum(counts)
        
        if (nReads < depth)
          return (NULL)
        
        set.seed(sum(cumsum(counts > 0) + seed + idx) %% 2147483647)
        
        breaks <- c(0, cumsum(counts))
        labels <- biom$counts$i[biom$counts$j == idx]
        retain <- sample.int(nReads, depth)
        retain <- table(cut(retain, breaks=breaks, labels=labels))
        retain <- retain[retain > 0]
        
        matrix(c(as.numeric(names(retain)), rep(idx, length(retain)), unname(retain)), ncol=3)
      }
    }
  }



  #--------------------------------------------------------------
  # Drop unobserved taxa and excluded samples
  #--------------------------------------------------------------

  TaxaIDs       <- biom$counts$dimnames[[1]][sort(unique(res[,1]))]
  SampleIDs     <- biom$counts$dimnames[[2]][sort(unique(res[,2]))]

  biom$taxonomy <- biom$taxonomy[TaxaIDs,,drop=FALSE]
  biom$metadata <- biom$metadata[SampleIDs,,drop=FALSE]

  if (!is.null(biom$phylogeny)) {
    dropped        <- setdiff(biom$phylogeny$tip.label, TaxaIDs)
    biom$phylogeny <- ape::drop.tip(biom$phylogeny, dropped)
  }



  #--------------------------------------------------------------
  # Construct a new simple triplet matrix
  #--------------------------------------------------------------

  biom$counts <- slam::simple_triplet_matrix(
    i = as.numeric(factor(res[,1])),
    j = as.numeric(factor(res[,2])),
    v = res[,3],
    dimnames = list(TaxaIDs, SampleIDs)
  )


  pb$close()

  return (biom)

}
