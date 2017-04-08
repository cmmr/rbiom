#' Provides a table of sample diversity values at multiple read depths.
#'
#' A useful function for making rarefaction curve plots.
#'
#' @param biom  A BIOM object, as returned from \link{read.biom}.
#' @param x  A numeric vector giving the depths to rarefy at. If the vector
#'     has only a length of one, the value is taken to mean the number of
#'     different depths to automatically select.
#' @param progressbar  Whether to display a progress bar and status messages
#'     (TRUE/FALSE). Will automatically tie in with \pkg{shiny} if run within a
#'     \pkg{shiny} session.
#' @return A \code{data frame} with column names \bold{Sample}, 
#'     \bold{Depth}, \bold{OTUs}, \bold{Shannon}, \bold{Simpson}, and 
#'     \bold{Chao1}.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     rc <- subsample(biom, x=5)
#'     head(rc)
#' 


subsample <- function (biom, x=10, progressbar=FALSE) {
  
  
  #--------------------------------------------------------------
  # Sanity Checks
  #--------------------------------------------------------------
  
  if (!is(biom, "BIOM"))
    stop(simpleError("subsample function needs a BIOM object."))
  
  if (is.null(x))
    stop(simpleError("subsample: 'x' cannot be NULL."))
  if (!is(x, "numeric"))
    stop(simpleError("subsample: 'x' must be numeric."))
  if (!all(x > 0))
    stop(simpleError("subsample: 'x' be positive."))
  
  
  
  #--------------------------------------------------------------
  # Choose the rarefaction depths to sample at
  #--------------------------------------------------------------
  
  if (length(x) == 1) {
    x            <- floor(unname(x))
    sample_sums  <- slam::col_sums(biom[['counts']])
    rarefyLevels <- unique(floor(10 ** (c(1:x) * (log10(fivenum(sample_sums)[[4]]) / x))))
    
  } else {
    rarefyLevels <- sort(unique(floor(unname(x))))
  }
  
  
  
  #--------------------------------------------------------------
  # Multithread the heavy-lifting
  #--------------------------------------------------------------
  
  pb <- progressBar(progressbar)
  cl <- configCluster(nTasks=max(biom$counts$j), pb, "Subsampling")

  pb$set(0, 'Computing rarefaction curve.')
  
  res <- {
    
    set <- NULL
    
    foreach(set=cl$sets, .combine='rbind', .options.snow=cl$opts) %dopar% {
      
      plyr::ldply(set, function (idx) {
        
        counts <- biom$counts$v[biom$counts$j == idx]
        nReads <- floor(sum(as.numeric(counts)))
        seed   <- set.seed(sum(cumsum(counts > 0) + idx) %% 2147483647)
        
        rLvls  <- rarefyLevels[rarefyLevels <= nReads]
        if (length(rLvls) == 0) return (NULL)
        
        breaks <- c(0, cumsum(counts))
      
        plyr::ldply(rLvls, function(rLvl) {
          
          set.seed(seed)
          otus   <- unname(table(cut(sample.int(sum(counts), rLvl), breaks=breaks)))
          otus   <- otus[otus > 0]
          otus_p <- otus / sum(otus)
          
          data.frame(
            'Sample'  = idx,
            'Depth'   = rLvl,
            'OTUs'    = length(otus),
            'Shannon' = -sum(otus_p * log(otus_p)),
            'Simpson' = 1 - sum(otus_p ** 2),
            'Chao1'   = length(otus) + ( (sum(otus == 1) ** 2) / (2 * sum(otus == 2)) )
          )
        })
      })
    }
  }
  
  res[['Sample']] <- factor(res[['Sample']], 1:cl$nTasks, biom$counts$dimnames[[2]])
  
  
  pb$close()
  
  return (res)
}
