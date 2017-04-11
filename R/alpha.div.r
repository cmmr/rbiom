#' Estimate the diversity of each sample.
#' 
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'     object, as returned from \link{read.biom}. For matrices, the rows and 
#'     columns are assumed to be the taxa and samples, respectively.
#' @param progressbar  Whether to display a progress bar and status messages
#'     (TRUE/FALSE). Will automatically tie in with \pkg{shiny} if run within
#'     a \pkg{shiny} session.
#' @return A numeric matrix of four diversity values for each sample in 
#'     \code{biom}. Row names are the sample names and column names are the 
#'     diversity metrics: \bold{OTUs}, \bold{Shannon}, \bold{Simpson}, and 
#'     \bold{Chao1}.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     ad <- alpha.div(biom)
#'     head(ad)
#'


alpha.div <- function (biom, progressbar=FALSE) {
  
  #--------------------------------------------------------------
  # Get the input into a simple_triplet_matrix
  #--------------------------------------------------------------
  if (is(biom, "simple_triplet_matrix")) { counts <- biom
  } else if (is(biom, "BIOM"))           { counts <- biom$counts
  } else if (is(biom, "matrix"))         { counts <- slam::as.simple_triplet_matrix(biom)
  } else {
    stop(simpleError("biom must be a matrix, simple_triplet_matrix, or BIOM object."))
  }
  
  
  samples <- counts$dimnames[[2]]
  metrics <- c('OTUs', 'Shannon', 'Simpson', 'Chao1')
  
  pb <- progressBar(progressbar = progressbar)
  cl <- configCluster(nTasks = length(samples), pb, "Calculating Alpha Diversity")
  
  matrix(
    ncol     = length(metrics),
    byrow    = TRUE,
    dimnames = list(samples, metrics),
    data     = {
    
      set <- idx <- NULL
      
      foreach (set=cl$sets, .combine='c', .options.snow=cl$opts) %dopar% {
        foreach (idx=set, .combine='c') %do% {
          
          otus    <- floor(counts$v[counts$j == idx])
          otus    <- otus[otus > 0]
          otus_p  <- otus / sum(otus)
        
          nOTUs   <- length(otus)
          Shannon <- -sum(otus_p * log(otus_p))
          Simpson <- 1 - sum(otus_p ** 2)
          Chao1   <- nOTUs + (sum(otus == 1) ** 2) / (2 * sum(otus == 2))
          
          return (c(nOTUs, Shannon, Simpson, Chao1))
        }
      }
    }
  )
}
