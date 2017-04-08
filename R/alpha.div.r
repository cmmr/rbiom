#' Estimate the diversity of each sample.
#' 
#' @param biom  A BIOM object, as returned from \link{read.biom}.
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
  
  samples <- biom$counts$dimnames[[2]]
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
          
          otus    <- floor(biom$counts$v[biom$counts$j == idx])
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
