#' Estimate the diversity of each sample.
#' 
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM}
#'     object, as returned from \link{read.biom}. For matrices, the rows and
#'     columns are assumed to be the taxa and samples, respectively.
#' @param rarefy  Control how/whether rarefactions are done prior to alpha 
#'     diversity computations. Options are: 
#'     \describe{
#'       \item{\code{FALSE}}{ Use each sample's current set of observations without
#'                       applying any rarefaction. (Default) }
#'       \item{\code{TRUE}}{ Automatically select and apply a single rarefaction. }
#'       \item{\code{"multi"}}{ Automatically select and apply multiple rarefactions.}
#'       \item{\emph{integer vector}}{ Rarefy at the specified depth(s). }
#'     }
#' @param metrics  Character vector of one or more of the following: \code{Depth}, 
#'     \code{OTUs}, \code{Shannon}, \code{Chao1}, \code{Simpson}, \code{InvSimpson}.
#'     Non-ambiguous abbreviations are also accepted. The default, \code{all},
#'     returns all of them.
#' @param md  How to include metadata in the output data frame. Options are: 
#'     \describe{
#'       \item{\code{NULL}}{ Don't include metadata (Default) }
#'       \item{\code{"wide"}}{ Add a separate column for each alpha diversity metric. }
#'       \item{\code{"long"}}{ Add two columns, "Metric" and "Diversity", and add more
#'                       rows to attain all combinations of samples x metric. }
#'     }
#' @return A data frame of four diversity values for each sample in
#'     \code{biom}. The column names are \bold{Sample}, \bold{Depth} and the 
#'     diversity metrics: \bold{OTUs}, \bold{Shannon}, \bold{Chao1}, \bold{Simpson}, 
#'     and \bold{InvSimpson}. The row names are the sample names, except when
#'     multiple rarefactions are done.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     ad <- alpha.div(biom)
#'     head(ad)
#'     
#'     biom <- subset(biom, `Body Site` == "Saliva" & Age < 26)
#'     ad   <- alpha.div(biom, "multi")
#'     boxplot(Shannon ~ Depth, data=ad, xlab="Reads", ylab="Diversity")
#'     


alpha.div <- function (biom, rarefy=FALSE, metrics="all", md=NULL) {
  
  #--------------------------------------------------------------
  # Enable abbreviations of metric names.
  #--------------------------------------------------------------
  if (!identical(metrics, 'all')) {
    metricList <- c("Depth", "OTUs", "Shannon", "Chao1", "Simpson", "InvSimpson")
    metrics    <- metricList[pmatch(tolower(metrics), tolower(metricList))]
    metrics    <- metrics[!is.na(metrics)]
    if (length(metrics) == 0) stop(simpleError("Invalid 'alpha.div(metric=' argument"))
  }
  
  
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
  # Define the rarefaction depths to sample at
  #--------------------------------------------------------------
  if (identical(rarefy, TRUE)) {
    counts <- rbiom::rarefy(counts)
    rLvls <- NA
    
  } else if (identical(rarefy, "multi")) {
    sample_sums <- slam::col_sums(counts)
    rLvls <- floor(10 ** (c(1:10) * (log10(fivenum(sample_sums)[[4]]) / 10)))
    remove("sample_sums")
    
  } else if (is.numeric(rarefy)) {
    rLvls <- sort(rarefy)
    
  } else {
    rLvls <- NA
  }
  
  
  result <- NULL
  
  for (rLvl in rLvls) {
    
    otus <- if (is.na(rLvl)) counts else rbiom::rarefy(counts, rLvl)
    
    df           <- rcpp_alpha_div(otus)
    rownames(df) <- df[['Sample']]
    
    
    if (!identical(metrics, 'all'))
      df <- df[,metrics,drop=F]
    
    
    if (identical(md, 'wide')) {
      wide   <- metadata(biom)
      colIDs <- setdiff(names(df), names(wide))
      rowIDs <- intersect(rownames(wide), rownames(df))
      df     <- df[rowIDs,colIDs,drop=F]
      df     <- cbind(df, wide[rowIDs,,drop=F])
    }
    
    if (identical(md, 'long')) {
      long           <- metadata(biom)[rownames(df),,drop=F]
      rownames(long) <- c()
      df_metric      <- setdiff(names(df), 'Sample')
      df_diversity   <- as.numeric(t(as.matrix(df[,df_metric,drop=F])))
      df             <- data.frame(
        check.names = FALSE,
        'Sample'    = df[['Sample']],
        long,
        'Metric'    = df_metric,
        'Diversity' = df_diversity
      )
    }
    
    
    if (length(rLvls) != 1) rownames(df) <- c()
    
    if (is.null(result)) { result <- df
    } else {               result <- rbind(result, df) }
    
  }
  
  
  return (result)
  
}
