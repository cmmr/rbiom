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
#' @param long  Pivot the returned data to long format?
#'     \describe{
#'       \item{\bold{FALSE}}{ Each metric has its own column. (Default) }
#'       \item{\bold{TRUE}}{ "Sample", "Metric" and "Diversity" are the columns 
#'                       returned. Rows are added to attain all combinations of 
#'                       samples x metrics. }
#'     }
#' @param md  Include metadata in the output data frame? Options are: 
#'     \describe{
#'       \item{\code{FALSE}}{ Don't include metadata. (Default) }
#'       \item{\code{TRUE}}{ Include all metadata. }
#'       \item{\emph{character vector}}{ Include only the specified metadata columns. }
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


alpha.div <- function (biom, rarefy=FALSE, metrics="all", long=FALSE, md=FALSE) {
  
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
    df   <- rcpp_alpha_div(otus)
    
    
    #--------------------------------------------------------------
    # Convert 'all' to actual names of all adiv metrics
    #--------------------------------------------------------------
    if ('all' %in%  metrics)
      metrics <- colnames(df)[-1]
    
    
    #--------------------------------------------------------------
    # Pivot Longer
    #--------------------------------------------------------------
    if (isTRUE(long)) {
      
      if (length(rLvls) > 1 && length(metrics) == 1) {
        
        df <- df[,unique(c('Sample', 'Depth', metrics)),drop=F]
        
      } else if (length(metrics) > 1) {
      
        mtx <- as.matrix(df[,metrics,drop=F])
        df  <- data.frame(
          stringsAsFactors = FALSE,
          Sample    = df[['Sample']][row(mtx)],
          Depth     = df[['Depth']][row(mtx)],
          Metric    = colnames(mtx)[col(mtx)],
          Diversity = as.numeric(mtx)
        )
        
        if (length(rLvls) == 1)
          df <- df[,-2,drop=F]
        
      } else {
        
        df <- df[,c('Sample', metrics),drop=F]
      }
      
    } else {
      
      #--------------------------------------------------------------
      # The metrics of interest in wide format
      #--------------------------------------------------------------
      if (length(rLvls) == 1) {
        df <- df[,c('Sample', metrics),drop=F]
      } else {
        df <- df[,unique(c('Sample', 'Depth', metrics)),drop=F]
      }
    }
    
    
    
    result <- rbind(result, df)
    
  }
  
  
  #--------------------------------------------------------------
  # Add Metadata
  #--------------------------------------------------------------
  if (identical(md, TRUE))  md <- colnames(rbiom::metadata(biom))
  if (identical(md, FALSE)) md <- c()
  for (i in unique(md))
    result[[i]] <- metadata(biom, i)[result[['Sample']]]
  
  
  
  if (length(rLvls) == 1 && !isTRUE(long))
    rownames(result) <- result[['Sample']]
  
  if (isTRUE(long))
    if (length(unique(result[['Metric']])) > 1)
      attr(result, 'facet') <- "Metric"
  
  if (length(metrics) == 1) { attr(result, 'response') <- metrics 
  } else                    { attr(result, 'response') <- "Diversity" }
  
  return (result)
  
}
