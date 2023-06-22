#' Estimate the diversity of each sample.
#' 
#' @name adiv_table
#' 
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM}
#'        object, as returned from \link{read_biom}. For matrices, the rows and
#'        columns are assumed to be the taxa and samples, respectively.
#'     
#' @param rarefy  Control how/whether rarefactions are done prior to alpha 
#'        diversity computations. Options are: 
#'        \describe{
#'        
#'          \item{\code{FALSE}}{
#'            Use each sample's current set of observations without applying
#'            any rarefaction. (Default) }
#'            
#'          \item{\code{TRUE}}{
#'            Automatically select and apply a single rarefaction. }
#'            
#'          \item{\code{"multi"}, \code{"multi_log"}, \code{"multi_even"}}{
#'            Automatically select and apply multiple rarefactions.
#'            \code{"multi"} provides \code{"multi_log"} at the low end 
#'            and \code{"multi_even"} at the high end. }
#'            
#'          \item{\emph{integer vector}}{
#'            Rarefy at the specified depth(s). }
#'        }
#'     
#' @param metrics  Character vector of one or more of the following: 
#'        \code{depth}, \code{otus}, \code{shannon}, \code{chao1}, 
#'        \code{simpson}, \code{inv_simp}. Non-ambiguous abbreviations are 
#'        also accepted. A ".depth" column will always be included when alpha 
#'        diversity is calculated at different rarefaction levels. 
#'        The default, \code{all}, returns all metrics.
#'     
#' @param long  Pivot the returned data to long format?
#'        \describe{
#'          \item{\bold{FALSE} (Default)}{
#'            Each metric always has its own column, named ".depth", ".otus", 
#'            ".shannon", etc. }
#'          \item{\bold{TRUE}}{
#'            The name of the metric is in a ".metric" column, and the 
#'            alpha diversity values are in a ".diversity" column. }
#'        }
#'     
#' @param md  Include metadata in the output data frame? Options are: 
#'     \describe{
#'       \item{\code{FALSE} (Default)}{ Don't include metadata. }
#'       \item{\code{TRUE}}{ Include all metadata. }
#'       \item{\emph{character vector}}{
#'         Include only the specified metadata columns. }
#'     }
#'        
#' @return A data frame of diversity values for each sample in \code{biom}. The 
#'         first column name is \bold{.sample}. Remaining column names are 
#'         dependent on the function's arguments.
#'     
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     ad <- adiv_table(hmp50)
#'     head(ad)
#'     
#'     biom <- subset(hmp50, `Body Site` == "Saliva" & Age < 26)
#'     ad   <- adiv_table(biom, "multi")
#'     boxplot(.shannon ~ .depth, data=ad, xlab="Reads", ylab="Diversity")
#'     


adiv_table <- function (biom, rarefy=FALSE, metrics="all", long=FALSE, md=FALSE) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("adiv_table", params)
  if (!is.null(cache_value <- get_cache_value(cache_file)))
    return (cache_value)
  
  
  #________________________________________________________
  # Enable abbreviations of metric names.
  #________________________________________________________
  if (!identical(metrics, 'all')) {
    metricList <- c("depth", "otus", "shannon", "chao1", "simpson", "inv_simp")
    metrics    <- metricList[pmatch(tolower(metrics), tolower(metricList))]
    metrics    <- metrics[!is.na(metrics)]
    if (length(metrics) == 0) stop(simpleError("Invalid 'adiv_table(metric=' argument"))
  }
  
  
  #________________________________________________________
  # Convert 'all' to actual names of all adiv metrics
  #________________________________________________________
  if ('all' %in% metrics)
    metrics <- c("depth", "otus", "shannon", "chao1", "simpson", "inv_simp")
  
  metrics <- paste0(".", metrics)
  
  
  
  #________________________________________________________
  # Get the input into a simple_triplet_matrix
  #________________________________________________________
  if (is(biom, "simple_triplet_matrix")) { counts <- biom
  } else if (is(biom, "BIOM"))           { counts <- biom$counts
  } else if (is(biom, "matrix"))         { counts <- slam::as.simple_triplet_matrix(biom)
  } else {
    stop(simpleError("biom must be a matrix, simple_triplet_matrix, or BIOM object."))
  }
  
  
  #________________________________________________________
  # Define the rarefaction depths to sample at
  #________________________________________________________
  rLvls <- NA
  
  if (identical(rarefy, TRUE)) {
    counts <- rarefy(counts)
    
  } else if (is.numeric(rarefy)) {
    rLvls <- sort(rarefy)
    
  } else if (is.character(rarefy)) {
    upper <- fivenum(slam::col_sums(counts))[[4]]
    
    if (identical(rarefy, "multi")) {
      
      # Log intervals until rLvl/2, then even intervals until rLvl*2
      rLvl  <- default_rarefaction_depth(counts)
      rLvls <- 10 ** (c(1,3,5,7,8,9) * (log10(rLvl / 2) / 10))
      rLvls <- c(rLvls, seq(from = rLvl / 2, to = rLvl * 2, length.out = 5))
      rLvls <- floor(rLvls)
      remove("rLvl")
      
    } else if (identical(rarefy, "multi_log")) {
      rLvls <- floor(10 ** (c(1:10) * (log10(upper) / 10)))
      
    } else if (identical(rarefy, "multi_even")) {
      rLvls <- floor(seq(from = 5, to = upper, length.out = 10))
    }
    remove("upper")
  }
  
  
  result <- NULL
  
  for (rLvl in rLvls) {
    
    otus <- if (is.na(rLvl)) counts else rarefy(counts, rLvl)
    df   <- rcpp_alpha_div(otus)
    
    
    #________________________________________________________
    # Pivot Longer
    #________________________________________________________
    if (isTRUE(long)) {
      
      mtx <- as.matrix(df[,metrics,drop=F])
      colnames(mtx) <- sub("^\\.", "", metrics)
      
      df <- data.frame(
        stringsAsFactors = FALSE,
        '.sample'    = df[['.sample']][row(mtx)],
        '.depth'     = df[['.depth']][row(mtx)],
        '.metric'    = colnames(mtx)[col(mtx)],
        '.diversity' = as.numeric(mtx)
      )
        
    } else {
      
      #________________________________________________________
      # The metrics of interest in wide format
      #________________________________________________________
      if (length(rLvls) == 1) {
        df <- df[,c('.sample', metrics),drop=F]
      } else {
        df <- df[,unique(c('.sample', '.depth', metrics)),drop=F]
      }
    }
    
    
    result <- rbind(result, df)
  }
  
  
  
  #________________________________________________________
  # Add rownames and facet attribute when appropriate
  #________________________________________________________
  if (isTRUE(long)) {
    attr(result, 'response') <- ".diversity"
    result[['.metric']] %<>% factor(levels = sub("^\\.", "", metrics))
    if (length(metrics) > 1) attr(result, 'facet') <- '.metric'
  }
  
  
  
  #________________________________________________________
  # Add Metadata
  #________________________________________________________
  if (identical(md, TRUE))  md <- colnames(metadata(biom))
  if (identical(md, FALSE)) md <- c()
  for (i in unique(md))
    result[[i]] <- metadata(biom, i)[result[['.sample']]]
  
  
  set_cache_value(cache_file, result)
  
  return (result)
}
