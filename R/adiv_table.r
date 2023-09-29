#' Calculate the alpha diversity of each sample.
#' 
#' @inherit adiv_boxplot params
#' 
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM}
#'        object, as returned from [read_biom()]. For matrices, the rows and
#'        columns are assumed to be the taxa and samples, respectively.
#'     
#' @param rarefy  Control how/whether rarefactions are done prior to alpha 
#'        diversity computations. Options are: 
#'        \itemize{
#'          \item{\code{FALSE} - }{
#'            Use each sample's current set of observations without applying
#'            any rarefaction.}
#'          \item{\code{TRUE} - }{
#'            Automatically select and apply a single rarefaction. }
#'          \item{\code{"multi"}, \code{"multi_log"}, \code{"multi_even"} - }{
#'            Automatically select and apply multiple rarefactions.
#'            \code{"multi"} provides \code{"multi_log"} at the low end 
#'            and \code{"multi_even"} at the high end. }
#'          \item{\emph{integer vector} - }{
#'            Rarefy at the specified depth(s). }
#'        }
#'        Default: \code{FALSE}
#'     
#' @param long  Pivot the returned data to long format?
#'        \itemize{
#'          \item{\bold{FALSE} - }{
#'            Each alpha diversity metric always has its own column, named 
#'            ".OTUs", ".Shannon", etc. }
#'          \item{\bold{TRUE} - }{
#'            The name of the metric is in an ".adiv" column, and the 
#'            alpha diversity values are in a ".diversity" column. }
#'        }
#'        Default: \code{FALSE}
#'     
#' @param md  Include metadata in the output data frame? Options are: 
#'        \itemize{
#'          \item{\code{FALSE} - }{ Don't include metadata. }
#'          \item{\code{TRUE} - }{ Include all metadata. }
#'          \item{\emph{character vector} - }{
#'            Include only the specified metadata columns. }
#'        }
#'        Default: \code{FALSE}
#'        
#' @return A data frame of diversity values for each sample in \code{biom}. The 
#'         first column name is \bold{.sample}. Remaining column names are 
#'         dependent on the function's arguments.
#'     
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- sample_select(hmp50, 1:6)
#'     
#'     adiv_table(biom, adiv="all")
#'     adiv_table(biom, long = TRUE, md = TRUE)

adiv_table <- function (biom, rarefy=FALSE, adiv="Shannon", long=FALSE, md=FALSE) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("adiv_table", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
      return (readRDS(cache_file))
  remove("params")
  
  
  
  #________________________________________________________
  # Enable abbreviations of adiv metric names.
  #________________________________________________________
  adiv <- paste0(".", validate_arg(adiv, NULL, 'adiv', 'adiv'))
  
  
  
  
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
  # Define the rarefaction depths to sample at
  #________________________________________________________
  rLvls <- NA
  
  if (identical(rarefy, TRUE)) {
    counts <- sample_rarefy(counts)
    
  } else if (is.numeric(rarefy)) {
    rLvls <- sort(rarefy)
    
  } else if (is.character(rarefy)) {
    upper <- fivenum(slam::col_sums(counts))[[4]]
    
    if (identical(rarefy, "multi")) {
      
      # Log intervals until rLvl/2, then even intervals until rLvl*2
      rLvl  <- rare_suggest(counts)
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
    
    otus <- if (is.na(rLvl)) counts else sample_rarefy(counts, rLvl)
    df   <- rcpp_alpha_div(otus)
    
    
    #________________________________________________________
    # Pivot Longer
    #________________________________________________________
    if (isTRUE(long)) {
      
      mtx <- as.matrix(df[,adiv,drop=F])
      colnames(mtx) %<>% sub(".", "", ., fixed = TRUE)
      
      df <- data.frame(
        stringsAsFactors = FALSE,
        '.sample'    = df[['.sample']][row(mtx)],
        '.depth'     = df[['.depth']][row(mtx)],
        '.adiv'      = colnames(mtx)[col(mtx)],
        '.diversity' = as.numeric(mtx)
      )
        
    } else {
      
      #________________________________________________________
      # The adiv metrics of interest in wide format
      #________________________________________________________
      if (length(rLvls) == 1) {
        df <- df[,c('.sample', adiv),drop=F]
      } else {
        df <- df[,unique(c('.sample', '.depth', adiv)),drop=F]
      }
    }
    
    
    result <- rbind(result, df)
  }
  
  
  
  #________________________________________________________
  # Add rownames and facet attribute when appropriate
  #________________________________________________________
  if (isTRUE(long)) {
    attr(result, 'response') <- ".diversity"
    result[['.adiv']] %<>% factor(levels = sub(".", "", adiv, fixed = TRUE))
    if (length(adiv) > 1) attr(result, 'facet') <- '.adiv'
  }
  
  
  
  #________________________________________________________
  # Add Metadata
  #________________________________________________________
  if (identical(md, TRUE))  md <- colnames(sample_metadata(biom))
  if (identical(md, FALSE)) md <- c()
  for (i in unique(md))
    result[[i]] <- sample_metadata(biom, i)[result[['.sample']]]
  
  
  set_cache_value(cache_file, result)
  
  return (result)
}
