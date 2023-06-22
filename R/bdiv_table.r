#' Make a data.frame of distances between samples.
#' 
#' @name bdiv_table
#' 
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'        object, as returned from \link{read_biom}. For matrices, the rows and 
#'        columns are assumed to be the taxa and samples, respectively.
#'        
#' @param method  The distance algorithm to use. Options are:
#'        \code{"bray-curtis"}, \code{"manhattan"}, \code{"euclidean"}, 
#'        \code{"jaccard"}, and \code{"unifrac"}. Non-ambiguous abbreviations 
#'        of the method names are also accepted. A phylogenetic tree must be 
#'        present in \code{biom} or explicitly provided via \code{tree=} to use
#'        the UniFrac methods. (Default: \code{"bray-curtis"})
#'     
#' @param weighted  Take relative abundances into account. When 
#'        \code{weighted=FALSE}, only presence/absence is considered.
#'         (Default: \code{TRUE})
#'        
#' @param tree  A \code{phylo} object representing the phylogenetic
#'        relationships of the taxa in \code{biom}. Will be taken from the tree
#'        embedded in the \code{biom} object if not explicitly specified. Only
#'        required for computing UniFrac distance matrices.
#'         (Default: \code{NULL})
#'        
#' @param md  Include metadata in the output data frame? Options are: 
#'        \describe{
#'          \item{\code{FALSE}}{ Don't include metadata. (Default) }
#'          \item{\code{TRUE}}{ Include all metadata. }
#'          \item{\emph{character vector}}{ Include only the specified metadata 
#'          columns. Column names can be prefixed with \bold{==} or \bold{!=} to 
#'          indicate that only within or between groupings, respectively, are to 
#'          be kept. See examples below. }
#'        }
#'        
#' @param stat.by  Specify a categorical metadata column name to compute adonis 
#'        statistics and return results as attributes named 'stats_raw' and 
#'        'stats_tbl'. (Default: \code{NULL})
#'        
#' @param perms  Number of random permutations to use in adonis calcuation.
#'        (Default: \code{999})
#'        
#' @param seed  Random seed for adonis permutations. (Default: \code{0})
#'        
#' @return A data.frame with first three columns named ".sample1", ".sample2", 
#'         and ".distance".
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     # Return in long format with metadata
#'     biom <- select(hmp50, 18:21)
#'     bdiv_table(biom, 'unifrac', md = c("Body Site", "Sex"))
#'     
#'     # Only look at distances among the stool sample
#'     bdiv_table(biom, 'unifrac', md = c("==Body Site", "Sex"))
#'     
#'     # Or between males and females
#'     bdiv_table(biom, 'unifrac', md = c("Body Site", "!=Sex"))
#'

bdiv_table <- function (
    biom, method="bray-curtis", weighted=TRUE, tree=NULL, md=FALSE, stat.by=NULL, 
    seed=0, perms=999) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("bdiv_table", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  history <- attr(biom, 'history')
  history %<>% c(sprintf("bdiv_table(%s)", as.args(params, fun = bdiv_table)))
  
  
  #________________________________________________________
  # Compute the distance matrix
  #________________________________________________________
  dm <- bdiv_distmat(biom, method, weighted, tree, stat.by, seed, perms)
  
  
  #________________________________________________________
  # Convert to long form
  #________________________________________________________
  df <- as.matrix(dm)
  df <- data.frame(
    stringsAsFactors = FALSE,
    '.sample1'       = rownames(df)[row(df)],
    '.sample2'       = colnames(df)[col(df)],
    '.distance'      = as.numeric(df)
  )
  df <- subset(df, .sample1 < .sample2)
  
  
  #________________________________________________________
  # Add metadata columns
  #________________________________________________________
  if (!isFALSE(md)) {
    
    if (isTRUE(md))        md <- names(metadata(biom))
    if (!is.character(md)) md <- names(metadata(biom))[md]
    
    for (i in which(!duplicated(md))) {
      
      col <- md[[i]]
      op  <- attr(md, 'op', exact = TRUE)[[i]]
      
      # Convert '==' or '!=' prefix to an attribute
      #________________________________________________________
      if (isTRUE(substr(col, 1, 2) %in% c("==", "!="))) {
        op  <- substr(col, 1, 2)
        col <- substr(col, 3, nchar(col))
        attr(col, 'op') <- op
      }
      
      map <- metadata(biom, col)
      
      # Limit to only within or between comparisons.
      #________________________________________________________
      if (!is_null(attr(col, 'op', exact = TRUE))) {
        op <- attr(col, 'op', exact = TRUE)
        df <- df[get(op)(map[df$.sample1], map[df$.sample2]),,drop=F]
      }
      
      v1 <- as.character(map[df$.sample1])
      v2 <- as.character(map[df$.sample2])
      
      
      # Change "Male vs Female" to "Female vs Male" (alphabetical).
      #________________________________________________________
      df[[col]] <- paste(
        ifelse(v1 < v2, v1, v2), 
        "vs", 
        ifelse(v1 < v2, v2, v1) )
      
      
      # Change "Male vs Male" to "Male".
      #________________________________________________________
      df[[col]] <- ifelse(v1 == v2, v1, df[[col]])
      
      
      # Keep factors as factors when possible
      #________________________________________________________
      if (is.factor(map)) {
        if (identical(op, "==")) {
          df[[col]] <- factor(df[[col]], levels = levels(map))
        } else {
          df[[col]] <- factor(df[[col]])
        }
      }
      
    }
  }
  
  if (!is_null(stat.by)) {
    attr(df, 'stats_raw') <- attr(dm, 'stats_raw', exact = TRUE)
    attr(df, 'stats_tbl') <- attr(dm, 'stats_tbl', exact = TRUE)
  }
  
  
  attr(df, 'response') <- ".distance"
  attr(df, 'history')  <- history
  
  
  set_cache_value(cache_file, df)
  return (df)
}



