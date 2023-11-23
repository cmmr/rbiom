

#' Taxa abundances per sample, at the specified taxonomic rank.
#' 
#' @inherit documentation_default
#' 
#' @family taxa_abundance
#' 
#' @return A numeric matrix with taxa as rows, and samples as columns.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     taxa_ranks(hmp50)
#'     
#'     phyla <- taxa_matrix(hmp50, 'Phylum')
#'     phyla[1:4,1:6]
#'

taxa_matrix <- function (
    biom, rank = -1, taxa = NULL, lineage = FALSE, 
    sparse = FALSE, unc = "singly", other = FALSE ) {
  
  validate_biom(clone = FALSE)
  params <- eval_envir(environment())
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  validate_taxa(null_ok = TRUE)
  validate_rank()
  validate_bool("sparse")
  
  
  #________________________________________________________
  # Fetch taxonomy map with ambiguous names corrected.
  #________________________________________________________
  map <- otu_taxonomy(
    biom    = biom, 
    rank    = rank, 
    unc     = unc, 
    lineage = lineage )
  
  
  
  #________________________________________________________
  # Compute abundance matrix for this rank.
  #________________________________________________________
  mtx <- tryCatch(
    error = function (e) stop("Unable to group by taxonomic level: ", e),
    expr  = local({
      
      counts <- biom[['counts']]
      mtx <- slam::rollup(counts[names(map),], 1L, map, sum)
      mtx <- mtx[order(tolower(rownames(mtx))), colnames(counts), drop=FALSE]
      
      stopifnot(is(mtx, "simple_triplet_matrix"))
      
      return (mtx)
    }))
  
  
  
  #________________________________________________________
  # Only retain taxa of interest (by abundance or name).
  #________________________________________________________
  if (!is_null(taxa)) {
    
    if (is.numeric(taxa) && length(taxa) == 1) {
      rel <- sort(slam::row_means(t(t(mtx) / slam::col_sums(mtx))), decreasing = TRUE)
      if (taxa >= 1) { taxa <- head(names(rel), taxa) 
      } else         { taxa <- names(rel)[rel >= taxa] }
      
    } else if (is.character(taxa)) {
      taxa <- intersect(as.character(taxa), rownames(mtx))
      
    } else {
      stop ("Invalid argument for `taxa`: ", capture.output(str(taxa)))
    }
    
    
    if (length(taxa) == 0)
      stop("No taxa match the criteria: ", capture.output(str(taxa)))
    
    
    if (isFALSE(other) || is.null(other)) {
      
      mtx <- mtx[taxa,,drop=FALSE]
      
    } else {
      
      if (!is_scalar_character(other) || is_na(other))
        other <- "Other"
      
      mtx <- mtx[setdiff(rownames(mtx), taxa),,drop=FALSE] %>% 
        slam::col_sums() %>%
        matrix(., nrow = 1, dimnames = list(other, names(.))) %>%
        slam::as.simple_triplet_matrix() %>%
        rbind(mtx[taxa,,drop=FALSE], .)
      
      taxa <- c(taxa, other)
    }
    
  }
  
  
  
  #________________________________________________________
  # Convert from slam to base matrix if sparse=FALSE
  #________________________________________________________
  if (isFALSE(sparse))
    mtx <- as.matrix(mtx)
  
  
  
  set_cache_value(cache_file, mtx)
  
  return (mtx)
}
