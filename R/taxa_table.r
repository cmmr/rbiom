

#' Taxa abundances per sample.
#' 
#' \itemize{
#'   \item{`taxa_matrix()` - }{ Accepts a single `rank` and returns a matrix. }
#'   \item{`taxa_table()` - }{ Can accept more than one `rank` and returns a tibble data.frame.  }
#' }
#' 
#' 
#' @inherit documentation_default
#' 
#' @family taxa_abundance
#' 
#' @return 
#' \itemize{
#'   \item{`taxa_matrix()` - }{
#'     A numeric matrix with taxa as rows, and samples as columns. }
#'   \item{`taxa_table()` - }{
#'     A tibble data frame with column names .sample, .taxa, .abundance, and any requested by `md`. }
#' }
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom$ranks
#'     
#'     taxa_matrix(hmp50, 'Phylum')[1:4,1:6]
#'     
#'     taxa_table(hmp50, 'Phylum')

taxa_matrix <- function (
    biom, rank = -1, taxa = NULL, lineage = FALSE, 
    sparse = FALSE, unc = "singly", other = FALSE, trans = "none" ) {
  
  biom <- as_rbiom(biom)
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
  validate_var_choices('trans', c("none", "rank", "log", "log1p", "sqrt", "percent"))
  
  
  #________________________________________________________
  # Fetch taxonomy map with ambiguous names corrected.
  #________________________________________________________
  map <- taxa_map(
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
      
      counts <- biom$counts
      if (eq(trans, 'percent')) counts <- rescale_cols(counts)
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
  # Optionally transform the computed abundance values.
  #________________________________________________________
  if (trans %in% c("rank", "log", "log1p", "sqrt"))
    mtx$v <- do.call(`::`, list('base', trans))(mtx$v)
  
  
  #________________________________________________________
  # Convert from slam to base matrix if sparse=FALSE
  #________________________________________________________
  if (isFALSE(sparse))
    mtx <- as.matrix(mtx)
  
  
  
  set_cache_value(cache_file, mtx)
  
  return (mtx)
}





#' @rdname taxa_matrix
#' @export
taxa_table <- function (
    biom, rank = -1, taxa = NULL, lineage = FALSE, 
    md = ".all", unc = "singly", other = FALSE, trans = "none" ) {
  
  biom <- as_rbiom(biom)
  
  
  
  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  validate_rank(max = Inf)
  validate_meta('md', max = Inf, null_ok = TRUE)
  
  
  
  #________________________________________________________
  # Return multiple ranks in a single table.
  #________________________________________________________
  tbl <- NULL
  
  for (r in rank)
    tbl %<>% dplyr::bind_rows(local({
      
      mtx <- taxa_matrix(
        biom    = biom, 
        rank    = r, 
        taxa    = taxa, 
        lineage = lineage, 
        sparse  = FALSE, 
        unc     = unc, 
        other   = other,
        trans   = trans )
      
      
      #________________________________________________________
      # Pivot Longer
      #________________________________________________________
      tibble(
        '.rank'      = r,
        '.sample'    = colnames(mtx)[col(mtx)],
        '.taxa'      = rownames(mtx)[row(mtx)],
        '.abundance' = as.numeric(mtx) )
      
    }))
  
  tbl[['.rank']]   %<>%  factor(., levels = rank)
  tbl[['.sample']] %<>% {factor(., levels = intersect(biom$samples, .))}
  tbl[['.taxa']]   %<>% {factor(., levels = unique(.))}
  
  
  
  #________________________________________________________
  # Add Metadata
  #________________________________________________________
  if (length(md) > 0)
    tbl %<>% left_join( 
      by = '.sample',
      y  = biom$metadata[,unique(c('.sample', md))] )
  
  
  
  attr(tbl, 'response') <- ".abundance"
  
  
  return (tbl)
}




#' Get summary taxa abundances.
#' 
#' @inherit documentation_default
#' 
#' @family taxa_abundance
#' 
#' @param rank  The taxonomic rank to return sums or means for. The default, 
#'        \code{0}, returns per-OTU summaries.
#'        
#' @return A named, sorted numeric vector.
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     taxa_sums(hmp50) %>% head(4)
#'     
#'     taxa_means(hmp50, 'Family') %>% head(5)
#'

taxa_sums <- function (biom, rank = 0) {
  
  taxa_matrix(biom, rank, sparse = TRUE) %>%
    slam::row_sums() %>%
    sort(decreasing = TRUE)
}



#' @rdname taxa_sums
#' @export

taxa_means <- function (biom, rank = 0) {
  
  taxa_matrix(biom, rank, sparse = TRUE) %>%
    slam::row_means() %>%
    sort(decreasing = TRUE)
}
