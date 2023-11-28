
#' Abundance table of taxa per sample at the specified rank.
#' 
#' @inherit documentation_default
#' 
#' @family taxa_abundance
#'        
#' @return A tibble data frame with column names \code{.sample}, \code{.taxa}, 
#'         \code{.abundance}, and any requested by \code{md}.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     taxa_ranks(hmp50)
#'     
#'     head(taxa_table(hmp50, 'Phylum'))
#'

taxa_table <- function (
    biom, rank = -1, taxa = NULL, lineage = FALSE, 
    md = ".all", unc = "singly", other = FALSE ) {
  
  validate_biom(clone = FALSE)
  
  params  <- eval_envir(environment())
  history <- append_history('tbl ', params)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  validate_rank(max = Inf)
  validate_meta('md', max = Inf)
  
  
  
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
        other   = other )
      
      
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
  tbl[['.sample']] %<>% {factor(., levels = intersect(sample_levels(biom), .))}
  tbl[['.taxa']]   %<>% {factor(., levels = unique(.))}
  
  
  
  #________________________________________________________
  # Add Metadata
  #________________________________________________________
  if (length(md) > 0)
    tbl %<>% left_join( 
      by = '.sample',
      y  = sample_metadata(biom)[,unique(c('.sample', md))] )
  
  
  
  attr(tbl, 'response') <- ".abundance"
  
  
  
  attr(tbl, 'history') <- history
  set_cache_value(cache_file, tbl)
  
  return (tbl)
}

