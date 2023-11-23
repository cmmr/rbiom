#' Make a data.frame of distances between samples.
#' 
#' @inherit documentation_cmp
#' @inherit documentation_default
#' 
#' @family beta_diversity
#'     
#' @param md  Include metadata in the output data frame? Options are: 
#'        \itemize{
#'          \item{\code{NULL} - }{ Don't include metadata. }
#'          \item{\code{TRUE} - }{ Include all metadata. }
#'          \item{\emph{character vector} - }{
#'            Include only the specified metadata columns. }
#'        }
#'        Default: \code{NULL}
#'        
#' @return A data.frame with columns named `.sample1`, `.sample2`, `.weighted`, 
#'         `.bdiv`, `.distance`, and any fields requested by `md`.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     # Return in long format with metadata
#'     biom <- sample_select(hmp50, 18:21)
#'     bdiv_table(biom, 'unifrac', md = c("Body Site", "Sex"))
#'     
#'     # Only look at distances among the stool samples
#'     bdiv_table(biom, 'unifrac', md = c("==Body Site", "Sex"))
#'     
#'     # Or between males and females
#'     bdiv_table(biom, 'unifrac', md = c("Body Site", "!=Sex"))
#'

bdiv_table <- function (
    biom, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, 
    md = NULL, within = NULL, between = NULL ) {
  
  validate_biom(clone = FALSE)
  validate_tree(null_ok = TRUE)
  
  params  <- eval_envir(environment())
  history <- append_history('tbl ', params)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  validate_bdiv(max = Inf)
  validate_bool('weighted', max = Inf)
  validate_meta('md', null_ok = TRUE, max = Inf, cmp = TRUE)
  validate_meta_cmp('md') # Validates and appends to `within` and `between`.
  
  
  
  #________________________________________________________
  # Multiple combinations of bdiv/weighted into one table.
  #________________________________________________________
  tbl <- NULL
  
  for (w in weighted)
    for (b in bdiv)
      tbl %<>% dplyr::bind_rows(local({
        
        #________________________________________________________
        # Compute the distance matrix
        #________________________________________________________
        dm <- bdiv_distmat(biom = biom, bdiv = b, weighted = w, tree = tree)
        
        
        #________________________________________________________
        # Convert to long form
        #________________________________________________________
        mtx <- as.matrix(dm)
        tibble(
          '.sample1'  = rownames(mtx)[row(mtx)],
          '.sample2'  = colnames(mtx)[col(mtx)],
          '.weighted' = w,
          '.bdiv'     = b,
          '.distance' = as.numeric(mtx) ) %>%
          dplyr::filter(.sample1 < .sample2)
      }))
  
  tbl[['.bdiv']] %<>% factor(levels = bdiv)
  
  
  
  
  #________________________________________________________
  # Limit to only within or between comparisons.
  #________________________________________________________
  for (col in c(within, between)) {
    map <- sample_metadata(biom, col)
    v1  <- map[tbl[['.sample1']]] %>% as.character()
    v2  <- map[tbl[['.sample2']]] %>% as.character()
    if (col %in% within)  tbl %<>% slice(which(v1 == v2))
    if (col %in% between) tbl %<>% slice(which(v1 != v2))
  }
  
  
  #________________________________________________________
  # Add metadata columns
  #________________________________________________________
  for (col in unique(c(md, within, between))) {
    
    map <- sample_metadata(biom, col)
    v1  <- map[tbl[['.sample1']]] %>% as.character()
    v2  <- map[tbl[['.sample2']]] %>% as.character()
    
    
    # Change "Male vs Male" --> "Male", and sort by factor
    # level, e.g. "Male vs Female" --> "Female vs Male".
    #________________________________________________________
    tbl[[col]] <- ifelse(
      test = (v1 == v2), 
      yes  = v1, 
      no   = ifelse(
        test = (v1 < v2), 
        yes  = paste(v1, "vs", v2), 
        no   = paste(v2, "vs", v1) ))
    
    
    # Keep factors as factors.
    #________________________________________________________
    if (is.factor(map)) {
      
      lvls <- levels(map)
      
      if (length(lvls) > 1)
        lvls %<>% {c(., apply(combn(., 2), 2L, paste, collapse=" vs "))}
      
      tbl[[col]] %<>% {factor(., levels = intersect(lvls, .))}
    }
    
  }
  
  
  lvls <- levels(biom[['metadata']][['.sample']])
  tbl[['.sample1']] %<>% {factor(., levels = intersect(lvls, .))}
  tbl[['.sample2']] %<>% {factor(., levels = intersect(lvls, .))}
  
  
  attr(tbl, 'response') <- ".distance"
  
  
  attr(tbl, 'history') <- history
  set_cache_value(cache_file, tbl)
  
  return (tbl)
}
