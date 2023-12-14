#' Distance / dissimilarity between samples.
#' 
#' @inherit documentation_cmp
#' @inherit documentation_default
#' 
#' @family beta_diversity
#'     
#' @param md  Include metadata in the output data frame? Options are: 
#'        \itemize{
#'          \item{`NULL` - }{ Don't include metadata. }
#'          \item{`TRUE` - }{ Include all metadata. }
#'          \item{\emph{character vector} - }{
#'            Include only the specified metadata columns.
#'            Prefix the column name(s) with `==` or `!=` to limit comparisons 
#'            to within or between groups, respectively. }
#'        }
#'        Default: `NULL`
#' 
#' @param weighted  Take relative abundances into account. When 
#'        `weighted=FALSE`, only presence/absence is considered.
#'        `bdiv_table()` can accept multiple values.
#'        Default: `TRUE`
#' 
#' @param bdiv  Beta diversity distance algorithm to use. Options are:
#'        `"Bray-Curtis"`, `"Manhattan"`, `"Euclidean"`, 
#'        `"Jaccard"`, and `"UniFrac"`. `bdiv_table()` can accept multiple 
#'        values. Default: `"Bray-Curtis"`
#' 
#' @return
#' \itemize{
#'   \item{`bdiv_distmat()` - }{ A dist-class distance matrix. }
#'   \item{`bdiv_table()` - }{
#'     A tibble data.frame with columns names .sample1, .sample2, .weighted, 
#'     .bdiv, .distance, and any fields requested by `md`. }
#' }
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     # Subset to four samples
#'     biom <- hmp50$clone()
#'     biom$counts <- biom$counts[,c("HMP18", "HMP19", "HMP20", "HMP21")]
#'     
#'     # Return in long format with metadata
#'     bdiv_table(biom, 'unifrac', md = c("Body Site", "Sex"))
#'     
#'     # Only look at distances among the stool samples
#'     bdiv_table(biom, 'unifrac', md = c("==Body Site", "Sex"))
#'     
#'     # Or between males and females
#'     bdiv_table(biom, 'unifrac', md = c("Body Site", "!=Sex"))
#'     
#'     # Distance matrix for all samples
#'     dm <- bdiv_distmat(biom, 'unifrac')
#'     as.matrix(dm)
#'     plot(hclust(dm))
#'

bdiv_distmat <- function (biom, bdiv="Bray-Curtis", weighted=TRUE, tree=NULL) {
  
  #________________________________________________________
  # Take care not to cache filepath to tree.
  #________________________________________________________
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- eval_envir(environment())
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  remove("params")
  
  
  #________________________________________________________
  # Check for valid arguments.
  #________________________________________________________
  validate_bool("weighted")
  validate_bdiv()
  
  
  counts <- biom$counts
  
  
  #________________________________________________________
  # Find the UniFrac tree
  #________________________________________________________
  if (bdiv == "UniFrac") {
    
    # Find a tree for the UniFrac algorithm
    if (is.null(tree)) {
      if (is.null(biom$tree))
        stop ("No tree provided to bdiv_distmat().")
      tree <- biom$tree
    }
    
    
    # Make sure the matrix has rownames set
    if (is_null(rownames(counts)))
      stop("The abundance matrix must have OTU names as row names.")
    
    
    # Abundance matrix's Taxa IDs aren't all in the tree
    if (length(missing <- setdiff(rownames(counts), tree$tip.label)) > 0) {
      
      # Try swapping spaces/underscores in tip labels
      if (any(grepl("_", missing, fixed = TRUE))) {
        tree$tip.label %<>% gsub(" ", "_", ., fixed = TRUE)
      } else if (any(grepl(" ", missing, fixed = TRUE))) {
        tree$tip.label %<>% gsub("_", " ", ., fixed = TRUE)
      }
      
      if (!all(rownames(counts) %in% tree$tip.label))
        missing %>%
          glue_collapse(sep = ", ", width = 30, last = " and ") %>%
          paste("OTUs missing from reference tree:", .) %>%
          stop()
    }
    
    # Prune the tree down to only what we need
    if (length(setdiff(tree$tip.label, rownames(counts))) > 0)
      tree <- tree_subset(tree, rownames(counts))
    
    counts <- counts[as.character(tree$tip.label),]
  }
  
  
  
  #________________________________________________________
  # Order the sparse matrix's values by sample, then by taxa
  #________________________________________________________
  
  ord      <- order(counts$j, counts$i)
  counts$i <- counts$i[ord]
  counts$j <- counts$j[ord]
  counts$v <- counts$v[ord]
  
  
  
  #________________________________________________________
  # Run C++ implemented dissimilarity algorithms multithreaded
  #________________________________________________________
  
  if (bdiv == "UniFrac") {
    
    dm <- par_unifrac(counts, tree, ifelse(weighted, 1L, 0L))
    
    
  } else {
    
    counts <- t(as.matrix(counts))
    dm <- par_beta_div(counts, tolower(bdiv), ifelse(weighted, 1L, 0L))
    dm <- as.dist(dm)
    attr(dm, 'Labels') <- rownames(counts)
    
  }
  
  
  
  #________________________________________________________
  # Return the dist object
  #________________________________________________________
  set_cache_value(cache_file, dm)
  return (dm)
}






#' @rdname bdiv_distmat
#' @export
bdiv_table <- function (
    biom, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, 
    md = NULL, within = NULL, between = NULL ) {
  
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)
  
  params <- eval_envir(environment())
  
  
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
    map <- pull(biom, col)
    v1  <- map[tbl[['.sample1']]] %>% as.character()
    v2  <- map[tbl[['.sample2']]] %>% as.character()
    if (col %in% within)  tbl %<>% slice(which(v1 == v2))
    if (col %in% between) tbl %<>% slice(which(v1 != v2))
  }
  
  
  #________________________________________________________
  # Add metadata columns
  #________________________________________________________
  for (col in unique(c(md, within, between))) {
    
    map <- pull(biom, col)
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
  
  
  lvls <- biom$metadata[['.sample']]
  tbl[['.sample1']] %<>% {factor(., levels = intersect(lvls, .))}
  tbl[['.sample2']] %<>% {factor(., levels = intersect(lvls, .))}
  
  
  attr(tbl, 'response') <- ".distance"
  
  
  set_cache_value(cache_file, tbl)
  
  return (tbl)
}

