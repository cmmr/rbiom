#' Make a distance matrix of samples vs samples.
#' 
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' 
#' @param weighted  Take relative abundances into account. When 
#'        \code{weighted=FALSE}, only presence/absence is considered.
#'        Default: \code{TRUE}
#' 
#' @param bdiv  Beta diversity distance algorithm(s) to use. Options are:
#'        \code{"Bray-Curtis"}, \code{"Manhattan"}, \code{"Euclidean"}, 
#'        \code{"Jaccard"}, and \code{"UniFrac"}. For \code{"UniFrac"}, a 
#'        phylogenetic tree must be present in \code{biom} or explicitly 
#'        provided via \code{tree=}. Default: \code{"Bray-Curtis"}
#' 
#' @return A \code{dist}-class distance matrix.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- sample_select(hmp50, 1:10)
#'     dm   <- bdiv_distmat(biom, 'unifrac')
#'     as.matrix(dm)[1:4,1:4]
#'     plot(hclust(dm))
#'

bdiv_distmat <- function (biom, bdiv="Bray-Curtis", weighted=TRUE, tree=NULL) {
  
  #________________________________________________________
  # Take care not to cache filepath to tree.
  #________________________________________________________
  validate_biom(clone = FALSE)
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
  
  
  counts <- otu_matrix(biom, sparse = TRUE)
  
  
  #________________________________________________________
  # Find the UniFrac tree
  #________________________________________________________
  if (bdiv == "UniFrac") {
    
    # Find a tree for the UniFrac algorithm
    if (is.null(tree)) {
      if (!has_tree(biom))
        stop ("No tree provided to bdiv_distmat().")
      tree <- otu_tree(biom)
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



