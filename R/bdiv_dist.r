#' Make a distance matrix of samples vs samples.
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
#' @param stat.by  Specify a categorical metadata column name to compute adonis 
#'        statistics and return results as attributes named 'stats_raw' and 
#'        'stats_tbl'. (Default: \code{NULL})
#'        
#' @param perms  Number of random permutations to use in adonis calcuation.
#'        (Default: \code{999})
#'        
#' @param seed  Random seed for adonis permutations. (Default: \code{0})
#'        
#' @return A distance matrix.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- select(hmp50, 1:10)
#'     dm <- bdiv_distmat(biom, 'unifrac')
#'     as.matrix(dm)[1:4,1:4]
#'     plot(hclust(dm))
#'

bdiv_distmat <- function (biom, method="bray-curtis", weighted=TRUE, tree=NULL, stat.by=NULL, seed=0, perms=999) {
  
  #--------------------------------------------------------------
  # Enable abbreviations of metric names.
  #--------------------------------------------------------------
  method <- validate_metrics(biom, method, mode="bdiv", tree=tree)
  
  
  #--------------------------------------------------------------
  # Sanity Checks
  #--------------------------------------------------------------
  
  if (!is.logical(weighted)) stop(simpleError("Weighted must be TRUE/FALSE."))
  if (length(method) != 1)   stop(simpleError("Invalid method for bdiv_distmat()"))
  if (is.na(method))         stop(simpleError("Invalid method for bdiv_distmat()"))
  
  
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
  # Sanity check the matrix
  #--------------------------------------------------------------
  if (!is.numeric(counts[['v']]))     stop("The abundance matrix must be numeric.")
  if (!all(is.finite(counts[['v']]))) stop("Non-finite values in abundance matrix.")
  
  
  #--------------------------------------------------------------
  # Find the UniFrac tree
  #--------------------------------------------------------------
  
  if (method == "UniFrac") {
    
    # Find a tree for the UniFrac algorithm
    if (!is(tree, "phylo")) {
      if (is(biom, "BIOM")) {
        if (is(biom$phylogeny, "phylo")) {
          tree <- biom$phylogeny
        }
      }
      if (is(tree, "character")) {
        if (file.exists(tree)) {
          tree <- rbiom::read_tree(tree)
        }
      }
      if (!is(tree, "phylo")) {
        stop(simpleError("No tree provided to bdiv_distmat()."))
      }
    }
    
    
    # Make sure the matrix has rownames set
    if (is.null(rownames(counts)))
      stop("The abundance matrix does not have rownames set to Taxa IDs.")
    
    
    # Abundance matrix's Taxa IDs aren't all in the tree
    if (length(missing <- setdiff(rownames(counts), tree$tip.label)) > 0) {
      
      # Try swapping spaces/underscores in tip labels
      if (any(grepl("_", missing, fixed = TRUE))) {
        tree$tip.label %<>% gsub(" ", "_", ., fixed = TRUE)
      } else if (any(grepl(" ", missing, fixed = TRUE))) {
        tree$tip.label %<>% gsub("_", " ", ., fixed = TRUE)
      }
      
      if (!all(rownames(counts) %in% tree$tip.label))missing %>%
        glue::glue_collapse(sep = ", ", width = 30, last = " and ") %>%
        paste("OTUs missing from reference tree:", .) %>%
        stop()
    }
    
    # Prune the tree down to only what we need
    if (length(setdiff(tree$tip.label, rownames(counts))) > 0)
      tree <- rbiom::subtree(tree, rownames(counts))
    
    counts <- counts[as.character(tree$tip.label),]
  }
  
  
  #--------------------------------------------------------------
  # Order the sparse matrix's values by sample, then by taxa
  #--------------------------------------------------------------
  
  ord      <- order(counts$j, counts$i)
  counts$i <- counts$i[ord]
  counts$j <- counts$j[ord]
  counts$v <- counts$v[ord]
  
  
  #--------------------------------------------------------------
  # Run C++ implemented dissimilarity algorithms multithreaded
  #--------------------------------------------------------------
  
  if (method == "UniFrac") {
    
    dm <- par_unifrac(counts, tree, ifelse(weighted, 1L, 0L))
    
    
  } else {
    
    counts <- t(as.matrix(counts))
    dm <- par_beta_div(counts, tolower(method), ifelse(weighted, 1L, 0L))
    dm <- as.dist(dm)
    attr(dm, 'Labels') <- rownames(counts)
    
  }
  
  
  
  #-----------------------------------------------
  # Adonis / PERMANOVA
  #-----------------------------------------------
  if (!is.null(g <- stat.by))
    attr(dm, 'stats_raw') <- try(silent=TRUE, local({
      if (length(g) == 1)     g <- unname(metadata(biom, g))
      if (!is.null(names(g))) g <- unname(g[rownames(counts)])
      set.seed(seed)
      return(vegan::adonis2(dm ~ g, permutations=perms))
    }))
  attr(dm, 'stats_tbl') <- adonis_table(attr(dm, 'stats_raw', exact = TRUE))
  
  
  
  #--------------------------------------------------------------
  # Return the dist object
  #--------------------------------------------------------------
  return (dm)
}


adonis_table <- function (x) {
  
  if (!is(x, 'anova'))
    return (NULL)
  
  try(silent=TRUE, local({
    with(x, data.frame(
      check.names = FALSE,
      '.test'     = "Adonis",
      '.p.val'    = `Pr(>F)`[[1]],
      '.r.sqr'    = `R2`[[1]],
      '.f.stat'   = `F`[[1]]
    ))
  }))
}



