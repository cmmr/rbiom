#' Make a distance matrix of samples vs samples.
#' 
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'     object, as returned from \link{read.biom}. For matrices, the rows and 
#'     columns are assumed to be the taxa and samples, respectively.
#' @param method  The distance algorithm to use. Options are:
#'     \bold{\dQuote{manhattan}}, \bold{\dQuote{euclidean}}, 
#'     \bold{\dQuote{bray-curtis}}, \bold{\dQuote{jaccard}}, and
#'     \bold{\dQuote{unifrac}}. Non-ambiguous abbrevations of the method 
#'     names are also accepted. A phylogentic tree must be present in 
#'     \code{biom} or explicitly provided via \code{tree=} to use the UniFrac methods.
#' @param weighted  Take relative abundances into account. When 
#'     \code{weighted=FALSE}, only presence/absence is considered.
#' @param tree  A \code{phylo} object representing the phylogenetic
#'     relationships of the taxa in \code{biom}. Will be taken from the tree
#'     embedded in the \code{biom} object if not explicitly specified. Only
#'     required for computing UniFrac distance matrices.
#' @return A distance matrix.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     biom <- select(biom, 1:10)
#'     
#'     dm <- beta.div(biom, 'unifrac')
#'     
#'     as.matrix(dm)[1:4,1:4]
#'     plot(hclust(dm))
#'

beta.div <- function (biom, method, weighted=TRUE, tree=NULL) {
  
  #--------------------------------------------------------------
  # Enable abbreviations of metric names.
  #--------------------------------------------------------------
  
  methodList <- c("manhattan", "euclidean", "bray-curtis", "jaccard", "unifrac")
  method <- methodList[pmatch(tolower(method), methodList)]
  
  
  #--------------------------------------------------------------
  # Sanity Checks
  #--------------------------------------------------------------
  
  if (!is.logical(weighted)) stop(simpleError("Weighted must be TRUE/FALSE."))
  if (length(method) != 1)   stop(simpleError("Invalid method for beta.div()"))
  if (is.na(method))         stop(simpleError("Invalid method for beta.div()"))
  
  
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
  # Find the UniFrac tree
  #--------------------------------------------------------------
  
  if (identical(method, "unifrac")) {
    if (!is(tree, "phylo")) {
      if (is(biom, "BIOM")) {
        if (is(biom$phylogeny, "phylo")) {
          tree <- biom$phylogeny
        }
      }
      if (is(tree, "character")) {
        if (file.exists(tree)) {
          tree <- rbiom::read.tree(tree)
        }
      }
      if (!is(tree, "phylo")) {
        stop(simpleError("No tree provided to beta.div()."))
      }
    }
    
    if (length(setdiff(rownames(counts), tree$tip.label)) > 0)
      stop(simpleError("OTUs missing from reference tree."))
    
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
  
  if (identical(method, "unifrac")) {
    
    par_unifrac(counts, tree, ifelse(weighted, 1L, 0L))
    
    
  } else {
    
    counts <- t(as.matrix(counts))
    dm <- par_beta_div(counts, method, ifelse(weighted, 1L, 0L))
    dm <- as.dist(dm)
    attr(dm, 'Labels') <- rownames(counts)
    
    return (dm)
    
    # rcpp_parallel_js_distance(as.matrix(counts))
    
  }
  
  
}





