#' Compute Weighted and Unweighted UniFrac distance matrices.
#' 
#' This is the function called internally by \link{beta.div}, but is made
#' visible here so you can use it with matrices and trees without having to
#' first convert them to \code{BIOM} objects.
#'
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'     object, as returned from \link{read.biom}. For matrices, the rows and 
#'     columns are assumed to be the taxa and samples, respectively.
#' @param weighted  Use weighted UniFrac, which takes abundance into account
#'     rather than simply presence/absence.
#' @param tree  A \code{phylo} object providing a phylogenetic tree for the
#'     taxa names in \code{biom}. If \code{tree=NULL}, then the tree will be
#'     loaded from \code{biom}, if encoded there.
#' @return A distance matrix of class \code{dist}.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     biom <- select(biom, 1:10)
#'     
#'     dm <- unifrac(biom)
#'     plot(hclust(dm), cex=.8)
#'     as.matrix(dm)[1:4,1:4]
#'     
#'     # Using a custom matrix and tree
#'     mtx <- matrix(sample.int(12*20), ncol=20)
#'     dimnames(mtx) <- list(LETTERS[1:12], letters[1:20])
#'     tree <- ape::as.phylo(hclust(dist(mtx)))
#'     
#'     dm <- unifrac(mtx, tree=tree)
#'     as.matrix(dm)[1:4,1:4]
#'     

unifrac <- function (biom, weighted=TRUE, tree=NULL) {
  
  rbiom::beta.div(biom, 'unifrac', weighted, tree)
  
}

