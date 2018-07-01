#' Names of a phylogenetic tree's tips/leafs.
#' 
#' @param x  A phylo object, as returned from \link{read.tree}..
#' @return A character vector with the leaf names.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree <- read.tree(infile)
#'     
#'     leafs   <- tips(tree)
#'     subtree <- subtree(tree, head(leafs))
#'
tips <- function (x) {
  x$tip.label
}
