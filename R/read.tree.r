#' Read a newick formatted phylogenetic tree.
#' 
#' A phylogenetic tree is required for computing UniFrac distance matrices.
#' You can load a tree either from a file or by providing the tree string
#' directly.
#'
#' @param file  Path to the newick file.
#' @param text  The newick tree as a string.
#' @return A \code{phylo} class object.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree <- read.tree(infile)
#'
read.tree <- function (file=NULL, text=NULL) {
  
  if (is.null(text)) {
    
    if (is.null(file))
      stop(simpleError("Please provide a value to read.tree()"))
    
    if (!file.exists(file))
      stop(simpleError(sprintf("Can't find file '%s'", file)))
    
    text <- readChar(file, file.size(file))
  }
  
  # Remove newlines and comments
  text <- gsub("\n",        "", text, fixed=TRUE)
  text <- gsub("\r",        "", text, fixed=TRUE)
  text <- gsub("\\[.*?\\]", "", text, perl=TRUE)
  
  tree <- rcpp_read_tree(text)
  
  if (all(nchar(tree$node.label) == 0))
    tree$node.label <- NULL
  
  if (all(tree$edge.length == 0))
    tree$edge.length <- NULL
  
  return (tree)
}
