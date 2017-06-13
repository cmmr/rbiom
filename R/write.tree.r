#' Write a newick formatted phylogenetic tree.
#' 
#' @param tree  A phylo object, as returned from \link{read.tree}.
#' @param file  Filename or connection to write the newick file to (optional).
#' @return If file is NULL, the newick string as a character vector. Otherwise,
#'         the return value from \link[base]{writeChar}, typically invsible(NULL).
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree   <- read.tree(infile)
#'     newick <- write.tree(tree)
#'
write.tree <- function (tree=NULL, file=NULL) {
  
  if (is.null(tree))
    stop(simpleError("Please provide a value for tree to write.tree()"))
  
  if (!is(tree, "phylo"))
    stop(simpleError("Provided tree is not a phylo class object."))
  
  
  rootNode <- setdiff(tree$edge[,1], tree$edge[,2])
  parentAt <- aggregate(1:nrow(tree$edge), by=list(tree$edge[,1]), c, simplify=FALSE)
  parentAt <- setNames(lapply(parentAt[,2], unlist), parentAt[,1])
  
  fx <- function (root=NULL) {
    
    nodes <- parentAt[[as.character(root)]]
    
    if (length(nodes) == 0)
      return (tree$tip.label[root])
    
    children <- tree$edge[nodes, 2]
    children <- sapply(children, fx)
    
    if (!is.null(tree$edge.length))
      children <- paste(sep=":", children, tree$edge.length[nodes])
    
    sprintf("(%s)", paste(collapse=",", children))
  }
  
  newick <- paste0(fx(rootNode), ";")
  
  
  if (!is.null(file))
    return (writeChar(object=newick, con=file))
  
  return (newick)
}
