#' Create a subtree by specifying tips to keep.
#' 
#' @param tree  A phylo object, as returned from \link{read.tree}.
#' @param tips  A character, numeric, or logical vector of tips to keep.
#' @return A \code{phylo} object for the subtree.
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
subtree <- function (tree, tips) {
  
  #--------------------------------------------------------------
  # Did we get a tree?
  #--------------------------------------------------------------
  
  if (!is(tree, "phylo"))
    stop(simpleError("Must provide a phylo object as the tree."))
  
  nTips <- length(setdiff(tree$edge[,2], tree$edge[,1]))
  
  
  #--------------------------------------------------------------
  # Tips specified as strings c("SteSpe52", "HrbSpe65", ...)
  #--------------------------------------------------------------
  
  if (is.character(tips)) {
    
    if (!all(tips %in% tree$tip.label)) {
      missing <- tips[!tips %in% tree$tip.label]
      if (length(missing) > 5) missing <- c(missing[1:5], "...")
      missing <- paste(collapse=", ", missing)
      msg <- sprintf("Tips missing from tree: %s", missing)
      stop(simpleError(msg))
    }
    
    tips <- which(tree$tip.label %in% tips)
  }
  
  
  #--------------------------------------------------------------
  # Tips specified as logicals c(TRUE, TRUE, FALSE, TRUE, ...)
  #--------------------------------------------------------------
  
  if (is.logical(tips)) {
    
    if (length(tips) != nTips)
      stop(simpleError("Logical vector must be same length as tips."))
    
    tips <- which(tips)
  }
  
  
  #--------------------------------------------------------------
  # Tips specified as numeric c(5, 1, 2, 31, 3, 12, ...)
  #--------------------------------------------------------------
  
  if (is.numeric(tips)) {
    
    if (max(tips) > nTips)
      stop(simpleError("There aren't that many tips in the tree."))
    
    if (min(tips) < 1)
      stop(simpleError("Tips with index < 1 are not allowed."))
    
    if (any(tips %% 1 > 0))
      stop(simpleError("Tip indices must be whole numbers."))
    
  } else {
    stop(simpleError("Tips must be given as character, logical, or numeric."))
  }
  
  
  #--------------------------------------------------------------
  # Keeping zero or all tips?
  #--------------------------------------------------------------
  
  tips <- unique(tips)
  
  if (length(tips) < 2)
    stop(simpleError("At least two tips must be specified."))
  
  if (length(tips) == nTips)
    return (tree)
  
  
  #--------------------------------------------------------------
  # Prune the tree
  #--------------------------------------------------------------
  
  root      <- setdiff(tree$edge[,1], tree$edge[,2])
  childAt   <- order(c(tree$edge[,2], root))
  nChildren <- c(rep(0, nTips), unname(table(tree$edge[,1])))
  
  eLength   <- tree$edge.length
  if (is.null(eLength))
    eLength <- rep(0, nrow(tree$edge))
  
  for (tip in setdiff(1:nTips, tips)) repeat {
    
    tipIdx <- childAt[tip]
    parent <- tree$edge[tipIdx,1]
    
    nChildren[parent] <- nChildren[parent] - 1
    eLength[tipIdx] <- NA
    
    if (nChildren[parent] > 0) break
    tip <- parent
  }
  
  
  #--------------------------------------------------------------
  # Merge branches
  #--------------------------------------------------------------
  
  for (tip in tips) repeat {
    
    if (tip == root) break
    
    tipIdx <- childAt[tip]
    parent <- tree$edge[tipIdx,1]
    
    if (nChildren[parent] == 1 & parent != root) {
      
      parentIdx                   <- childAt[parent]
      tree$edge[parentIdx, 2]     <- tree$edge[tipIdx, 2]
      eLength[parentIdx] <- eLength[parentIdx] + eLength[tipIdx]
      eLength[tipIdx]    <- NA
      
      nChildren[parent] <- 0
    }
    
    tip <- parent
  }
  
  
  #--------------------------------------------------------------
  # Move the root
  #--------------------------------------------------------------
  
  if (nChildren[root] == 1)
    eLength[which(tree$edge[,1] == root)] <- NA
  
  
  #--------------------------------------------------------------
  # Discard nodes flagged with NA; make the new subtree
  #--------------------------------------------------------------
  
  tree$edge  <- tree$edge[!is.na(eLength),]
  tree$Nnode <- length(unique(tree$edge[,1]))
  
  if (!is.null(tree$edge.length)) tree$edge.length <- eLength[!is.na(eLength)]
  if (!is.null(tree$tip.label))   tree$tip.label   <- tree$tip.label[tips]
  if (!is.null(tree$node.label))  tree$node.label  <- tree$tip.label[unique(sort(tree$edge[,1]))]
  
  tree$edge <- matrix(as.numeric(factor(as.vector(tree$edge))), ncol=2)
  
  return (tree)
  
}

