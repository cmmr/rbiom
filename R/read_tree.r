#' Read a newick formatted phylogenetic tree.
#' 
#' A phylogenetic tree is required for computing UniFrac distance matrices.
#' You can load a tree from a file or by providing the tree string directly. 
#' This tree must be in Newick format, also known as parenthetic format and
#' New Hampshire format.
#' 
#' @family phylogeny
#'
#' @param src   Input data as either a file path, URL, or Newick string. 
#'        Compressed (gzip or bzip2) files are also supported.
#'        
#' @return A \code{phylo} class object representing the tree.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree <- read_tree(infile)
#'     print(tree)
#'     
#'     tree <- read_tree("
#'         (A:0.99,((B:0.87,C:0.89):0.51,(((D:0.16,(E:0.83,F:0.96)
#'         :0.94):0.69,(G:0.92,(H:0.62,I:0.85):0.54):0.23):0.74,J:0.1
#'         2):0.43):0.67);")
#'     plot(tree)
#'
read_tree <- function (src) {
  
  stopifnot(is_scalar_character(src) && !is_na(src))
  
  src <- trimws(src)
  stopifnot(isTRUE(nchar(src) > 0))
  
  
  #________________________________________________________
  # Shortened src string for use in error messages.
  #________________________________________________________
  text <- src
    
  if (nchar(src) > 20)
    src <- paste0(
      substr(src, 1, 10), 
      "...", 
      substr(src, nchar(src)-10, nchar(src)) )
  
  
  #________________________________________________________
  # Get the Newick data into a string
  #________________________________________________________
  if (!startsWith(text, '(') && length(text) == 1)
    if (file.exists(text) || grepl('^(http|ftp)s{0,1}://', tolower(text)))
      text <- tryCatch(
        error = function (e) cli_abort("Can't read file {.file {src}}: {e}"), 
        expr  = local({
          con <- gzfile(text)
          on.exit(close(con))
          paste0(collapse = "", readLines(con, warn = FALSE)) }))
  
  stopifnot(is_scalar_character(text) && !is_na(text))
  
  
  #________________________________________________________
  # Remove newlines, comments, and leading whitespace
  #________________________________________________________
  text <- gsub("[\ \t]*[\r\n]+[\ \t]*", "", text)
  text <- gsub("\\[.*?\\]", "",             text, perl=TRUE)
  text <- sub("^[\ \t]+",  "",              text)
  
  stopifnot(isTRUE(nchar(text) >= 2))
  stopifnot(isTRUE(startsWith(text, '(')))
  
  
  #________________________________________________________
  # Parse the Newick string into a phylo object
  #________________________________________________________
  
  tree        <- .Call(C_read_tree, text)
  names(tree) <- c('edge', 'Nnode', 'tip.label', 'edge.length', 'node.label')
  attr(tree, 'class') <- 'phylo'
  attr(tree, 'order') <- 'cladewise'
  
  if (all(nchar(tree$node.label) == 0))
    tree$node.label <- NULL
  
  if (all(tree$edge.length == 0))
    tree$edge.length <- NULL
  
  return (tree)
}



#' Create a subtree by specifying tips to keep.
#' 
#' @name tree_subset
#' 
#' @family phylogeny
#' 
#' @param tree  A phylo object, as returned from [read_tree()].
#' @param tips  A character, numeric, or logical vector of tips to keep.
#' @return A \code{phylo} object for the subtree.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree <- read_tree(infile)
#'     tree
#'     
#'     subtree <- tree_subset(tree, tips = head(tree$tip.label))
#'     subtree
#'
tree_subset <- function (tree, tips) {
  
  validate_tree()
  params <- eval_envir(environment())
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('tree_subset', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Did we get a tree?
  #________________________________________________________
  
  if (!inherits(tree, "phylo"))
    cli_abort("Can't convert {.type {tree}} to a `phylo` object.")
  
  nTips <- length(setdiff(tree$edge[,2], tree$edge[,1]))
  
  
  #________________________________________________________
  # Tips specified as strings c("SteSpe52", "HrbSpe65", ...)
  #________________________________________________________
  
  if (is.character(tips)) {
    
    if (length(x <- setdiff(tips, tree$tip.label)))
      cli_abort("Tips missing from tree: {x}.")
    
    tips <- which(tree$tip.label %in% tips)
  }
  
  
  #________________________________________________________
  # Tips specified as logicals c(TRUE, TRUE, FALSE, TRUE, ...)
  #________________________________________________________
  
  if (is.logical(tips)) {
    
    if (length(tips) != nTips)
      cli_abort("Logical vector `tips` must be same length as `tree$tip.label`.")
    
    tips <- which(tips)
  }
  
  
  #________________________________________________________
  # Tips specified as numeric c(5, 1, 2, 31, 3, 12, ...)
  #________________________________________________________
  
  if (is.numeric(tips)) {
    if (length(tips) > 0) {
      if (max(tips) > nTips)  cli_abort("`tips` index {.val {max(tips)}} doesn't exist.")
      if (min(tips) < 1)      cli_abort("`tips` index {.val {min(tips)}} doesn't exist.")
      if (any(tips %% 1 > 0)) cli_abort("`tips` indices must be whole numbers.")
    }
    
  } else {
    cli_abort("`tips` must character, logical, or numeric, not {.type {tips}}.")
  }
  
  
  #________________________________________________________
  # Keeping zero or all tips?
  #________________________________________________________
  
  tips <- unique(tips)
  
  if (length(tips) < 2)
    cli_abort("At least two tips must be specified.")
  
  if (length(tips) == nTips) {
    set_cache_value(cache_file, tree)
    return (tree)
  }
  
  
  #________________________________________________________
  # Prune the tree
  #________________________________________________________
  
  root      <- setdiff(tree$edge[,1], tree$edge[,2])
  childAt   <- order(c(tree$edge[,2], root))
  nChildren <- c(rep(0, nTips), unname(table(tree$edge[,1])))
  
  eLength   <- tree$edge.length
  if (is_null(eLength))
    eLength <- rep(0, nrow(tree$edge))
  
  for (tip in setdiff(1:nTips, tips)) repeat {
    
    tipIdx <- childAt[tip]
    parent <- tree$edge[tipIdx,1]
    
    nChildren[parent] <- nChildren[parent] - 1
    eLength[tipIdx] <- NA
    
    if (nChildren[parent] > 0) break
    tip <- parent
  }
  
  
  #________________________________________________________
  # Merge branches
  #________________________________________________________
  
  for (tip in tips) repeat {
    
    if (tip == root) break
    
    tipIdx <- childAt[tip]
    parent <- tree$edge[tipIdx,1]
    
    if (nChildren[parent] == 1 && parent != root) {
      
      parentIdx                   <- childAt[parent]
      tree$edge[parentIdx, 2]     <- tree$edge[tipIdx, 2]
      eLength[parentIdx] <- eLength[parentIdx] + eLength[tipIdx]
      eLength[tipIdx]    <- NA
      
      nChildren[parent] <- 0
    }
    
    tip <- parent
  }
  
  
  #________________________________________________________
  # Move the root
  #________________________________________________________
  
  if (nChildren[root] == 1)
    eLength[which(tree$edge[,1] == root)] <- NA
  
  
  #________________________________________________________
  # Discard nodes flagged with NA; make the new subtree
  #________________________________________________________
  
  tree$edge  <- tree$edge[!is.na(eLength),]
  tree$Nnode <- length(unique(tree$edge[,1]))
  
  if (!is_null(tree$edge.length)) tree$edge.length <- eLength[!is.na(eLength)]
  if (!is_null(tree$tip.label))   tree$tip.label   <- tree$tip.label[tips]
  if (!is_null(tree$node.label))  tree$node.label  <- tree$tip.label[unique(sort(tree$edge[,1]))]
  
  tree$edge <- matrix(as.numeric(factor(as.vector(tree$edge))), ncol=2)
  
  
  set_cache_value(cache_file, tree)
  return (tree)
}
