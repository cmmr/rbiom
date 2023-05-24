#' Read a newick formatted phylogenetic tree.
#' 
#' A phylogenetic tree is required for computing UniFrac distance matrices.
#' You can load a tree from a file or by providing the tree string directly. 
#' This tree must be in Newick format, also known as parenthetic format and
#' New Hampshire format.
#'
#' @param src Input data as either a file path, URL, or Newick string. 
#'        Compressed (gzip or bzip2) files are also supported.
#'        
#' @return A \code{phylo} class object representing the tree.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree <- read_tree(infile)
#'     
#'     tree <- read_tree("
#'         (t9:0.99,((t5:0.87,t2:0.89):0.51,(((t10:0.16,(t7:0.83,t4:0.96)
#'         :0.94):0.69,(t6:0.92,(t3:0.62,t1:0.85):0.54):0.23):0.74,t8:0.1
#'         2):0.43):0.67);")
#'
read_tree <- function (src) {
  
  stopifnot(is_scalar_character(src))
  
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
  if (!startsWith(text, '('))
    text <- tryCatch(
      error = function (e) stop("Can't read file ", src, ".\n", e), 
      expr  = local({
        con <- gzfile(text)
        on.exit(close(con))
        paste0(collapse = "", readLines(con, warn = FALSE)) }))
  
  stopifnot(is_scalar_character(text))
  
  
  #________________________________________________________
  # Remove newlines, comments, and leading whitespace
  #________________________________________________________
  text <- gsub("[\ \t]*[\r\n]+[\ \t]*", "", text)
  text <- gsub("\\[.*?\\]", "",             text, perl=TRUE)
  text <- sub("^[\ \t]+",  "",              text)
  
  stopifnot(nchar(text) >= 2)
  stopifnot(startsWith(text, '('))
  
  
  #________________________________________________________
  # Parse the Newick string into a phylo object
  #________________________________________________________
  
  tree <- rcpp_read_tree(text)
  
  if (all(nchar(tree$node.label) == 0))
    tree$node.label <- NULL
  
  if (all(tree$edge.length == 0))
    tree$edge.length <- NULL
  
  return (tree)
}
