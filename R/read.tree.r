#' Read a newick formatted phylogenetic tree.
#' 
#' A phylogenetic tree is required for computing UniFrac distance matrices.
#' You can load a tree either from a file or by providing the tree string
#' directly. This tree must be in Newick format, also known as parenthetic
#' format and New Hampshire format.
#'
#' @param src Input data as either a file path, URL, or Newick string. URLs 
#'     must begin with \kbd{http://}, \kbd{https://}, \kbd{ftp://}, or 
#'     \kbd{ftps://}. Newick strings must have \code{(} as their first 
#'     non-whitespace character. Compressed (gzip or bzip2) Newick files 
#'     are also supported.
#' @return A \code{phylo} class object representing the tree.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree <- read.tree(infile)
#'     
#'     tree <- read.tree("
#'         (t9:0.99,((t5:0.87,t2:0.89):0.51,(((t10:0.16,(t7:0.83,t4:0.96)
#'         :0.94):0.69,(t6:0.92,(t3:0.62,t1:0.85):0.54):0.23):0.74,t8:0.1
#'         2):0.43):0.67);")
#'
read.tree <- function (src) {
  
  #--------------------------------------------------------------
  # Sanity check src value
  #--------------------------------------------------------------
  errmsg <- "In read.tree(), 'src' must be a character vector of length 1, not"
  
  if (length(src) != 1)      stop(simpleError(paste(errmsg, "length", length(src))))
  if (is.null(src))          stop(simpleError(paste(errmsg, "value 'NULL'.")))
  if (is.na(src))            stop(simpleError(paste(errmsg, "value 'NA'.")))
  if (!is(src, "character")) stop(simpleError(paste(errmsg, "class", class(src))))
  if (!nchar(src))           stop(simpleError(paste(errmsg, "empty string ''.")))
  
  
  #--------------------------------------------------------------
  # Shortened src string for use in error messages
  #--------------------------------------------------------------
  if (nchar(src) <= 20) {
    shortsrc <- src
  } else {
    shortsrc <- paste0(substr(src, 1, 10), "...", substr(src, nchar(src)-10, nchar(src)))
  }
  
  
  #--------------------------------------------------------------
  # Get the Newick data into a string
  #--------------------------------------------------------------
  
  if (length(grep("^[ \t\n]*\\(", src)) == 1) {
    
    text <- src
    
  } else {
    
    if (length(grep("^(ht|f)tps{0,1}://.+", src)) == 1) {
      
      fp <- tempfile(fileext=basename(src))
      on.exit(unlink(fp), add=TRUE)
      if (!identical(0L, try(download.file(src, fp, quiet=TRUE), silent=TRUE)))
        stop(simpleError(sprintf("Cannot retrieve URL %s", src)))
      
    } else {
      
      fp <- try(normalizePath(src, mustWork = TRUE), silent = TRUE)
      
      if (is(fp, "try-error"))
        stop(simpleError(sprintf("Cannot locate file '%s': %s", shortsrc, fp)))
      
    }
    
    
    #--------------------------------------------------------------
    # Uncompress files that are in gzip or bzip2 format
    #--------------------------------------------------------------
    
    file_con   <- file(fp)
    file_class <- summary(file_con)$class
    close.connection(file_con)
    
    if (file_class %in% c("gzfile", "bzfile")) {
      
      if (identical(file_class, "gzfile"))
        fp <- R.utils::gunzip(fp, destname=tempfile(), remove=FALSE)
      
      if (identical(file_class, "bzfile"))
        fp <- R.utils::bunzip2(fp, destname=tempfile(), remove=FALSE)
      
      on.exit(unlink(fp), add=TRUE)
    }
    
    remove("file_con", "file_class")
    
    
    #--------------------------------------------------------------
    # Read file contents into 'text' variable
    #--------------------------------------------------------------
    text <- try(readLines(con = fp, warn = FALSE), silent = TRUE)
    
  }
  
  
  #--------------------------------------------------------------
  # Sanity check the data string
  #--------------------------------------------------------------
  errmsg <- "Error in read.tree(): Unable to load data from '%s'."
  errmsg <- sprintf(errmsg, shortsrc)
  
  if (length(text) != 1)                    stop(simpleError(errmsg))
  if (is.null(text))                        stop(simpleError(errmsg))
  if (is.na(text))                          stop(simpleError(errmsg))
  if (is(text, "error"))                    stop(simpleError(paste(errmsg, text)))
  if (is(text, "try-error"))                stop(simpleError(paste(errmsg, text)))
  if (!identical(class(text), "character")) stop(simpleError(errmsg))
  
  # Remove newlines, comments, and leading whitespace
  text <- gsub("[\ \t]*[\r\n]+[\ \t]*", "", text)
  text <- gsub("\\[.*?\\]", "",             text, perl=TRUE)
  text <-  sub("^[\ \t]+",  "",             text)
  
  if (nchar(text) < 2)                      stop(simpleError(errmsg))
  if (!identical(substr(text, 1, 1), "("))  stop(simpleError(errmsg))
  
  
  #--------------------------------------------------------------
  # Parse the Newick string into a phylo object
  #--------------------------------------------------------------
  
  tree <- rcpp_read_tree(text)
  
  if (all(nchar(tree$node.label) == 0))
    tree$node.label <- NULL
  
  if (all(tree$edge.length == 0))
    tree$edge.length <- NULL
  
  return (tree)
}
