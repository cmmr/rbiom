#' Get the taxonomy table.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' 
#' @param rank  Return just a named character vector of this rank. (Default: 
#'        \code{NULL})
#' 
#' @param fix.names  Clean up taxa names. See example below. (Default: 
#'        \code{FALSE})
#'                   
#' @return A character matrix of the named taxonomies in \code{biom}.
#'         When \code{fix.names = TRUE}, a second matrix is attached as an
#'         attribute called 'fixed' and notes which taxa names were
#'         rewritten.
#'         
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     taxonomy(hmp50)[1:4,]
#'     
#'     # Sometimes taxonomic names are incomplete
#'     taxonomy(hmp50)[c(107,170,194), 'Genus', drop=FALSE]
#'     
#'     # rbiom can insert more descriptive placeholders
#'     taxonomy(hmp50, fix.names = TRUE)[c(107,170,194), 'Genus', drop=FALSE]
#'

taxonomy <- function (biom, rank = NULL, fix.names = FALSE) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In taxonomy(), biom must be a BIOM-class object.'))
  
  if (!identical(fix.names, TRUE) || ncol(biom[['taxonomy']]) == 0)
    return (biom[['taxonomy']][,rank])
  
  
  # Fix all taxa names that are ambiguous keywords or not an 
  # upper-case letter followed by one or more lower case letters.
  #----------------------------------------------------------------
  canonical <- "^[A-Z][a-z\ ]+$"
  ambiguous <- "unknown|uncultured|unclassified|unidentified|group|subsection|family|lineage|candidate"
  
  x <- biom[['taxonomy']]
  x <- sub("^.__", "", x) # Remove leading p__ c__ etc
  x <- sub("^_+", "", x)  # Remove leading underscores
  x <- sub(";$", "", x)   # Remove trailing semicolons
  
  isCanonical <- grepl(canonical, x, ignore.case = FALSE)
  isAmbiguous <- grepl(ambiguous, x, ignore.case = TRUE)
  isValid     <- matrix(nrow = nrow(x), isCanonical & !isAmbiguous)
  
  for (row in seq_len(nrow(x))) {
    
    lastValid <- "Root"
    parenText <- NULL
    
    for (col in seq_len(ncol(x))) {
      
      if (isValid[row, col]) {
        lastValid <- paste(colnames(x)[col], x[row, col])
        parenText <- NULL
        
      } else {
        parenText   <- paste(c(parenText, x[row, col]), collapse = " ")
        x[row, col] <- sprintf("%s (%s)", lastValid, parenText)
      }
    }
  }
  
  x <- x[,rank]
  attr(x, 'fixed') <- isValid
  
  return (x)
}
