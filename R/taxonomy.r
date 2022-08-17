#' Get the taxonomy table.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read_biom}.
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
#' @family accessors
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
    stop('In taxonomy(), biom must be a BIOM-class object.')
  
  if (length(missing <- setdiff(NULL, taxa_ranks(biom))) > 0)
    stop("Invalid ranks(s): ", paste(collapse = ", ", missing))
  
  
  x <- biom[['taxonomy']]
  
  
  
  
  # Fix all taxa names that are ambiguous keywords or not an 
  # upper-case letter followed by one or more lower case letters.
  #----------------------------------------------------------------
  if (isTRUE(fix.names)) {
    
    canonical <- "^[A-Z][a-z\ ]+$"
    ambiguous <- "unknown|uncultured|unclassified|unidentified|group|subsection|family|lineage|candidate"
    
    x <- sub("^.__", "", x) # Remove leading p__ c__ etc
    x <- sub("^_+", "", x)  # Remove leading underscores
    x <- sub(";$", "", x)   # Remove trailing semicolons
    
    isCanonical <- grepl(canonical, x, ignore.case = FALSE)
    isAmbiguous <- grepl(ambiguous, x, ignore.case = TRUE)
    isValid     <- matrix(nrow = nrow(x), isCanonical & !isAmbiguous)
    
    for (row in seq_len(nrow(x))) {
      
      lastValid     <- "Root"
      lastValidRank <- "Root"
      lastValidTaxa <- "Root"
      parenText     <- NULL
      
      for (col in seq_len(ncol(x))) {
        
        if (isValid[row, col]) {
          lastValidRank <- colnames(x)[col]
          lastValidTaxa <- x[row, col]
          
          lastValid <- paste(lastValidRank, lastValidTaxa)
          parenText <- NULL
          
        } else {
          
          parenText <- paste(c(parenText, x[row, col]), collapse = " ")
          
          # Prevent duplicate taxa names in final string
          parenText <- sub(
            pattern     = sprintf("^%s($|\\ |\\_)", lastValidTaxa), 
            replacement = "", 
            x           = parenText, 
            ignore.case = TRUE )
          
          # Remove '_group' suffix from names; remove incertae_sedis
          parenText <- sub("(_group|Incertae_Sedis)$", "", parenText)
          
          # Replace "g", etc
          parenText <- sub("(\\ [a-z])+$", "", parenText)
          if (isTRUE(nchar(parenText) <= 1))
            parenText <- "uncultured"
          
          x[row, col] <- sprintf("%s (%s)", lastValid, parenText)
        }
      }
    }
  }
  
  x <- cbind(x, 'OTU' = taxa_names(biom))
  if (!is.null(rank)) x <- x[,rank]
  
  return (x)
}




#' Set the taxonomy table.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read_biom}.
#' 
#' @param value  A character matrix with rownames \code{taxa_names(x)}. If
#'        there are more rownames than taxa names, the matrix will be subset.
#'         
#' @family setters
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

`taxonomy<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.matrix(value))
  stopifnot(typeof(value) == "character")
  stopifnot(all(taxa_names(x) %in% rownames(value)))
  
  x[['taxonomy']] <- value[taxa_names(x),]
  
  return (x) 
}


