


#' Get the sample metadata.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read_biom}.
#' 
#' @param field  The name of a single metadata column to retrieve. If provided,
#'        a named vector will be returned instead of a data.frame.
#'        (Default: NULL)
#'        
#' @param id  Copy the sample names (rownames) into this column. 
#'        (Default: NULL)
#'        
#' @param cleanup  Change character columns to factor or dates. Date formats 
#'        recognized are: "%Y-%m-%d" or "%m/%d/%Y" (four-digit year). Also
#'        '.' with '_' in column names when '.' is the first character in order
#'        to avoid conflicts with auto-generated columns. (Default: FALSE)
#'        
#' @return A data frame of the metadata in \code{biom}. If \code{field} is
#'    given, will return a named vector of the field's values.
#'    
#' @family accessors
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     metadata(hmp50)[1:4,1:3]
#'     
#'     head(metadata(hmp50, "Body Site"))
#'     
#'     head(metadata(hmp50, "Body Site", cleanup=TRUE))
#'

metadata <- function (biom, field=NULL, id=NULL, cleanup=FALSE) {
  
  if (!is(biom, 'BIOM'))
    stop (simpleError('In metadata(), biom must be a BIOM-class object.'))
  
  md <- biom[['metadata']]
  
  
  #________________________________________________________
  # Add 'SampleID', '.id', etc as the first column
  #________________________________________________________
  if (!is_null(id)) {
    md[[id]] <- rownames(md)
    md       <- md[,unique(c(id, names(md))),drop=FALSE]
  }
  
  if (isTRUE(cleanup)){
    
    #________________________________________________________
    # Rename columns starting with a '.'
    #________________________________________________________
    names(md) <- sub("^\\.", "_", names(md))
  
    #________________________________________________________
    # Set columns to factor/numeric/date
    #________________________________________________________
    for (i in seq_len(ncol(md))) {
      
      vals <- as.vector(na.omit(md[[i]]))
      
      if (is.numeric(vals)) {
        md[[i]] <- as.numeric(md[[i]])
        
      } else if (all(sapply(vals, grepl, pattern="^\\d{4}\\-\\d{1,2}\\-\\d{1,2}$"))) {
        md[[i]] <- as.Date(md[[i]], "%Y-%m-%d")
        
      } else if (all(sapply(vals, grepl, pattern="^\\d{1,2}\\/\\d{1,2}\\/\\d{4}$"))) {
        md[[i]] <- as.Date(md[[i]], "%m/%d/%Y")
        
      } else {
        md[[i]] <- as.factor(md[[i]])
        
      }
    }
  }
  
  if (is_null(field))
    return (md)
  
  if (!field %in% names(md))
    stop(paste0("Field '", field, "' is not present in the metadata."))
  
  return (setNames(md[[field]], rownames(md)))
}





#' Set the \code{BIOM} object's metadata.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read_biom}.
#' 
#' @param value  A data.frame with the metadata. All the \code{rownames()} must 
#'        be in \code{sample_names(biom)}. If there are fewer rows of 
#'        \code{data} than samples, then the \code{biom} object will be subset.
#'    
#' @family setters
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     md <- metadata(hmp50)
#'     md <- md[,c('Sex', 'Body Site')]
#'     metadata(hmp50) <- md
#'     head(metadata(hmp50))
#'
`metadata<-` <- function(x, value) {
  
  stopifnot(is(x, 'BIOM'))
  
  
  #________________________________________________________
  # Parse a file/URL.
  #________________________________________________________
  if (is_scalar_character(value))
    value <- import_table(value, row.names = 1, header = TRUE)
  
  stopifnot(is.data.frame(value))
  
  
  #________________________________________________________
  # All current sample names must be in the new metadata.
  #________________________________________________________
  missing <- setdiff(sample_names(x), rownames(value))
  n       <- length(missing)
  if (n > 0) {
    if (n > 4) missing <- c(head(missing, 4), "...")
    missing <- paste(collapse = ", ", missing)
    msg <- "%i Sample ID%s missing from the new metadata: %s"
    stop(sprintf(msg, n, ifelse(n == 1, " is", "s are"), missing))
  }
  
  
  #________________________________________________________
  # Ignore IDs that aren't currently in the BIOM object.
  #________________________________________________________
  ignored <- setdiff(rownames(value), sample_names(x))
  n       <- length(ignored)
  if (n > 0) {
    if (n > 4) ignored <- c(head(ignored, 4), "...")
    ignored <- paste(collapse = ", ", ignored)
    msg <- "%i Sample ID%s have no matching sample in the BIOM: %s"
    warning(sprintf(msg, n, ifelse(n == 1, " is", "s are"), ignored))
  }
  
  
  
  
  x[['metadata']] <- value[sample_names(x),,drop=FALSE]
  
  return (x)
}
