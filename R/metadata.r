


#' Get the sample metadata.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
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
#' @family accessor functions
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
  
  
  #-----------------------------------------------
  # Add 'SampleID', '.id', etc as the first column
  #-----------------------------------------------
  if (!is.null(id)) {
    md[[id]] <- rownames(md)
    md       <- md[,unique(c(id, names(md))),drop=FALSE]
  }
  
  if (isTRUE(cleanup)){
    
    #-----------------------------------------------
    # Rename columns starting with a '.'
    #-----------------------------------------------
    names(md) <- sub("^\\.", "_", names(md))
  
    #-----------------------------------------------
    # Set columns to factor/numeric/date
    #-----------------------------------------------
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
  
  if (is.null(field))
    return (md)
  
  if (!field %in% names(md))
    stop(paste0("Field '", field, "' is not present in the metadata."))
  
  return (setNames(md[[field]], rownames(md)))
}