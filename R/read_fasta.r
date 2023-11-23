
#' Parse a fasta file into a named character vector.
#'
#' @param file   A file/URL with fasta-formatted sequences. Can optionally be 
#'        compressed with gzip, bzip2, xz, or lzma.
#' 
#' @param ids   Character vector of IDs to retrieve. The default, \code{NULL}, 
#'        will retrieve everything.
#' 
#' @return A named character vector in which names are the fasta headers and 
#'         values are the sequences.
#' 
#' @export
#'

read_fasta <- function (file, ids = NULL) {
  
  if (!file.exists(file))
    return(simpleError(paste0("File '", file, "' does not exist.")))
  
  
  #________________________________________________________
  # Open uncompressed files as well as gzip/bzip2/xz/lzma.
  #________________________________________________________
  con <- gzfile(file, "r")
  
  
  #________________________________________________________
  # Parse headers and potentially multi-line sequences.
  #________________________________________________________
  sid <- NULL
  res <- c()
  
  while (length(line <- readLines(con, 1)) > 0) {
    
    if (eq(substr(line, 1, 1), ">")) {
      sid <- substr(line, 2, nchar(line))
      
      if (sid %in% names(res)) {
        close(con)
        return(simpleError(paste0("Sequence ID '", sid, "' is used more than once in '", file, "'.")))
      }
      
      if (!is_null(ids))
        if (!sid %in% ids)
          sid <- NULL
      
    } else if (!is_null(sid)) {
      
      if (sid %in% names(res)) {
        res[[sid]] <- paste0(res[[sid]], line)
      } else {
        res[[sid]] <- line
      }
    }
    
  }
  
  close(con)
  
  
  #________________________________________________________
  # Ensure we return as many IDs as requested
  #________________________________________________________
  if (!is_null(ids)) {
    res <- res[ids]
    names(res) <- ids
  }
  
  
  return (res)
}
