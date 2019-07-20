#' Parse a fasta file into a named character vector.
#'
#' @param file  A file with fasta-formatted sequences. Can optionally be compressed with gzip, bzip2, xz, or lzma.
#' @param ids  Character vector of IDs to retrieve. The default, \code{NULL}, will retrieve everything.
#' @return A named character vector in which names are the fasta headers and values are the sequences.
#' @export
#'


read.fasta <- function (file, ids=NULL) {
  
  if (!file.exists(file))
    return(simpleError(paste0("File '", file, "' does not exist.")))
  
  
  #--------------------------------------------------------------
  # Open uncompressed files as well as gzip, bzip2, xz, or lzma.
  #--------------------------------------------------------------
  con <- gzfile(file, "r")
  
  
  #--------------------------------------------------------------
  # Parse headers and potentially multiline sequences
  #--------------------------------------------------------------
  sid <- NULL
  res <- c()
  
  while (length(line <- readLines(con, 1)) > 0) {
    
    if (identical(substr(line, 1, 1), ">")) {
      sid <- substr(line, 2, nchar(line))
      
      if (sid %in% names(res)) {
        close(con)
        return(simpleError(paste0("Sequence ID '", sid, "' is used more than once in '", file, "'.")))
      }
      
      if (!is.null(ids))
        if (!sid %in% ids)
          sid <- NULL
      
    } else if (!is.null(sid)) {
      
      if (sid %in% names(res)) {
        res[[sid]] <- paste0(res[[sid]], line)
      } else {
        res[[sid]] <- line
      }
    }
    
  }
  
  close(con)
  
  
  #--------------------------------------------------------------
  # Ensure we return as many IDs as requested
  #--------------------------------------------------------------
  if (!is.null(ids)) {
    res <- res[ids]
    names(res) <- ids
  }
  
  
  return (res)
}
