
#' Parse a fasta file into a named character vector.
#'
#' @param file   A file/URL with fasta-formatted sequences. Can optionally be 
#'        compressed with gzip, bzip2, xz, or lzma.
#' 
#' @param ids   Character vector of IDs to retrieve. The default, `NULL`, 
#'        will retrieve everything.
#' 
#' @return A named character vector in which names are the fasta headers and 
#'         values are the sequences.
#' 
#' @export
#'

read_fasta <- function (file, ids = NULL) {
  
  if (!file.exists(file))
    cli_abort("Fasta file {.file {file}} does not exist.")
  
  
  #________________________________________________________
  # Open uncompressed files as well as gzip/bzip2/xz/lzma.
  #________________________________________________________
  con <- gzfile(file, "r")
  on.exit(close(con), add = TRUE)
  
  
  #________________________________________________________
  # Parse headers and potentially multi-line sequences.
  #________________________________________________________
  sid <- NULL
  res <- c()
  
  while (length(line <- readLines(con, 1)) > 0) {
    
    if (eq(substr(line, 1, 1), ">")) {
      sid <- substr(line, 2, nchar(line))
      
      if (sid %in% names(res))
        cli_abort(c(
          '!' = "Fasta file {.file {file}}", 
          'x' = 'Has multiple entries for {.val {sid}}.' ))
      
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
  
  
  #________________________________________________________
  # Ensure we return as many IDs as requested
  #________________________________________________________
  if (!is_null(ids)) {
    
    if (length(x <- setdiff(ids, names(res))))
      cli_abort(c(
        '!' = "In fasta file {.file {file}}", 
        'x' = 'Expected IDs are missing: {.val {x}}.' ))
    
    res <- res[ids]
  }
  
  
  return (res)
}
