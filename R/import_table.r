
#' Read a data.frame from a filename/URL. Auto-detects separator.
#' 
#' @noRd
#' 
#' @param filename  A file name or URL. Must be tab- or comma-separated 
#'        content. Either all or no fields should be quoted.
#' 
#' @param matrix  Set to \code{TRUE} to convert data.frame to a matrix. Set to
#'        \code{"character"}, \code{"integer"}, etc to coerce to a matrix of
#'        that type. The default, \code{FALSE}, does no conversion.
#' 
#' @param ...  Additional arguments for \code{read.table()}.
#' 
#' @return A data.frame or matrix.
#' 

import_table <- function (filename, matrix = FALSE, ...) {
  
  read_table_args <- list(...)
  stopifnot(is_scalar_character(filename) && !is.na(filename))
    
    
  #________________________________________________________
  # Download from a URL.
  #________________________________________________________
  if (grepl("^.{3,5}://", filename)) {
    
    tmp <- tempfile()
    on.exit(unlink(tmp), add=TRUE)
    
    res <- tryCatch(
      expr  = download.file(filename, tmp, quiet=TRUE),
      error = function (e) stop("Can't download ", filename,"\n", e) )
    if (!identical(res, 0L) || !file.exists(tmp))
      stop("Download failed for ", filename)
    
    filename <- tmp
  }
  
  
  
  #________________________________________________________
  # Import from a file.
  #________________________________________________________
  if (!file.exists(filename))
    stop("File not found: ", filename)
  
  
  
  #________________________________________________________
  # See if first line has more tabs or commas.
  #________________________________________________________
  line  <- strsplit(readLines(filename, 1L), '')[[1]]
  sep   <- ifelse(sum(line == "\t") > sum(line == ","), "\t", ',')
  
  read_table_args[['file']]         <-  filename
  read_table_args[['sep']]         %<>% if.null(sep)
  read_table_args[['check.names']] %<>% if.null(FALSE)
  
  df <- tryCatch(
    expr  = do.call(read.table, read_table_args),
    error = function (e) stop("Can't parse file ", filename, "\n", e) )
  
  
  
  #________________________________________________________
  # Enforce unique row/column names.
  #________________________________________________________
  
  if (any(x <- duplicated(rownames(df)))) {
    x <- unique(rownames(df)[x])
    if (length(x) > 4) x <- c(head(x, 4), "...")
    msg <- "Duplicated row names in %s: %s"
    stop(sprintf(msg, filename, paste(collapse = ", ", x)))
  }
  
  if (any(x <- duplicated(colnames(df)))) {
    x <- unique(colnames(df)[x])
    if (length(x) > 4) x <- c(head(x, 4), "...")
    msg <- "Duplicated column names in %s: %s"
    stop(sprintf(msg, filename, paste(collapse = ", ", x)))
  }
  
  
  
  
  #________________________________________________________
  # Coerce to matrix.
  #________________________________________________________
  
  if (isFALSE(matrix)) return (df)
  
  mtx <- tryCatch(
      expr  = as(df, 'matrix'),
      error = function (e)
        stop("Can't convert ", filename, " to matrix.\n", e) )
  
  if (is_scalar_character(matrix) && !typeof(mtx) == matrix)
    mtx[] <- tryCatch(
      expr  = as(mtx, matrix),
      error = function (e)
        stop("Can't coerce ", filename, " to ", matrix, ".\n", e) )
  
  if (!isTRUE(read_table_args[['header']]))    colnames(mtx) <- NULL
  if (is.null(read_table_args[['row.names']])) rownames(mtx) <- NULL
  
  
  
  return (mtx)
}




