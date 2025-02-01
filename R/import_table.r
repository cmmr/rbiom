
#' Read a data.frame from a file/URL/string. Auto-detects field separator.
#' 
#' @noRd
#' @keywords internal
#' 
#' @param src  A file name, URL, or literal data string.
#'        File types supported: csv, tsv, json, xls, and xlsx.
#' 
#' @return A tibble data.frame.
#' 

import_table <- function (src) {
  
  #________________________________________________________
  # Get the data into a file.
  #________________________________________________________
  fpath  <- as_filepath(src)
  format <- fpath$format
  path   <- fpath$path
  on.exit(fpath$cleanup(), add = TRUE)
  
  
  #________________________________________________________
  # Read json/text/excel data into a data frame.
  #________________________________________________________
  tbl <- if (format == 'excel') {
    readxl::read_excel(path = path)
  
  } else if (format == 'json') {
    as_tibble(jsonlite::read_json(path = path, simplifyVector = TRUE))
    
  } else if (format == 'text') {
    
    # read_delim auto-detects encoding and separator. 
    # Don't use read.table - it chokes on UTF-8 symbols.
    readr::read_delim(
      file           = path,  
      trim_ws        = TRUE, 
      name_repair    = trimws,
      comment        = '#q2:', # QIIME 2 numeric or categorical
      show_col_types = FALSE )
  
    
  } else {
    cli_abort("`src` data type not recognized: {src}")
  }
  
  
  #________________________________________________________
  # Enforce unique column names.
  #________________________________________________________
  
  if (any(x <- duplicated(colnames(tbl)))) {
    x <- unique(colnames(tbl)[x])
    cli_abort("Duplicated column names in {.file {src}}: {.val {x}}")
  }
  
  
  return (tbl)
}


#' Determine if a file path, URL, or literal data string.
#' is 'text', 'json', 'hdf5', or 'excel' format.
#' 
#' @noRd
#' @keywords internal
#' 
#' @param src  A file path, URL, or literal data string.
#' @param filename  The string to examine for file extensions.
#' 
#' @return 'text', 'json', 'hdf5', or 'excel'. 
#'         'AsIs' class when `src` is literal data rather than a file.
#' 

file_type <- function (src, filename = src) {
  
  if (!is_nz_string(src)) cli_abort('`src` must be a string, not {.type {src}}.')
  if (is_url(src))        cli_abort('`src` cannot be a URL: {.url {src}}.')
  
  if (startsWith(src, '{') || startsWith(src, '[')) return (I('json'))
  if (grepl('\n', src, fixed = TRUE))               return (I('text'))
  
  if (!file.exists(src)) cli_abort('File not found: {filename}.')
  
  x <- sub('\\.(gz|bz2)$', '', tolower(filename))
  if (grepl('\\.(xls|xlsx)$', x))        return ('excel')
  if (grepl('\\.(txt|tsv|csv|tab)$', x)) return ('text')
  if (grepl('\\.json$', x))              return ('json')
  if (grepl('\\.(hdf|hdf5)$', x))        return ('hdf5')
  
  con  <- base::file(src)
  cls  <- summary(con)$class
  peek <- suppressWarnings(readChar(con, 1L))
  close(con)
  
  if (cls %in% c('gzfile', 'bzfile'))
    return (ifelse(peek %in% c('{', '['), 'json', 'text'))
  
  con  <- base::file(src, 'rb', raw = TRUE)
  peek <- readBin(con, 'raw', 8L)
  close(con)
  
  if (peek[[1]] %in% charToRaw('{['))                                      return ('json')
  if (identical(peek[1:4], charToRaw('\x89HDF')))                          return ('hdf5')
  if (identical(peek[1:4], charToRaw('\x50\x4B\x03\x04')))                 return ('excel')
  if (identical(peek[1:8], charToRaw('\xD0\xCF\x11\xE0\xA1\xB1\x1A\xE1'))) return ('excel')
  
  return ('text')
}



#' Pull URLs to a local file. Save literal data to a file.
#' `$rm()` will remove any files created by this function.
#' 
#' @noRd
#' @keywords internal
#' 
#' @param src  A file path, URL, or literal data string.
#' 
#' @return list(id, path, format, cleanup)
#' 

as_filepath <- function (src) {
    
  if (is_url(src)) {
    
    id      <- if (nchar(src) > 100) basename(src) else src
    path    <- tempfile()
    cleanup <- function () unlink(path)
    if (!eq(0L, x <- try(download.file(url = src, path, quiet=TRUE), silent=TRUE)))
      cli_abort("Cannot retrieve URL {.url {src}}: {x}")
    
    format <- file_type(src = path, filename = src)
    
    
  } else {
    
    format <- file_type(src = src)
    
    if (inherits(format, 'AsIs')) {
      
      # Save literal data string to a temp file.
      id      <- paste('Parsed from', format, 'string')
      path    <- tempfile()
      cleanup <- function () unlink(path)
      writeChar(object = src, con = path, eos = NULL)
      class(format) <- NULL
      
    } else {
      path    <- normalizePath(src, winslash = '/')
      id      <- if (nchar(path) > 100) basename(path) else path
      cleanup <- function () NULL
    }
  }
  
  result <- list(id = id, path = path, format = format, cleanup = cleanup)
  
  
  return (result)
}




#' @noRd
#' @keywords internal
import_metadata <- function (value, sids) {
  
  #________________________________________________________
  # Erasing all the metadata.
  #________________________________________________________
  if (is.null(value)) return (tibble(.sample = sids))
  
  #________________________________________________________
  # Load metadata from a file.
  #________________________________________________________
  if (is_scalar_character(value)) value <- import_table(src = value)
  
  #________________________________________________________
  # Convert data.frame/matrix to tibble.
  #________________________________________________________
  value <- tibble::as_tibble(value, rownames = NA, .name_repair = trimws)
  
  
  #________________________________________________________
  # Auto-detect the column with sample names.
  #________________________________________________________
  if (!tibble::has_rownames(value) && !hasName(value, '.sample'))
    for (i in seq_len(ncol(value)))
      if (!as.logical(anyDuplicated(vals <- trimws(value[[i]]))))
        if (any(vals %in% sids)) {
          colnames(value)[i] <- '.sample'
          break
        }
  
  
  #________________________________________________________
  # Make '.sample' the first column.
  #________________________________________________________
  if (!tibble::has_rownames(value) && !hasName(value, '.sample'))
    cli_abort(c(x = "Metadata table must have row names or a '.sample' column."))
  
  if (tibble::has_rownames(value) && hasName(value, '.sample')) {
    if (!all(rownames(value) == paste0('', seq_len(nrow(value)))))
      if (!all(rownames(value) == value$.sample))
        cli_abort(c(x = "Row names are not allowed when a '.sample' column is present."))
    value %<>% tibble::remove_rownames()
  }
  
  if (tibble::has_rownames(value)) { value %<>% tibble::rownames_to_column('.sample')
  } else                           { value %<>% relocate('.sample') }
  
  
  
  #________________________________________________________
  # Disallow column names starting with a "."
  #________________________________________________________
  if (length(i <- grep("^\\.", colnames(value)[-1], value = TRUE)) != 0)
    cli_abort(c(x = "Metadata field{?s} can't start with a '.': {.val {i}}."))
  
  
  #________________________________________________________
  # Convert all character cols (except .sample) to factor.
  #________________________________________________________
  value[['.sample']] %<>% as.character() %>% trimws()
  for (i in colnames(value)[-1])
    if (is.character(value[[i]])) value[[i]] %<>% trimws() %>% as.factor()
  
  
  #________________________________________________________
  # Intersect with sample names from $counts.
  #________________________________________________________
  expected <- sids
  provided <- unique(value[['.sample']])
  
  if (length(intersect(expected, provided)) == 0)
    cli_abort(c(
      'i' = "No matching sample names.", 
      '*' = "Expected: {.val {expected}}", 
      '*' = "Provided: {.val {provided}}", 
      'x' = "Can't subset to zero samples." ))
  
  if (length(i <- setdiff(expected, provided)) > 0)
    cli_warn(c('i' = paste(
      "Dropping {length(i)} sample{?s} from biom object",
      "since they are not in the new metadata: {.val {i}}." )))

  if (length(i <- setdiff(provided, expected)) > 0)
    cli_warn(c('i' = paste(
      "Ignoring metadata for {length(i)} sample{?s}", 
      "not currently in biom object: {.val {i}}." )))
  
  value <- value[value[['.sample']] %in% expected,,drop=FALSE]
  remove('expected', 'provided', 'i')
  
  
  #________________________________________________________
  # Drop duplicate sample mappings. Warn if non-identical.
  #________________________________________________________
  if (as.logical(anyDuplicated(value[['.sample']]))) {
    
    not_equal <- function (x) (length(unique(as.character(x))) > 1)
    to_string <- function (i, x) toJSON(as.vector(x[i,,drop=FALSE]), auto_unbox = TRUE)
    
    plyr::d_ply(value, '.sample', function (x) {
      
      if (nrow(x) == 1) return (NULL)
      sid <- x[1,1]
      x   <- x[,apply(x, 2L, not_equal),drop=FALSE]
      if (ncol(x) == 0) return (NULL)
      
      rlang::warn(
        message = c(i = paste0("Discarding subsequent metadata for sample ", sid)),
        body    = setNames(
          object = sapply(1:nrow(x), to_string, x), 
          nm     = c('v', rep('x', nrow(x) - 1)) ))
    })
   
    value <- value[!duplicated(value[['.sample']]),,drop=FALSE]
    remove("not_equal", "to_string")
  }
  
  
  return (value)
}


#' @noRd
#' @keywords internal
import_taxonomy <- function (value, otus) {
  
  #________________________________________________________
  # Erasing all the taxonomy.
  #________________________________________________________
  if (is.null(value)) return (tibble(.otu = otus))
  
  
  #________________________________________________________
  # Load taxonomy from a file.
  #________________________________________________________
  if (is_scalar_character(value)) value <- import_table(src = value)
  
  #________________________________________________________
  # Convert data.frame/matrix to tibble.
  #________________________________________________________
  value <- tibble::as_tibble(value, rownames = NA, .name_repair = trimws)
  
  
  #________________________________________________________
  # Auto-detect the column with OTU names.
  #________________________________________________________
  if (!tibble::has_rownames(value) && !hasName(value, '.otu'))
    for (i in seq_len(ncol(value)))
      if (!as.logical(anyDuplicated(vals <- trimws(value[[i]]))))
        if (any(vals %in% otus)) {
          colnames(value)[i] <- '.otu'
          break
        }
  
  
  #________________________________________________________
  # Make '.otu' the first column.
  #________________________________________________________
  if (!tibble::has_rownames(value) && !hasName(value, '.otu'))
    cli_abort(c(x = "Taxonomy table must have row names or an '.otu' column."))
  
  if (tibble::has_rownames(value) && hasName(value, '.otu')) {
    if (!all(rownames(value) == paste0('', seq_len(nrow(value)))))
      if (!all(rownames(value) == value$.otu))
        cli_abort(c(x = "Row names are not allowed when an '.otu' column is present."))
    value %<>% tibble::remove_rownames()
  }
  
  if (tibble::has_rownames(value)) { value %<>% tibble::rownames_to_column('.otu')
  } else                           { value %<>% relocate('.otu') }
  
  
  #________________________________________________________
  # Disallow column names starting with a "."
  #________________________________________________________
  if (length(i <- grep("^\\.", colnames(value)[-1], value = TRUE)) != 0)
    cli_abort(c(x = "Taxonomy rank{?s} can't start with a '.': {.val {i}}."))
  
  
  #________________________________________________________
  # Convert all cols (except .otu) to factor.
  #________________________________________________________
  value[['.otu']] %<>% as.character() %>% trimws()
  for (i in seq_len(ncol(value))[-1])
    if (!is.factor(value[[i]])) value[[i]] %<>% trimws() %>% as.factor()
  
  
  #________________________________________________________
  # Intersect with OTU names from $counts.
  #________________________________________________________
  expected <- otus
  provided <- unique(value[['.otu']])
  
  if (length(intersect(expected, provided)) == 0)
    cli_abort(c(
      'i' = "No matching OTU names between counts and taxonomy.", 
      '*' = "Expected: {.val {expected}}", 
      '*' = "Provided: {.val {provided}}", 
      'x' = "Can't subset to zero OTUs." ))
  
  if (length(i <- setdiff(expected, provided)) > 0)
    cli_warn(c('i' = paste(
      "Dropping {length(i)} OTU{?s} from biom object",
      "since they are not in the new taxonomy: {.val {i}}." )))
  
  if (length(i <- setdiff(provided, expected)) > 0)
    cli_warn(c('i' = paste(
      "Ignoring taxonomy for {length(i)} OTU{?s}", 
      "not currently in biom object: {.val {i}}." )))
  
  value <- value[value[['.otu']] %in% expected,,drop=FALSE]
  remove('expected', 'provided', 'i')
  
  
  #________________________________________________________
  # Drop duplicate OTU mappings. Warn if non-identical.
  #________________________________________________________
  if (as.logical(anyDuplicated(value[['.otu']]))) {
    
    if (ncol(value) == 1) {
      keep <- which(!duplicated(value[['.otu']]))
      
    } else {
      
      keep <- sapply(split(1:nrow(value), value[['.otu']]), function (i) {
        
        if (length(i) == 1) return (i)
        
        lineages <- apply(value[i,-1,drop=FALSE], 1L, paste, collapse = "; ")
        longest  <- which.max(nchar(lineages))
        
        if (length(unique(lineages)) > 1)
          rlang::warn(
            message = c(i = paste0("Discarding less verbose mapping for OTU ", value[i[[1]],1])),
            body    = setNames(
              object = lineages, 
              nm     = ifelse(seq_along(lineages) == longest, 'v', 'x') ))
        
        return (i[[longest]])
      })
    }
    
    value <- value[sort(unname(keep)),,drop=FALSE]
    remove("keep")
  }
  
  
  #________________________________________________________
  # Drop 'Confidence' column (from qiime2).
  #________________________________________________________
  if (identical(tail(colnames(value), 1), 'Confidence'))
    value <- value[,seq_len(ncol(value) - 1)]
  
  
  #________________________________________________________
  # Split semicolon-separated values into columns.
  #________________________________________________________
  if (ncol(value) == 2) {
    ssv <- paste0(collapse = '', as.character(value[[2]]), '\n')
    ssv <- readr::read_delim(file = ssv, delim = ';', quote = '', col_names = FALSE, col_types = 'f')
    colnames(ssv) <- default_taxa_ranks(ncol(ssv))
    value <- dplyr::bind_cols(value[,'.otu'], ssv)
    remove('ssv')
  }
  
  
  #________________________________________________________
  # Drop unused factor levels.
  #________________________________________________________
  for (i in seq_len(ncol(value))[-1])
    value[[i]] %<>% { factor(., levels = intersect(levels(.), unique(.))) }
  
  
  return (value)
}

