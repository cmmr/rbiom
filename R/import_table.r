
#' Read a data.frame from a file/URL/string. Auto-detects field separator.
#' 
#' @noRd
#' @keywords internal
#' 
#' @param file  A file name, URL, or literal data string.
#'        File types supported: csv, tsv, json, xls, and xlsx.
#' 
#' @param ...  Not used.
#' 
#' @return A tibble data.frame.
#' 

import_table <- function (file, ...) {
  
  dots <- list(...)
  
  stopifnot(is_scalar_character(file))
  stopifnot(!is.na(file))
  stopifnot(nzchar(trimws(file)))
  
  
  #________________________________________________________
  # Import first sheet from an Excel file.
  #________________________________________________________
  if (grepl("\\.(xls|xlsx)$", tolower(file))) {
    df <- readxl::read_excel(file) # from tidyverse pkg
  }
  
  
  #________________________________________________________
  # Import from a text file, URL, or literal data string.
  #________________________________________________________
  else {
    
    # Heuristics for flagging as a json file or string.
    is_json <- endsWith(tolower(file), '.json')
    is_json <- is_json || startsWith(file, '[')
    if (!is_json && file.exists(file))
      is_json <- is_json || identical(readChar(file, 1L), '[')
    
    
    if (is_json) {
      df <- as_tibble(jsonlite::fromJSON(file))
      
    } else {
      
      # read_delim will auto-detect the field separator
      df <- readr::read_delim(
        file           = file,  
        trim_ws        = TRUE, 
        name_repair    = trimws,
        show_col_types = FALSE )
    }
  }
  
  
  #________________________________________________________
  # Enforce unique column names.
  #________________________________________________________
  
  if (any(x <- duplicated(colnames(df)))) {
    x <- unique(colnames(df)[x])
    cli_abort("Duplicated column names in {.file {file}}: {.val {x}}")
  }
  
  
  
  return (df)
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
  if (is_scalar_character(value)) value <- import_table(file = value)
  
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
    to_string <- function (i, x) jsonlite::toJSON(as.vector(x[i,,drop=FALSE]), auto_unbox = TRUE)
    
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
  if (is_scalar_character(value)) value <- import_table(file = value)
  
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
  # Drop unused factor levels.
  #________________________________________________________
  for (i in seq_len(ncol(value))[-1])
    value[[i]] %<>% { factor(., levels = intersect(levels(.), unique(.))) }
  
  
  return (value)
}

