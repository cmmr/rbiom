#' Get or set the sample metadata.
#' 
#' @inherit documentation_default
#' 
#' @family setters
#' 
#' @param field   The name of a single metadata column to retrieve. When 
#'        provided, a named vector will be returned instead of a data.frame.
#'        Default: \code{NULL}
#' 
#' @param value  A data.frame with the metadata. Its row names or '.sample' 
#'        column must match the IDs in \code{sample_names(x)}. The returned
#'        rbiom object will be subset to just the samples in the new metadata.
#'        With the exception of '.sample', all character columns will be 
#'        converted to factor.
#'        
#' @return A tibble (data.frame) of the metadata in \code{biom}. '.sample' will 
#'         always be the first column. If \code{field} is given, will instead 
#'         return a named vector of the field's values.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     sample_metadata(hmp50) %>% head(4)
#'     
#'     sample_metadata(hmp50, "Body Site") %>% head(4)
#'     
#'     md <- sample_metadata(hmp50)
#'     md <- md[,c('.sample', 'Sex', 'Body Site')]
#'     sample_metadata(hmp50) <- md
#'     sample_metadata(hmp50) %>% head(4)
#'

sample_metadata <- function (biom, field=NULL) {
  
  validate_biom(clone = FALSE)
  validate_meta("field", null_ok = TRUE)
  
  
  if (is.null(field)) {
    result <- biom[['metadata']] %>% add_class('rbiom_tbl')
    
  } else {
    result <- setNames(
      object = biom[['metadata']][[field]], 
      nm     = as.character(biom[['metadata']][['.sample']]) )
  }
  
  
  return (result)
}


#' @rdname sample_metadata
#' @export

`sample_metadata<-` <- function(biom, value) {
  
  validate_biom(clone = TRUE)
  
  cmd     <- attr(value, 'cmd', exact = TRUE)
  samples <- sample_names(biom)
  
  
  #________________________________________________________
  # Parse a file/URL.
  #________________________________________________________
  if (is_scalar_character(value)) {
    
    value <- import_table(value, row.names = 1, header = TRUE)
    
    if (hasName(value, '.sample') || !any(rownames(value) %in% samples))
      stop ("The first column must have the sample IDs.")
  }
  
  
  #________________________________________________________
  # Convert data.frame/matrix to tibble.
  #________________________________________________________
  value <- tibble::as_tibble(
    x            = value, 
    rownames     = NA,
    .name_repair = trimws )

    
  #________________________________________________________
  # Make sure there's a single .sample column, at pos 1.
  #________________________________________________________
  if (sum(colnames(value) == '.sample') > 1)
    stop ("Only one column named '.sample' is allowed.")
  
  if (!(has_rownames(value) || hasName(value, '.sample')))
    stop ("Metadata must have row names or a '.sample' column.")
  
  if (has_rownames(value) && hasName(value, '.sample'))
    stop ("Row names are not allowed when a '.sample' column is present.")
  
  if (has_rownames(value)) { value %<>% rownames_to_column('.sample')
  } else                   { value %<>% relocate(.sample) }
  
  
  
  #________________________________________________________
  # Disallow duplicated column names / starting with a "."
  #________________________________________________________
  local({
    cols <- colnames(value)
    dups <- unique(cols[duplicated(cols)])
    dots <- unique(setdiff(cols[startsWith(cols, ".")], '.sample'))
    if (length(dups) > 0) stop ("Duplicated column names: ", paste(dups, collapse = ", "))
    if (length(dots) > 0) stop ("Column names can't start with a '.': ", paste(dups, collapse = ", "))
  })
  
  
  #________________________________________________________
  # Convert all character cols to factor.
  #________________________________________________________
  value %<>% mutate(across(where(is.character), as.factor))
  
  
  
  #________________________________________________________
  # Ignore IDs that aren't currently in the rbiom object.
  #________________________________________________________
  ignored <- setdiff(as.character(value[['.sample']]), samples)
  n       <- length(ignored)
  if (n > 0) {
    if (n > 4) ignored <- c(head(ignored, 4), "...")
    ignored <- paste(collapse = ", ", ignored)
    msg <- "Ignoring %i extra sample ID%s: %s"
    message(sprintf(msg, n, ifelse(n == 1, "", "s"), ignored))
  }
  
  
  
  #________________________________________________________
  # Reorder and subset to match incoming sample ids.
  #________________________________________________________
  samples <- intersect(as.character(value[['.sample']]), samples)
  biom[['counts']]   <- biom[['counts']][,samples]
  biom[['metadata']] <- left_join(
    x  = tibble(.sample = samples),
    y  = value, 
    by = ".sample" )
  
  
  biom <- biom_repair(biom)
  
  
  if (!is.null(cmd)) {
    
    if (startsWith(cmd, "%<>%")) { lhs <- "sample_metadata(biom)"
    } else                       { lhs <- "sample_metadata(biom) <- " }
    
    attr(biom, 'history') <- paste0(
      collapse = "\n", c(
        attr(biom, 'history', exact = TRUE),
        paste(lhs, cmd)))
  }
  
  
  invalidate_biom()
  return (biom)
}


#' Modify sample metadata.
#' 
#' @inherit documentation_default
#' 
#' @family metadata
#' 
#' @param expr   A single or compound expression to evaluate.
#'        
#' @return An \code{rbiom} object with the updated metadata.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     sample_metadata(biom)
#'     
#'     biom %<>% within(Sex2 <- as.numeric(Sex))
#'     sample_metadata(biom)
#'     
#'     biom %<>% within({
#'       `Body Site` <- substr(`Body Site`, 1, 4)
#'       remove("BMI", "Age", "Sex2")
#'     })
#'     sample_metadata(biom)
#'

within.rbiom <- function (biom, expr) {
  
  validate_biom(clone = FALSE)
  
  md <- biom[['metadata']]
  
  #________________________________________________________
  # NextMethod() doesn't work here.
  # Copy/pasted contents of within.data.frame:
  #________________________________________________________
  parent <- parent.frame()
  e <- evalq(environment(), md, parent)
  eval(substitute(expr), e)
  l <- as.list(e, all.names = TRUE)
  l <- l[!vapply(l, is.null, NA, USE.NAMES = FALSE)]
  nl <- names(l)
  del <- setdiff(names(md), nl)
  md[nl] <- l
  md[del] <- NULL
  sample_metadata(biom) <- md
  
  
  attr(biom, 'history') %<>% c(format(match.call())) %>% paste(collapse = "\n")
  
  invalidate_biom()
  return (biom)
}



#' Subset by sample metadata.
#' 
#' @inherit documentation_default
#' 
#' @family metadata
#' 
#' @param subset   Logical expression indicating elements or rows to keep.
#' 
#' @param select   An expression, indicating columns to select from a data frame.
#'        
#' @return An \code{rbiom} object with only the specified samples.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     sample_metadata(biom) %>% head()
#'     
#'     biom %<>% subset(Sex == "Male")
#'     sample_metadata(biom) %>% head()
#'

subset.rbiom <- function (biom, subset, select) {
  
  validate_biom(clone = FALSE)
  
  md <- biom[['metadata']]
  
  #________________________________________________________
  # NextMethod() doesn't work here.
  # Copy/pasted contents of subset.data.frame:
  #________________________________________________________
  r <- if (missing(subset)) 
    rep_len(TRUE, nrow(md))
  else {
    e <- substitute(subset)
    r <- eval(e, md, parent.frame())
    if (!is.logical(r)) 
      stop("'subset' must be logical")
    r & !is.na(r)
  }
  vars <- if (missing(select)) 
    rep_len(TRUE, ncol(md))
  else {
    nl <- as.list(seq_along(md))
    names(nl) <- names(md)
    eval(substitute(select), nl, parent.frame())
  }
  
  sample_metadata(biom) <- md[r, vars]
  
  
  attr(biom, 'history') %<>% c(format(match.call())) %>% paste(collapse = "\n")
  
  invalidate_biom()
  return (biom)
}

