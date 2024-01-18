
#' Print code blocks neatly.
#' 
#' @noRd
#' @keywords internal
#' @export
print.rbiom_code <- function (x) {
  
  if (nzchar(system.file(package = "prettycode"))) {
    
    code_style <- list(
      reserved = crayon::red,
      number   = crayon::magenta,
      null     = crayon::combine_styles(crayon::magenta, crayon::bold),
      operator = crayon::green,
      call     = crayon::cyan,
      string   = crayon::yellow,
      comment  = crayon::combine_styles(crayon::make_style("darkgrey"), crayon::italic),
      bracket  = c(crayon::yellow, crayon::magenta, crayon::cyan)
    )
    
    x <- strsplit(x, '\n')[[1]] %>%
      prettycode::highlight(style = code_style) %>%
      paste(collapse = '\n')
  }
  
  cat(x)
  
  return (invisible(NULL))
}


#' Also search attributes.
#' 
#' @noRd
#' @keywords internal
#' @export
`$.rbiom_tbl` <- function(obj, nm) {
  
  if (!nm %in% c('code', 'stats', 'taxa_coords', 'taxa_stats')) {
    NextMethod()
    
  } else if (hasName(obj, nm)) {
    NextMethod()
    
  } else {
    
    val <- attr(obj, nm, exact = TRUE)
    
    if (is.null(val) && hasName(obj, 'data'))
      val <- attr(obj[['data']], nm, exact = TRUE)
    
    return (val)
  }
  
}

`$.rbiom_plot` <- `$.rbiom_tbl`



#' Adds attr(,'tbl_sum') to tibble header.
#' 
#' @noRd
#' @keywords internal
#' @export
tbl_sum.rbiom_tbl <- function (x) {
  c(attr(x, 'tbl_sum'), NextMethod())
}






#' Convert an rbiom object to a base R list.
#' 
#' @inherit documentation_biom.rbiom
#' 
#' @family conversion
#' 
#' @return A list with names
#'         `c('counts', 'metadata', 'taxonomy', 'tree', 'sequences', 'id', 'comment', 'date', 'generated_by')`.
#' 
#' @export
#' 
as.list.rbiom <- function (biom) {
  list(
    'counts'       = biom$counts, 
    'metadata'     = biom$metadata, 
    'taxonomy'     = biom$taxonomy, 
    'tree'         = biom$tree, 
    'sequences'    = biom$sequences, 
    'id'           = biom$id, 
    'comment'      = biom$comment, 
    'date'         = biom$date, 
    'generated_by' = biom$generated_by )
}




#' Map sample names to metadata field values.
#' 
#' @inherit documentation_biom.rbiom
#' 
#' @family samples
#' 
#' @param field   The metadata field name specified as:
#' \itemize{
#'   \item{The metadata field name to retrieve. Can be abbreviated.}
#'   \item{A positive integer, giving the position counting from the left.}
#'   \item{A negative integer, giving the position counting from the right.}
#' }
#' Default: `-1`
#' 
#' @param name   The column to be used as names for a named vector. 
#'        Specified in a similar manner as var. Default: `".sample"`
#' 
#' @param ...   Passed on to [dplyr::pull()].
#' 
#' @return A vector of metadata values, named with sample names.
#' 
#' @seealso `taxa_map()`
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     pull(hmp50, 'Age') %>% head()
#'     
#'     pull(hmp50, 'bod') %>% head(4)
#'     
pull.rbiom <- function (biom, field = -1, name = ".sample", ...) {
  if (!is_integerish(field)) validate_meta("field")
  if (!is_integerish(name))  validate_meta("name", null_ok = TRUE)
  dplyr::pull(.data = biom$metadata, var = field, name = name, ...)
}




#' Create, modify, and delete metadata fields.
#' 
#' mutate() creates new fields in `$metadata` that are functions of existing 
#' metadata fields. It can also modify (if the name is the same as an existing 
#' field) and delete fields (by setting their value to NULL).
#' 
#' @inherit documentation_biom.rbiom
#' @inherit documentation_return.biom return
#' @inherit documentation_default
#' 
#' @name mutate
#' @family transformations
#' 
#' @param ...   Passed on to [dplyr::mutate()] or [dplyr::rename()].
#' 
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     biom <- slice_max(hmp50, BMI, n = 6)
#'     biom$metadata
#'     
#'     # Add a new field to the metadata
#'     biom <- mutate(biom, Obsese = BMI >= 30)
#'     biom$metadata
#'     
#'     # Rename a metadata field
#'     biom <- rename(biom, 'Age (years)' = "Age")
#'     biom$metadata
#' 
mutate.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(dplyr::mutate(.data = biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}



#' @rdname mutate
#' @export
rename.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(dplyr::rename(.data = biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}





#' Subset an rbiom object by sample names or metadata.
#'
#' @inherit documentation_biom.rbiom
#' @inherit documentation_return.biom return
#' @inherit documentation_default
#' 
#' @name subset
#' @family transformations
#' 
#' @param ...   Passed on to [base::subset()].
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     # Subset to specific samples
#'     biom <- hmp50[c('HMP20', 'HMP42', 'HMP12')]
#'     biom$metadata
#'     
#'     # Subset according to metadata
#'     biom <- subset(hmp50, `Body Site` %in% 'Saliva' & Age < 25)
#'     biom$metadata
#' 

subset.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(base::subset(x = biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname subset
#' @export
`[.rbiom` <- function(biom, i) {
  
  #________________________________________________________
  # Sanity checks.
  #________________________________________________________
  if (is.null(i) || !(is.character(i) || is_integerish(i) || is.logical(i)))
    cli_abort("Expected a character, integer, or logical vector, not {.type {i}}.")
  
  if (anyNA(i))
    cli_abort("Can't subset rbiom samples: NA in vector")
  
  if (is.character(i) && length(x <- setdiff(i, biom$samples)) > 0)
    cli_abort("Sample{?s} not in `biom`: {x}")
  
  if (is.numeric(i) && length(x <- i[i < 1 | i > biom$n_samples]) > 0)
    cli_abort("Invalid sample indices: {x}")
  
  if (is.logical(i) && length(i) != biom$n_samples)
    cli_abort("Logical vector must have {biom$n_samples} items, not {length(i)}.")
  
  
  biom <- biom$clone()
  biom$counts <- biom$counts[,i]
  return (biom)
}





#' Subset to a specific number of samples.
#'
#' @inherit documentation_biom.rbiom
#' @inherit documentation_return.biom return
#' @inherit documentation_default
#' 
#' @name slice
#' @family transformations
#' 
#' @param ...   Passed on to [dplyr::slice()].
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     # The last 3 samples in the metadata table.
#'     biom <- slice_tail(hmp50, n = 3)
#'     biom$metadata
#'     
#'     # The 3 oldest subjects sampled.
#'     biom <- slice_max(hmp50, Age, n = 3)
#'     biom$metadata
#'     
#'     # Pick 3 samples at random.
#'     biom <- slice_sample(hmp50, n = 3)
#'     biom$metadata
#' 
slice.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(dplyr::slice(.data = biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname slice
#' @export
slice_head.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(dplyr::slice_head(.data = biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname slice
#' @export
slice_tail.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(dplyr::slice_tail(.data = biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname slice
#' @export
slice_min.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(dplyr::slice_min(.data = biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname slice
#' @export
slice_max.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(slice_max(.data = biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname slice
#' @export
slice_sample.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(slice_sample(.data = biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


