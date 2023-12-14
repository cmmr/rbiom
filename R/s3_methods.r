

#' Map sample names to metadata field values.
#' 
#' @inherit documentation_biom.rbiom
#' 
#' @family samples
#' @family rarefaction
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




#' Print code blocks neatly.
#' 
#' @noRd
#' @keywords internal
#' @export
print.rbiom_code <- function (biom) {
  
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





#' Transform an rbiom object by its metadata.
#'
#' @inherit documentation_biom.rbiom
#' @inherit documentation_return.biom return
#' @inherit documentation_default
#' 
#' @name subset
#' @family transformations
#' 
#' @param ...   Passed on to [base::subset()], [dplyr::slice()], or [dplyr::mutate()].
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     subset(hmp50, `Body Site` %in% 'Saliva' & Age < 25)$metadata
#'     
#'     slice(hmp50, 1:4)$metadata
#'     
#'     mutate(hmp50, Obsese = BMI >= 30)$metadata
#' 

subset.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(base::subset(biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname subset
#' @export
slice.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(dplyr::slice(biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname subset
#' @export
mutate.rbiom <- function (biom, ..., clone = TRUE) {
  if (isTRUE(clone)) biom <- biom$clone()
  biom$metadata <- eval.parent(dplyr::mutate(biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


