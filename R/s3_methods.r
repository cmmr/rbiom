
#' Print code blocks neatly.
#' 
#' @param x   An object of class `rbiom_code`.
#' @param ...   Not used.
#' 
#' @noRd
#' @keywords internal
#' @export
print.rbiom_code <- function (x, ...) {
  
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
`$.rbiom_tbl` <- function (obj, nm) {
  
  if (!nm %in% c('cmd', 'code', 'stats', 'taxa_coords', 'taxa_stats')) {
    NextMethod()
    
  } else if (hasName(obj, nm)) {
    NextMethod()
    
  } else {
    
    val <- attr(obj, nm, exact = TRUE)
    
    if (!is.null(val) && nm %in% c('cmd', 'code'))
      return (add_class(val, 'rbiom_code'))
    
    if (is.null(val) && hasName(obj, 'data'))
      val <- attr(obj[['data']], nm, exact = TRUE)
    
    return (val)
  }
  
}


#' @export
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
#' @param ...   Not used.
#' 
#' @return A list with names
#'         `c('counts', 'metadata', 'taxonomy', 'tree', 'sequences', 'id', 'comment', 'date', 'generated_by')`.
#' 
#' @export
#' 
as.list.rbiom <- function (x, ...) {
  biom <- x
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





#' Convert an rbiom object to a simple count matrix.
#' 
#' Identical to running `as.matrix(biom$counts)`.
#' 
#' @inherit documentation_biom.rbiom
#' 
#' @family conversion
#' 
#' @param ...   Not used.
#' 
#' @return A base R matrix with OTUs as rows and samples as columns.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     as.matrix(hmp50)[1:5,1:5]
#' 
as.matrix.rbiom <- function (x, ...) {
  as.matrix(x$counts)
}




#' Map sample names to metadata field values.
#' 
#' @inherit documentation_biom.rbiom
#' 
#' @family samples
#' 
#' @param var   The metadata field name specified as:
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
#' @param ...   Not used.
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
pull.rbiom <- function (.data, var = -1, name = ".sample", ...) {
  biom <- .data
  if (!is_integerish(var))  validate_biom_field("var")
  if (!is_integerish(name)) validate_biom_field("name", null_ok = TRUE)
  dplyr::pull(.data = biom$metadata, var = var, name = name)
}




#' Get a glimpse of your metadata.
#' 
#' @inherit documentation_biom.rbiom
#' 
#' @family metadata
#' 
#' @param width   Width of output. See [pillar::glimpse()] documentation. 
#'        Default: `NULL`
#' 
#' @param ...   Not used.
#' 
#' @return The original `biom`, invisibly.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     glimpse(hmp50)
#'     
glimpse.rbiom <- function (x, width = NULL, ...) {
  biom <- x
  eval.parent(pillar::glimpse(x = biom$metadata, width = width))
  return (invisible(biom))
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
#' @name modify_metadata
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
mutate.rbiom <- function (.data, ..., clone = TRUE) {
  biom <- if (isTRUE(clone)) .data$clone() else .data
  biom$metadata <- eval.parent(dplyr::mutate(.data = biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}



#' @rdname modify_metadata
#' @export
rename.rbiom <- function (.data, ..., clone = TRUE) {
  biom <- if (isTRUE(clone)) .data$clone() else .data
  biom$metadata <- eval.parent(dplyr::rename(.data = biom$metadata, ...))
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}





#' Evaluate expressions on metadata.
#' 
#' `with()` will return the result of your expression. `within()` will return 
#' an rbiom object.
#' 
#' @inherit documentation_biom.rbiom
#' @inherit documentation_default
#' 
#' @name with
#' @family transformations
#' 
#' @param expr   Passed on to [base::with()] or [base::within()].
#' 
#' @param ...   Not used.
#' 
#' @return See description.
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     with(hmp50, table(`Body Site`, Sex))
#'     
#'     biom <- within(hmp50, {
#'       age_bin = cut(Age, 5)
#'       bmi_bin = cut(BMI, 5)
#'     })
#'     biom$metadata
#' 
with.rbiom <- function (data, expr, ...) {
  eval(expr = substitute(expr), envir = data$metadata, enclos = parent.frame())
}



#' @rdname with
#' @export
within.rbiom <- function (data, expr, clone = TRUE, ...) {
  
  biom <- if (isTRUE(clone)) data$clone() else data
  data <- biom$metadata
  
  # copied from base::within.data.frame
  parent <- parent.frame()
  e <- evalq(environment(), data, parent)
  eval(substitute(expr), e)
  l <- as.list(e, all.names = TRUE)
  l <- l[!vapply(l, is.null, NA, USE.NAMES = FALSE)]
  nl <- names(l)
  del <- setdiff(names(data), nl)
  data[nl] <- l
  data[del] <- NULL
  
  biom$metadata <- data
  
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}





#' Subset an rbiom object by sample names, OTU names, metadata, or taxonomy.
#' 
#' Dropping samples or OTUs will lead to observations being removed from the 
#' OTU matrix (`biom$counts`). OTUs and samples with zero observations are 
#' automatically removed from the rbiom object.
#'
#' @inherit documentation_biom.rbiom
#' @inherit documentation_return.biom return
#' @inherit documentation_default
#' 
#' @name subset
#' @family transformations
#' 
#' @param subset   Logical expression for rows to keep. See [base::subset()].
#' 
#' @param i,j   The sample or OTU names to keep. Or a logical/integer vector 
#'        indicating which sample names from `biom$samples` or `biom$otus` to 
#'        keep. Subsetting with `[i]` takes `i` as samples, whereas `[i,j]` 
#'        takes `i` as otus and `j` as samples (corresponding to `[rows, cols]`
#'        in the underlying `biom$counts` matrix).
#' 
#' @param fields   Which metadata field(s) to check for `NA`s, or `".all"` to
#'        check all metadata fields.
#' 
#' @param drop   Not used
#' 
#' @param ...   Not used.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     library(dplyr)
#'     
#'     # Subset to specific samples
#'     biom <- hmp50[c('HMP20', 'HMP42', 'HMP12')]
#'     biom$metadata
#'     
#'     # Subset to specific OTUs
#'     biom <- hmp50[c('LtbAci52', 'UncO2012'),] # <- Trailing ,
#'     biom$taxonomy
#'     
#'     # Subset to specific samples and OTUs
#'     biom <- hmp50[c('LtbAci52', 'UncO2012'), c('HMP20', 'HMP42', 'HMP12')]
#'     as.matrix(biom)
#'     
#'     # Subset samples according to metadata
#'     biom <- subset(hmp50, `Body Site` %in% c('Saliva') & Age < 25)
#'     biom$metadata
#'     
#'     # Subset OTUs according to taxonomy
#'     biom <- subset_taxa(hmp50, Phylum == 'Cyanobacteria')
#'     biom$taxonomy
#'     
#'     # Remove samples with NA metadata values
#'     biom <- mutate(hmp50, BS2 = na_if(`Body Site`, 'Saliva'))
#'     biom$metadata
#'     biom <- na.omit(biom)
#'     biom$metadata
#' 

subset.rbiom <- function (x, subset, clone = TRUE, ...) {
  biom <- if (isTRUE(clone)) x$clone() else x
  keep <- eval(expr = substitute(subset), envir = biom$metadata, enclos = parent.frame())
  biom$counts <- biom$counts[,keep]
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname subset
#' @export
`[.rbiom` <- function (x, i, j, ..., clone = TRUE, drop = FALSE) {
  biom        <- if (isTRUE(clone)) x$clone() else x
  biom$counts <- if (nargs() == 2) biom$counts[,i,drop] else biom$counts[i,j,drop]
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}




#' @rdname subset
#' @export
na.omit.rbiom <- function (object, fields = ".all", clone = TRUE, ...) {
  
  biom <- if (isTRUE(clone)) object$clone() else object
  
  if (length(fields) == 0) return (biom)
  
  if (eq(fields, ".all")) fields <- biom$fields
  validate_biom_field('fields', max = Inf)
  
  keep <- stats::complete.cases(biom$metadata[,fields,drop=FALSE])
  biom$counts <- biom$counts[,keep]
  
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}





#' Subset to a specific number of samples.
#'
#' @inherit documentation_biom.rbiom
#' @inherit documentation_return.biom return
#' @inherit documentation_default
#' @inheritParams dplyr::slice
#' 
#' @name slice_metadata
#' @family transformations
#' 
#' @param ...   For `slice()`, integer row indexes. For other `slice_*()` 
#'        functions, not used. See [dplyr::slice()].
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
slice.rbiom <- function (.data, ..., .by = NULL, .preserve = FALSE, clone = TRUE) {
  biom <- if (isTRUE(clone)) .data$clone() else .data
  df   <- biom$metadata
  df   <- eval.parent(dplyr::slice(df, ..., .by = {{.by}}, .preserve = .preserve))
  suppressWarnings(biom$metadata <- df)
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname slice_metadata
#' @export
slice_head.rbiom <- function (.data, n, prop, by = NULL, clone = TRUE, ...) {
  biom <- if (isTRUE(clone)) .data$clone() else .data
  df   <- biom$metadata
  df   <- eval.parent(dplyr::slice_head(df, n = n, prop = prop, by = {{by}}))
  suppressWarnings(biom$metadata <- df)
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname slice_metadata
#' @export
slice_tail.rbiom <- function (.data, n, prop, by = NULL, clone = TRUE, ...) {
  biom <- if (isTRUE(clone)) .data$clone() else .data
  df   <- biom$metadata
  df   <- eval.parent(dplyr::slice_tail(df, n = n, prop = prop, by = {{by}}))
  suppressWarnings(biom$metadata <- df)
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname slice_metadata
#' @export
slice_min.rbiom <- function (.data, order_by, n, prop, by = NULL, with_ties = TRUE, na_rm = FALSE, clone = TRUE, ...) {
  
  biom <- if (isTRUE(clone)) .data$clone() else .data
  
  df <- eval.parent(
    expr = dplyr::slice_min(
      .data     = biom$metadata,
      order_by  = {{order_by}},
      n         = n,
      prop      = prop,
      by        = {{by}},
      with_ties = with_ties,
      na_rm     = na_rm ))
  
  suppressWarnings(biom$metadata <- df)
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname slice_metadata
#' @export
slice_max.rbiom <- function (.data, order_by, n, prop, by = NULL, with_ties = TRUE, na_rm = FALSE, clone = TRUE, ...) {
  
  biom <- if (isTRUE(clone)) .data$clone() else .data
  
  df <- eval.parent(
    expr = slice_max(
      .data     = biom$metadata, 
      order_by  = {{order_by}}, 
      n         = n, 
      prop      = prop, 
      by        = {{by}}, 
      with_ties = with_ties, 
      na_rm     = na_rm ))
  
  suppressWarnings(biom$metadata <- df)
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


#' @rdname slice_metadata
#' @export
slice_sample.rbiom <- function (.data, n, prop, by = NULL, weight_by = NULL, replace = FALSE, clone = TRUE, ...) {
  biom <- if (isTRUE(clone)) .data$clone() else .data
  df   <- biom$metadata
  df   <- eval.parent(slice_sample(df, n = n, prop = prop, by = {{by}},  weight_by = {{weight_by}}, replace = replace))
  suppressWarnings(biom$metadata <- df)
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}


