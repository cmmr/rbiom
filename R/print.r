#' Summarize the contents of an rbiom object
#' 
#' @noRd
#' @keywords internal
#' @export
print.rbiom <- function (x, ...) {
  
  stopifnot(is(x, 'rbiom'))
  
  w  <- getOption("width") * 0.75
  gc <- function (vals) {
    last <- ifelse(length(vals) == 2, " and ", ", and ")
    glue_collapse(vals, sep=", ", width=w - 20, last=last)
  }
  i  <- x[['info']]
  
  cat(paste(collapse="\n", c(
    if (isTRUE(nchar(i[['id']])      > 0)) glue('{i$id} ({substr(i$date, 1, 10)})')     else NULL,
    if (isTRUE(nchar(i[['comment']]) > 0)) c(strwrap(i[['comment']], w), "-----------") else NULL,
    sprintf("%7.0f Samples:  (%s)", n_samples(x),  gc(sample_names(x))),
    sprintf("%7.0f Taxa:     (%s)", n_otus(x),     gc(otu_names(x))),
    sprintf("%7.0f Ranks:    (%s)", n_ranks(x),    gc(taxa_ranks(x))),
    sprintf("%7.0f Metadata: (%s)", n_metadata(x), gc(metadata_names(x))),
    sprintf("        Tree:     %s", ifelse(has_tree(x),  "Present", "Absent")), 
    "\n"
  )))
  
  invisible(NULL)
}



#' Print code blocks neatly.
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
`$.rbiom` <- function(obj, nm) {
  
  if (!nm %in% c('code', 'history', 'stats', 'taxa_coords', 'taxa_stats')) {
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

`$.rbiom_plot` <- `$.rbiom`
`$.rbiom_tbl`  <- `$.rbiom`



#' Adds attr(,'tbl_sum') to tibble header.
#' 
#' @noRd
#' @keywords internal
#' @export
tbl_sum.rbiom_tbl <- function (x) {
  c(attr(x, 'tbl_sum'), NextMethod())
}


