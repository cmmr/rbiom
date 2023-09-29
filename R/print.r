#' Summarize the contents of a BIOM object
#' 
#' @name print
#' @param x  A BIOM object, as returned from [read_biom()].
#' @param ...  Not used.
#' @return NULL (invisibly)
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     print(hmp50)
#'

print.BIOM <- function (x, ...) {
  
  w  <- getOption("width") * 0.75
  gc <- function (vals) {
    last <- ifelse(length(vals) == 2, " and ", ", and ")
    glue_collapse(vals, sep=", ", width=w - 20, last=last)
  }
  i  <- x[['info']]
  
  cat(paste(collapse="\n", c(
    if (isTRUE(nchar(i[['id']])      > 0)) glue('{i$id} ({substr(i$date, 1, 10)})')     else NULL,
    if (isTRUE(nchar(i[['comment']]) > 0)) c(strwrap(i[['comment']], w), "-----------") else NULL,
    sprintf("%7.0f Samples:  (%s)", n_samples(x),          gc(sample_names(x))),
    sprintf("%7.0f Taxa:     (%s)", n_otus(x),             gc(otu_names(x))),
    sprintf("%7.0f Ranks:    (%s)", length(taxa_ranks(x)),    gc(taxa_ranks(x))),
    sprintf("%7.0f Metadata: (%s)", ncol(sample_metadata(x)), gc(colnames(sample_metadata(x)))),
    sprintf("        Tree:     %s", ifelse(has_tree(x),  "Present", "Absent")), 
    "\n"
  )))
  
  invisible(NULL)
}

