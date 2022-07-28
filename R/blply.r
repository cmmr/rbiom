#' Split BIOM by metadata, apply function, and return results in a list.
#' 
#' Calls \code{plyr::dlply} internally.
#' 
#' @name blply
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#' 
#' @param vars   A character vector of metadata fields. Each unique combination
#'        of values in these columns will be used to create a subsetted BIOM
#'        object to pass to \code{FUN}. Unambiguous abbreviations of metadata
#'        fields are also accepted.
#'           
#' @param FUN   The function to execute on each BIOM subset. \code{FUN} may
#'        return any object, all of which will be returned in a named list.
#'                 
#' @param ...   Additional arguments to pass on to \code{FUN}.
#'        
#' @param fast   If \code{TRUE} (the default), the subsetted BIOM objects will
#'        still contain the full taxa table and phylogenetic tree. Set 
#'        \code{fast = FALSE} to run the slow steps of subsetting these 
#'        elements as well.
#'        
#' @return A list with the function outputs.
#' 
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     blply(hmp50, "Sex", nsamples)
#'     
#'     blply(hmp50, c("Body Site", "Sex"), function (b) {
#'       ad <- adiv_table(b)[,c("Shannon", "Simpson")]
#'       apply(ad, 2L, mean)
#'     })
#'     
#'     blply(hmp50, "Body Site", function (b) {
#'       r <- range(bdiv_dist(b, "bray"))
#'       data.frame(bray.min = r[[1]], bray.max = r[[2]])
#'     })
#'     
#'     
#'
blply <- function (biom, vars, FUN, ..., fast = TRUE) {
  
  #--------------------------------------------------------------
  # Sanity checks
  #--------------------------------------------------------------
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  if (!is(FUN, 'function')) stop("Please provide a function to FUN.")
  vars <- validate_metrics(biom, vars, mode = "meta", multi = TRUE)
  
  
  #--------------------------------------------------------------
  # Pass-through for case of no variables to ply by
  #--------------------------------------------------------------
  if (length(vars) < 1)
    vars <- NULL
  
  
  #--------------------------------------------------------------
  # Run user's function on subsetted BIOM objects
  #--------------------------------------------------------------
  .data <- metadata(biom, id = ".id")
  .vars <- if (is(vars, 'quoted')) vars else ply_cols(vars)
  .fun  <- function (df, ...) {
    subBIOM <- select(biom, df[['.id']], fast = fast)
    do.call(FUN, c(list(subBIOM), list(...)))
  }
  
  plyr::dlply(.data, .vars, .fun, ...)
  
}

