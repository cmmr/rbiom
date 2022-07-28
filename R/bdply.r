#' Split BIOM by metadata, apply function, and return results in a data frame.
#' 
#' Calls \code{plyr::ddply} internally.
#' 
#' @name bdply
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#' 
#' @param vars   A character vector of metadata fields. Each unique combination
#'        of values in these columns will be used to create a subsetted BIOM
#'        object to pass to \code{FUN}. If \code{NULL}, \code{biom} will be
#'        passed to \code{FUN} unaltered. Unambiguous abbreviations of metadata
#'        fields are also accepted.
#'           
#' @param FUN   The function to execute on each BIOM subset. \code{FUN} should
#'        return a data.frame, all of which will be \code{rbind}-ed together
#'        before being returned from \code{bdply}.
#'                 
#' @param ...   Additional arguments to pass on to \code{FUN}.
#'        
#' @param fast   If \code{TRUE} (the default), the subsetted BIOM objects will
#'        still contain the full taxa table and phylogenetic tree. Set 
#'        \code{fast = FALSE} to run the slow steps of subsetting these 
#'        elements as well.
#'        
#' @return A data.frame comprising the merged outputs of \code{FUN}, along with
#'         the columns specified by \code{vars}.
#' 
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     bdply(hmp50, "Sex", nsamples)
#'     
#'     bdply(hmp50, c("Body Site", "Sex"), function (b) {
#'       ad <- adiv_table(b)[,c("Shannon", "Simpson")]
#'       apply(ad, 2L, mean)
#'     })
#'     
#'     bdply(hmp50, "Body Site", function (b) {
#'       r <- range(bdiv_dist(b, "bray"))
#'       data.frame(bray.min = r[[1]], bray.max = r[[2]])
#'     })
#'     
#'     
#'
bdply <- function (biom, vars, FUN, ..., fast = TRUE) {
  
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
    return (do.call(FUN, c(list(biom), list(...))))
  
  
  #--------------------------------------------------------------
  # Run user's function on subsetted BIOM objects
  #--------------------------------------------------------------
  .data <- metadata(biom, id = ".id")
  .vars <- if (is(vars, 'quoted')) vars else ply_cols(vars)
  .fun  <- function (df, ...) {
    subBIOM <- select(biom, df[['.id']], fast = fast)
    do.call(FUN, c(list(subBIOM), list(...)))
  }
  
  plyr::ddply(.data, .vars, .fun, ...)
  
}

