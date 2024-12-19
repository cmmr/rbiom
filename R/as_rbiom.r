

#' Convert a variety of data types to an rbiom object.
#' 
#' Construct an rbiom object. The returned object is an R6 reference class. 
#' Use `b <- a$clone()` to create copies, not `b <- a`.
#' 
#' @inherit documentation_return.biom return
#' 
#' @param biom   Object which can be coerced to an rbiom-class object.
#'        For example:
#'        \describe{
#'          \item{\emph{file} - }{ Filepath or URL to a biom file. }
#'          \item{\emph{matrix} - }{ An abundance matrix with OTUs in rows and samples in columns. }
#'          \item{`phyloseq`-class object - }{ From the phyloseq Bioconductor R package. }
#'          \item{\emph{list} - }{ With `counts` and optionally `metadata`, `taxonomy`, `tree`, etc (see details). }
#'        }
#' 
#' @param ...   Properties to overwrite in biom: `metadata`, `taxonomy`, `tree`, etc (see details).
#'
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     # create a simple matrix ------------------------
#'     mtx <- matrix(
#'       data     = floor(runif(24) * 1000), 
#'       nrow     = 6, 
#'       dimnames = list(paste0("OTU", 1:6), paste0("Sample", 1:4)) )
#'     mtx
#'     
#'     # and some sample metadata ----------------------
#'     df <- data.frame(
#'       .sample   = paste0("Sample", 1:4),
#'       treatment = c("A", "B", "A", "B"),
#'       days      = c(12, 3, 7, 8) )
#'     
#'     # convert data set to rbiom ---------------------
#'     biom <- as_rbiom(mtx, metadata = df, id = "My BIOM")
#'     biom
#' 
as_rbiom <- function (biom, ...) {
  UseMethod("as_rbiom")
}



#' @export
as_rbiom.rbiom <- function (biom, ...) {
  
  #________________________________________________________
  # Already an rbiom/R6 object.
  #________________________________________________________
  if (!identical(biom$pkg_version, packageVersion("rbiom"))) {
    
    args <- c(list(...), as.list(biom))
    args <- args[!duplicated(names(args))]
    biom <- do.call(rbiom$new, args)
  }
  
  return (biom)
}



#' @export
as_rbiom.default <- function (biom, ...) {
  
  #________________________________________________________
  # Starting from a matrix-type.
  #________________________________________________________
  if (inherits(biom, c('matrix', 'simple_triplet_matrix')))
    return (rbiom$new(counts = biom, ...))
  
  
  #________________________________________________________
  # A list with 'counts' and possibly other elements.
  #________________________________________________________
  if (is.list(biom) && hasName(biom, "counts")) {
    args <- c(list(counts = biom$counts), list(...), biom)
    args <- args[!duplicated(names(args))]
    biom <- do.call(rbiom$new, args)
    return (biom)
  }
  
  
  #________________________________________________________
  # Read new rbiom object from filename / URL / JSON.
  #________________________________________________________
  if (is_scalar_character(biom) && !is.na(biom))
    return (read_biom(src = biom, ...))
  
  
  #________________________________________________________
  # Convert from phyloseq.
  #________________________________________________________
  if (inherits(biom, 'phyloseq')) {
    
    dots <- list(...)
    
    args <- list(
      counts    = as.simple_triplet_matrix(biom@otu_table), 
      metadata  = if (!is.null(biom@sam_data))  data.frame(biom@sam_data), 
      taxonomy  = if (!is.null(biom@tax_table)) data.frame(biom@tax_table), 
      sequences = if (!is.null(biom@refseq))    as.vector(biom@refseq), 
      tree      = if (!is.null(biom@phy_tree))  biom@phy_tree, 
      id        = "Imported PhyloSeq Data" )
    
    args <- c(dots, args)
    args <- args[!duplicated(names(args))]
    biom <- do.call(rbiom$new, args)
    
    return (biom)
  }
  
  
  cli_abort(c(
    'x' = "Unable to convert {.type {biom}} to {.cls rbiom/R6}.",
    'i' = "See {.fun as_rbiom} for a list of accepted data types for `biom`." ))
}




