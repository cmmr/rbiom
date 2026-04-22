

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
#' @param ...   Properties to overwrite in biom: `metadata`, `taxonomy`, 
#'        `tree`, etc (see details). Setting `underscores` here will pass it to 
#'        `read_tree()`.
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
as_rbiom.phyloseq <- function (biom, ...) {
  
  dots <- list(...)
  
  counts <- as(biom@otu_table@.Data, "sparseMatrix")
  if (!biom@otu_table@taxa_are_rows)
    counts <- t(counts) # nocov
  
  args <- list(
    counts    = counts, 
    metadata  = if (!is.null(biom@sam_data))  data.frame(biom@sam_data), 
    taxonomy  = if (!is.null(biom@tax_table)) data.frame(biom@tax_table), 
    sequences = if (!is.null(biom@refseq))    as.vector(biom@refseq), 
    tree      = if (!is.null(biom@phy_tree))  biom@phy_tree, 
    id        = "Imported phyloseq Data" )
  
  args <- c(dots, args)
  args <- args[!duplicated(names(args))]
  biom <- do.call(rbiom$new, args)
  
  return (biom)
}



#' @export
as_rbiom.SummarizedExperiment <- function (biom, ...) {

  dots <- list(...)
  
  require_package('SummarizedExperiment', paste("to convert", substitute(biom)))
  assay   <- getFromNamespace('assay',   'SummarizedExperiment')
  colData <- getFromNamespace('colData', 'SummarizedExperiment')
  rowData <- getFromNamespace('rowData', 'SummarizedExperiment')
  
  args <- list(
    id       = "Imported SummarizedExperiment Data",
    counts    = as(assay(biom), 'sparseMatrix'), 
    metadata  = if (!is.null(colData(biom))) data.frame(colData(biom)),
    taxonomy  = if (!is.null(rowData(biom))) data.frame(rowData(biom)) )
  
  args <- c(dots, args)
  args <- args[!duplicated(names(args))]
  biom <- do.call(rbiom$new, args)
  
  return (biom)
}



#' @export
as_rbiom.MultiAssayExperiment <- function (biom, ...) {
  as_rbiom(biom@ExperimentList[[1]], ...)
}



#' @export
as_rbiom.TreeSummarizedExperiment <- function (biom, ...) {
  
  require_package('TreeSummarizedExperiment', paste("to convert", substitute(biom)))
  assay        <- getFromNamespace('assay',            'SummarizedExperiment')
  colData      <- getFromNamespace('colData',          'SummarizedExperiment')
  rowData      <- getFromNamespace('rowData',          'SummarizedExperiment')
  rowTree      <- getFromNamespace('rowTree',      'TreeSummarizedExperiment')
  referenceSeq <- getFromNamespace('referenceSeq', 'TreeSummarizedExperiment')
  
  dots <- list(...)
  
  args <- list(
    id        = "Imported TreeSummarizedExperiment Data", 
    counts    = as(assay(biom), 'sparseMatrix'), 
    metadata  = if (!is.null(colData(biom)))      data.frame(colData(biom)),
    tree      = if (!is.null(rowTree(biom)))      rowTree(biom),
    sequences = if (!is.null(referenceSeq(biom))) as.character(referenceSeq(biom)),
    taxonomy  = if (!is.null(rowData(biom)))      data.frame(rowData(biom)) )
  
  args <- c(dots, args)
  args <- args[!duplicated(names(args))]
  biom <- do.call(rbiom$new, args)
  
  return (biom)
}



#' @export
as_rbiom.biom <- function (biom, ...) {
  
  dots <- list(...)
  
  require_package('biomformat', 'reading this object.')
  
  header               <- getFromNamespace('header',               'biomformat')
  biom_data            <- getFromNamespace('biom_data',            'biomformat')
  observation_metadata <- getFromNamespace('observation_metadata', 'biomformat')
  sample_metadata      <- getFromNamespace('sample_metadata',      'biomformat')
  
  info <- header(biom)
  
  args <- list(
    id           = info[['id']], 
    date         = info[['date']], 
    generated_by = info[['generated_by']],
    counts       = biom_data(biom), 
    metadata     = sample_metadata(biom),
    taxonomy     = observation_metadata(biom) )
  
  args <- c(dots, args)
  args <- args[!duplicated(names(args))]
  biom <- do.call(rbiom$new, args)
  
  return (biom)
}



#' @export
as_rbiom.default <- function (biom, ...) {
  
  #________________________________________________________
  # Starting from a matrix-type.
  #________________________________________________________
  if (inherits(biom, c('matrix', 'Matrix')))
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
  
  
  cli_abort(c(
    'x' = "Unable to convert {.type {biom}} to {.cls rbiom/R6}.",
    'i' = "See {.fun as_rbiom} for a list of accepted data types for `biom`." ))
}
