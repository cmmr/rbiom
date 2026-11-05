

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
#'          \item{\emph{file} - }{ Filepath or URL to a BIOM file (JSON, HDF5, or tabular), or a JSON string. }
#'          \item{\emph{matrix} - }{ An abundance matrix (base or `Matrix`) with OTUs in rows and samples in columns. }
#'          \item{\emph{data.frame} - }{ Treated as an abundance matrix: OTUs in rows (row names) and samples in columns. }
#'          \item{\emph{list} - }{ With `counts` and optionally `metadata`, `taxonomy`, `tree`, etc (see details). }
#'          \item{`rbiom` - }{ Returned as-is, unless `...` requests changes, in which case a modified clone is returned. }
#'          \item{`phyloseq` - }{ From the phyloseq Bioconductor R package. }
#'          \item{`SummarizedExperiment` - }{ From the SummarizedExperiment Bioconductor R package. }
#'          \item{`TreeSummarizedExperiment` - }{ From the TreeSummarizedExperiment Bioconductor R package. }
#'          \item{`MultiAssayExperiment` - }{ From the MultiAssayExperiment Bioconductor R package (its first experiment is used). }
#'          \item{`biom` - }{ From the biomformat Bioconductor R package. }
#'        }
#'
#' @param ...   Properties to overwrite in biom: `metadata`, `taxonomy`,
#'        `tree`, `sequences`, `id`, `comment`, or `date` (see details).
#'        Setting `underscores` here will pass it to `read_tree()`.
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

  dots <- list(...)

  #________________________________________________________
  # Objects from an older rbiom version are rebuilt.
  #________________________________________________________
  if (!identical(biom$pkg_version, packageVersion("rbiom"))) {

    args    <- c(dots, as.list(biom))
    args    <- args[!duplicated(names(args))]
    history <- tryCatch(prov_get(biom), error = function (e) NULL)
    biom    <- do.call(rbiom$new, args)
    if (!is.null(history)) prov_set(biom, history)

    return (biom)
  }


  #________________________________________________________
  # Apply any requested changes to a clone.
  #________________________________________________________
  if (length(dots) > 0) {

    fields <- c('counts', 'metadata', 'taxonomy', 'tree', 'sequences', 'id', 'comment', 'date')

    if (length(bad <- setdiff(names(dots), fields)) > 0)
      cli_abort(c(
        'x' = "Can't set {.val {bad}} on an rbiom object.",
        'i' = "Settable fields: {.val {fields}}." ))

    biom <- biom$clone()
    for (field in names(dots)) biom[[field]] <- dots[[field]]
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
    metadata  = if (!is.null(biom@sam_data))  data.frame(biom@sam_data, check.names = FALSE, stringsAsFactors = FALSE),
    taxonomy  = if (!is.null(biom@tax_table)) as(biom@tax_table, "matrix"),
    sequences = if (!is.null(biom@refseq))    as.character(biom@refseq),
    tree      = if (!is.null(biom@phy_tree))  biom@phy_tree,
    id        = "Imported phyloseq Data" )

  args <- c(dots, args)
  args <- args[!duplicated(names(args))]
  obj  <- do.call(rbiom$new, args)
  prov_init(obj, prov_call('as_rbiom', c(list(biom), dots)))

  return (obj)
}



#' @export
as_rbiom.SummarizedExperiment <- function (biom, ...) {

  dots <- list(...)

  require_package('SummarizedExperiment', 'to convert SummarizedExperiment objects')

  args <- list(
    id       = "Imported SummarizedExperiment Data",
    counts    = as(SummarizedExperiment::assay(biom), 'sparseMatrix'),
    metadata  = se_data_frame(SummarizedExperiment::colData(biom)),
    taxonomy  = se_data_frame(SummarizedExperiment::rowData(biom)) )

  args <- c(dots, args)
  args <- args[!duplicated(names(args))]
  obj  <- do.call(rbiom$new, args)
  prov_init(obj, prov_call('as_rbiom', c(list(biom), dots)))

  return (obj)
}



#' @export
as_rbiom.MultiAssayExperiment <- function (biom, ...) {
  require_package('MultiAssayExperiment', 'to convert MultiAssayExperiment objects')
  as_rbiom(biom@ExperimentList[[1]], ...)
}



#' @export
as_rbiom.TreeSummarizedExperiment <- function (biom, ...) {

  require_package('SummarizedExperiment',     'to convert TreeSummarizedExperiment objects')
  require_package('TreeSummarizedExperiment', 'to convert TreeSummarizedExperiment objects')

  dots <- list(...)

  seqs <- TreeSummarizedExperiment::referenceSeq(biom)

  args <- list(
    id        = "Imported TreeSummarizedExperiment Data",
    counts    = as(SummarizedExperiment::assay(biom), 'sparseMatrix'),
    metadata  = se_data_frame(SummarizedExperiment::colData(biom)),
    tree      = TreeSummarizedExperiment::rowTree(biom),
    sequences = if (!is.null(seqs)) as.character(seqs),
    taxonomy  = se_data_frame(SummarizedExperiment::rowData(biom)) )

  args <- c(dots, args)
  args <- args[!duplicated(names(args))]
  obj  <- do.call(rbiom$new, args)
  prov_init(obj, prov_call('as_rbiom', c(list(biom), dots)))

  return (obj)
}



#________________________________________________________
# A base data.frame (row names kept, column names untouched) from a
# Bioconductor DataFrame, or NULL when there are no columns.
#________________________________________________________
se_data_frame <- function (x) {
  if (is.null(x) || ncol(x) == 0) return (NULL)
  df <- as.data.frame(x, optional = TRUE, stringsAsFactors = FALSE)
  if (!identical(colnames(df), colnames(x))) colnames(df) <- colnames(x)
  if (!is.null(rownames(x))) rownames(df) <- rownames(x)
  return (df)
}



#' @export
as_rbiom.biom <- function (biom, ...) {

  dots <- list(...)

  require_package('biomformat', 'to convert biom objects')

  info     <- biomformat::header(biom)
  metadata <- biomformat::sample_metadata(biom)
  taxonomy <- biomformat::observation_metadata(biom)

  if (!is.data.frame(metadata) || nrow(metadata) == 0 || ncol(metadata) == 0)
    metadata <- NULL

  if (is.list(taxonomy) && !is.data.frame(taxonomy)) {
    # One vector of rank names per OTU; pad to a rectangular matrix.
    n   <- lengths(taxonomy)
    mtx <- matrix(NA_character_, nrow = length(taxonomy), ncol = max(n, 0L))
    mtx[cbind(rep(seq_along(taxonomy), n), sequence(n))] <- as.character(unlist(taxonomy))
    rownames(mtx) <- names(taxonomy)
    colnames(mtx) <- default_taxa_ranks(ncol(mtx))
    taxonomy <- mtx
  }

  if (!(is.data.frame(taxonomy) || is.matrix(taxonomy)) || nrow(taxonomy) == 0 || ncol(taxonomy) == 0)
    taxonomy <- NULL

  args <- list(
    id           = info[['id']],
    date         = info[['date']],
    generated_by = info[['generated_by']],
    counts       = as(biomformat::biom_data(biom), 'sparseMatrix'),
    metadata     = metadata,
    taxonomy     = taxonomy )

  args <- c(dots, args)
  args <- args[!duplicated(names(args))]
  obj  <- do.call(rbiom$new, args)
  prov_init(obj, prov_call('as_rbiom', c(list(biom), dots)))

  return (obj)
}



#' @export
as_rbiom.default <- function (biom, ...) {

  dots <- list(...)

  #________________________________________________________
  # Starting from a matrix-type (or data.frame of counts).
  #________________________________________________________
  if (inherits(biom, c('matrix', 'Matrix')) || is.data.frame(biom)) {
    obj <- rbiom$new(counts = if (is.data.frame(biom)) as.matrix(biom) else biom, ...)
    prov_init(obj, prov_call('as_rbiom', c(list(biom), dots)))
    return (obj)
  }


  #________________________________________________________
  # A list with 'counts' and possibly other elements.
  #________________________________________________________
  if (is.list(biom) && hasName(biom, "counts")) {
    args <- c(list(counts = biom$counts), dots, biom)
    args <- args[!duplicated(names(args))]
    obj  <- do.call(rbiom$new, args)
    prov_init(obj, prov_call('as_rbiom', c(list(biom), dots)))
    return (obj)
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
