#' Parse counts, metadata, taxonomy, and phylogeny from a BIOM file.
#'
#' @inherit documentation_return.biom return
#'
#' @param src   Input data as either a file path, URL, or JSON string.
#'        BIOM files can be formatted according to
#'        version 1.0 (JSON) or 2.1 (HDF5)
#'        \href{http://biom-format.org/documentation/}{specifications}, or as
#'        classical tabular format. URLs must begin with \code{http://},
#'        \code{https://}, \code{ftp://}, or \code{ftps://}. JSON strings are
#'        recognized by their first non-whitespace character being \code{\{}.
#'        Compressed (gzip or bzip2) BIOM files are also supported. NOTE: to
#'        read HDF5 formatted BIOM files, the \code{h5lite} R package must be
#'        installed.
#'
#' @param ...   Properties to set in the new rbiom object, for example,
#'        `metadata`, `id`, `comment`, or `tree`. These override the
#'        corresponding content of the file.
#'
#' @details
#'
#' The file's `generated_by` and date fields are kept in the returned object's
#' `$generated_by` and `$date`.
#'
#' For classic tabular files, the first column must hold OTU IDs and a column
#' named `taxonomy` (or `lineage`) is parsed into taxonomic ranks. All other
#' columns are abundance counts: blank cells are read as `0` with a warning,
#' and columns containing non-numeric text are ignored with a warning.
#' Empty taxonomic ranks are `NA`.
#'
#' @seealso `as_rbiom()`
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read_biom(infile)
#'
#'     print(biom)
#'
#'     # Taxa Abundances
#'     biom$counts[1:4,1:10] %>% as.matrix()
#'
#'     biom$taxonomy %>% head()
#'
#'     # Metadata
#'     biom$metadata %>% head()
#'
#'     table(biom$metadata$Sex, biom$metadata$`Body Site`)
#'
#'     sprintf("Mean age: %.1f", mean(biom$metadata$Age))
#'
#'     # Phylogenetic tree (plotting requires the 'ape' package)
#'     tree <- tree_subset(biom$tree, 1:10)
#'     tree
#'     if (requireNamespace("ape", quietly = TRUE))
#'       ape::plot.phylo(tree)
#'

read_biom <- function (src, ...) {

  dots <- list(...)

  if (is.logical(dots$tree) || identical(dots$tree, 'auto')) {
    details <- "`tree` must be a newick string, phylo object, or NULL."
    lifecycle::deprecate_warn("2.0.0", "read_biom()", details = details)
    if (identical(dots$tree, FALSE)) { dots['tree'] <- list(NULL) }
    else                             { dots$tree    <- NULL       }
  }


  #________________________________________________________
  # Get the data into a file.
  #________________________________________________________
  fpath  <- as_filepath(src)
  format <- fpath$format
  fp     <- fpath$path
  on.exit(fpath$cleanup(), add = TRUE)


  #________________________________________________________
  # Process the file according to its internal format.
  #________________________________________________________

  if (format == "hdf5") {

    #___________________#
    # HDF5 file format  #
    #___________________#

    require_package('h5lite', 'to read HDF5 formatted BIOM files')

    if (!h5lite::h5_exists(fp))
      cli_abort("HDF5 file not recognized by `h5lite`: {.file {src}}")

    h5 <- read_biom_hdf5(fp)
    on.exit(try(h5$close(), silent = TRUE), add = TRUE)

    counts    <- parse_hdf5_counts(h5)
    info      <- parse_hdf5_info(h5)
    sequences <- if (!hasName(dots, 'sequences')) parse_hdf5_sequences(h5)
    taxonomy  <- if (!hasName(dots, 'taxonomy'))  parse_hdf5_taxonomy(h5)
    metadata  <- if (!hasName(dots, 'metadata'))  parse_hdf5_metadata(h5)
    phylogeny <- if (!hasName(dots, 'tree'))      parse_hdf5_tree(h5)

    info$date         <- info[['creation-date']]
    info$generated_by <- info[['generated-by']]

    remove("h5")

  } else if (format == "json") {

    #___________________#
    # JSON file format  #
    #___________________#

    json      <- read_biom_json(fp)
    counts    <- parse_json_counts(json)
    info      <- parse_json_info(json)
    sequences <- if (!hasName(dots, 'sequences')) parse_json_sequences(json)
    taxonomy  <- if (!hasName(dots, 'taxonomy'))  parse_json_taxonomy(json)
    metadata  <- if (!hasName(dots, 'metadata'))  parse_json_metadata(json)
    phylogeny <- if (!hasName(dots, 'tree'))      parse_json_tree(json)

    remove("json")

  } else if (format == "text") {

    #___________________#
    # TSV file format   #
    #___________________#

    tsv       <- read_biom_tsv(fp, src = fpath$id)
    counts    <- tsv$counts
    info      <- list(id = fpath$id, type = 'OTU table')
    taxonomy  <- if (!hasName(dots, 'taxonomy')) tsv$taxonomy
    metadata  <- NULL
    sequences <- NULL
    phylogeny <- NULL

    remove("tsv")

  } else {
    cli_abort("`src` data type not recognized: {src}")
  }


  #________________________________________________________
  # Assemble everything as an rbiom class object.
  #________________________________________________________

  args <- c(
    list(
      'counts' = counts),
    dots,
    list(
      'metadata'     = metadata,
      'taxonomy'     = taxonomy,
      'sequences'    = sequences,
      'tree'         = phylogeny,
      'id'           = info$id,
      'date'         = info$date,
      'comment'      = info$comment,
      'generated_by' = info$generated_by ))

  args <- args[!duplicated(names(args))]
  biom <- do.call(rbiom$new, args)


  #________________________________________________________
  # Provenance: record how this object was created.
  #________________________________________________________
  prov_init(biom, prov_call('read_biom', c(list(src), dots)))


  return (biom)
}



#________________________________________________________
# Classic tabular format: OTU IDs in the first column, one column per
# sample, and optionally a 'taxonomy' column.
#________________________________________________________
read_biom_tsv <- function (fp, src = fp) {

  df <- read_delim_text(fp, src = src, header_pattern = "^#\\s*(Sample|OTU)")

  if (nrow(df) == 0) cli_abort("No abundance counts in {.file {src}}.")


  #________________________________________________________
  # OTU IDs (first column) and sample names (other columns).
  #________________________________________________________
  ids <- trimws(df[[1]])

  if (length(x <- which(!nzchar(ids))) > 0)
    cli_abort("Missing OTU ID on {qty(length(x))}data line{?s} {x} of {.file {src}}.")

  if (length(x <- unique(ids[duplicated(ids)])) > 0)
    cli_abort('Duplicated {qty(length(x))}OTU ID{?s} in {.file {src}}: {.val {x}}')

  cols <- trimws(colnames(df)[-1])
  if (length(x <- unique(cols[duplicated(cols)])) > 0)
    cli_abort('Duplicated {qty(length(x))}column name{?s} in {.file {src}}: {.val {x}}')


  #________________________________________________________
  # Classify the remaining columns as counts, taxonomy, or other text.
  #________________________________________________________
  tax_names <- c('taxonomy', 'taxon', 'lineage', 'consensus lineage',
                 'consensuslineage', 'consensus_lineage', 'classification')
  num_re    <- "^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$"

  tax_col  <- NULL
  text_col <- character(0)
  mtx_cols <- list()
  n_blank  <- integer(0)

  for (j in seq_along(cols)) {

    col <- cols[[j]]
    x   <- df[[j + 1]]

    if (is.null(tax_col) && tolower(col) %in% tax_names) {
      tax_col <- x
      next
    }

    blank <- x %in% c("", "NA")
    num   <- grepl(num_re, x)
    bad   <- !num & !blank

    if (all(bad | blank) && any(bad)) {
      text_col <- c(text_col, col)
      next
    }

    if (any(bad)) {
      k <- which(bad)[[1]]
      cli_abort(c(
        'x' = "Non-numeric value {.val {x[[k]]}} for OTU {.val {ids[[k]]}} in sample column {.val {col}} of {.file {src}}.",
        'i' = "Count columns must contain only numbers (or blanks, which are read as 0)." ))
    }

    if (any(blank)) {
      n_blank[[col]] <- sum(blank)
      x[blank]       <- "0"
    }

    mtx_cols[[col]] <- as.numeric(x)
  }

  if (length(text_col) > 0)
    cli_warn("Ignoring {qty(length(text_col))}non-numeric column{?s} in {.file {src}}: {.val {text_col}}.")

  if (length(n_blank) > 0)
    cli_warn(c(
      '!' = "Read {sum(n_blank)} blank cell{?s} as 0 in {.file {src}}.",
      'i' = "Affected sample{?s}: {.val {names(n_blank)}}." ))

  if (length(mtx_cols) == 0)
    cli_abort("No abundance counts in {.file {src}}.")


  #________________________________________________________
  # Sparse counts matrix.
  #________________________________________________________
  mtx <- matrix(
    data     = unlist(mtx_cols, use.names = FALSE),
    nrow     = length(ids),
    dimnames = list(ids, names(mtx_cols)) )

  if (anyNA(mtx))   cli_abort("Count values can't be NA in {.file {src}}.")
  if (sum(mtx) == 0) cli_abort("No abundance counts in {.file {src}}.")

  counts <- as(mtx, "sparseMatrix")
  remove("mtx", "mtx_cols")


  #________________________________________________________
  # Taxonomy from the 'taxonomy' column, or from the OTU IDs.
  #________________________________________________________
  taxonomy <- parse_tsv_taxonomy(ids, tax_col)


  return (list(counts = counts, taxonomy = taxonomy))
}



read_biom_json <- function (fp) {

  json <- tryCatch(
    expr  = {
      txt <- readLines(fp, warn = FALSE, encoding = "UTF-8")
      txt <- paste(strip_bom(txt), collapse = "\n")
      jsonlite::parse_json(txt, simplifyVector = FALSE)
    },
    error = function (e) cli_abort("Can't read JSON file {.file {fp}}: {conditionMessage(e)}") )

  if (!is.list(json) || length(x <- setdiff(c('data', 'matrix_type', 'rows', 'columns'), names(json))))
    cli_abort("BIOM file is missing {.var {x}}")

  return (json)
}

read_biom_hdf5 <- function (fp) {

  require_package('h5lite', 'to read HDF5 formatted BIOM files')

  h5 <- h5lite::h5_open(file = fp)

  entries  <- h5$ls(full.names = TRUE)
  expected <- c( "/observation/matrix/indptr", "/observation/matrix/indices",
                 "/observation/ids", "/observation/matrix/data", "/sample/ids" )

  missing <- setdiff(expected, entries)
  if (length(missing) > 0) {
    try(h5$close(), silent = TRUE)
    cli_abort(sprintf("BIOM file requires: %s", paste(collapse=",", missing)))
  }

  return (h5)
}



#________________________________________________________
# The 'id' of a JSON row/column entry, or NA when absent.
#________________________________________________________
json_id <- function (x) {
  id <- unlist(x[['id']], use.names = FALSE)
  if (length(id) == 0) NA_character_ else as.character(id[[1]])
}


#________________________________________________________
# Reduce a list of JSON scalars (one per sample) to an atomic vector.
# NULL/empty values become NA; multi-valued entries are joined with "; ".
#________________________________________________________
json_scalars <- function (vals) {

  vals <- lapply(vals, function (v) {
    v <- unlist(v, use.names = FALSE)
    if (length(v) == 0) return (NA)
    if (length(v) == 1) return (v)
    paste(v, collapse = "; ")
  })

  types <- unique(vapply(vals, function (v) class(v)[[1]], ""))

  if (all(types == "logical"))                           return (as.logical(unlist(vals)))
  if (all(types %in% c("logical", "integer", "numeric"))) return (as.numeric(unlist(vals)))

  as.character(unlist(vals))
}

parse_json_metadata <- function (json) {

  cols <- json$columns
  ids  <- vapply(cols, json_id, "")

  if (length(x <- unique(ids[duplicated(ids)])) > 0)
    cli_abort("Duplicated sample ID{?s}: {.val {x}}")

  for (x in cols)
    if (as.logical(anyDuplicated(names(x$metadata))))
      cli_abort("Metadata column names are not unique.")

  fields <- unique(unlist(lapply(cols, function (x) names(x$metadata))))
  fields <- fields[nzchar(fields)]

  df <- data.frame(row.names = ids, check.names = FALSE)
  for (f in fields)
    df[[f]] <- json_scalars(lapply(cols, function (x) x[['metadata']][[f]]))

  df <- convert_col_types(df)

  return (df)
}

parse_hdf5_metadata <- function (h5) {

  ids <- h5$read("/sample/ids")
  md  <- if (h5$exists("/sample/metadata")) h5$read("/sample/metadata") else list()
  if (!is.list(md)) md <- list()

  # Restore the column order written by rbiom (HDF5 groups are unordered).
  if (!is.null(ord <- attr(md, 'column_order')))
    md <- md[c(intersect(ord, names(md)), setdiff(names(md), ord))]

  md <- lapply(md, function (x) {
    cls <- attr(x, 'rbiom_class')
    attr(x, 'rbiom_class') <- NULL
    if (identical(cls, 'logical')) x <- as.logical(x)
    if (identical(cls, 'Date'))    x <- as.Date(x, format = "%Y-%m-%d")
    if (is.integer(x))             x <- as.numeric(x) # as in the JSON reader
    return (x)
  })

  md <- md[vapply(md, function (x) is.null(dim(x)) && length(x) == length(ids), TRUE)]

  df <- data.frame(row.names = ids, check.names = FALSE)
  for (f in names(md)) df[[f]] <- md[[f]]

  df <- convert_col_types(df)

  return (df)
}



parse_json_info <- function (json) {
  info <- list()
  for (i in setdiff(names(json), c("rows", "columns", "data", "phylogeny")))
    info[[i]] <- if (is_null(unlist(json[[i]]))) NA else unlist(json[[i]])
  return (info)
}

parse_hdf5_info <- function (h5) {
  info <- list()
  for (attr_name in h5$attr_names("/")) {
    info[[attr_name]] <- h5$read(name = "/", attr = attr_name)
  }
  return (info)
}



parse_json_counts <- function (json) {

  if (!any(c('sparse', 'dense') %in% json$matrix_type))
    stop("BIOM file's matrix_type must be either 'sparse' or 'dense'")

  if (length(json$data) == 0)
    stop("BIOM file does not have any count data.")

  TaxaIDs   <- unname(sapply(json$rows,    function (x) unlist(x$id)))
  SampleIDs <- unname(sapply(json$columns, function (x) unlist(x$id)))

  if (json$matrix_type == "sparse")
    counts <- sparseMatrix(
      i        = sapply(json$data, function (x) x[[1]]) + 1,
      j        = sapply(json$data, function (x) x[[2]]) + 1,
      x        = sapply(json$data, function (x) x[[3]]),
      dims     = c(length(TaxaIDs), length(SampleIDs)),
      dimnames = list(TaxaIDs, SampleIDs) )

  if (json$matrix_type == "dense")
    counts <- Matrix(
      data     = unlist(json$data),
      byrow    = TRUE,
      nrow     = length(TaxaIDs),
      ncol     = length(SampleIDs),
      dimnames = list(TaxaIDs, SampleIDs),
      sparse   = TRUE )

  return (counts)
}

parse_hdf5_counts <- function (h5) {
  indptr <- h5$read("/observation/matrix/indptr")

  i      <- rep(seq_len(length(indptr) - 1), times = diff(indptr))
  j      <- h5$read("/observation/matrix/indices") + 1L
  x      <- h5$read("/observation/matrix/data")
  sids   <- h5$read("/sample/ids")
  oids   <- h5$read("/observation/ids")

  sparseMatrix(
    i = i, j = j, x = x,
    dims = c(length(oids), length(sids)),
    dimnames = list(oids, sids) )
}



#________________________________________________________
# Finalize a taxonomy matrix: empty ranks become NA, and a single rank
# that merely repeats the OTU IDs is dropped.
#________________________________________________________
finish_taxonomy <- function (mtx, ids) {

  if (!is.matrix(mtx)) mtx <- matrix(mtx, ncol = 1)
  storage.mode(mtx) <- "character"
  mtx[which(!nzchar(mtx))] <- NA

  rownames(mtx) <- ids
  colnames(mtx) <- default_taxa_ranks(ncol(mtx))

  if (ncol(mtx) == 1 && eq(unname(mtx[,1]), ids)) return (NULL)

  return (mtx)
}

parse_tsv_taxonomy <- function (ids, tax_col = NULL) {

  # Taxonomy column if present, otherwise try to parse the OTU IDs.
  mtx <- split_lineages(if (is.null(tax_col)) ids else tax_col)

  finish_taxonomy(mtx, ids)
}

parse_json_taxonomy <- function (json) {

  ids <- vapply(json$rows, json_id, "")

  # JSON null inside an array must stay as NA (unlist() would drop it).
  as_chr <- function (v) {
    if (is.list(v)) vapply(v, function (e) if (length(e) == 0) NA_character_ else as.character(e[[1]]), "")
    else as.character(v)
  }

  lineages <- lapply(json$rows, function (x) {

    taxaNames <- as_chr(x$metadata$taxonomy)

    # The taxa names aren't where they're supposed to be
    if (length(taxaNames) == 0) {

      # MicrobiomeDB puts the taxa string in the 'ID' field, e.g, {"metadata":null,"id":"Archaea;Euryarchaeota;...
      if (eq(json[['generated_by']], "MicrobiomeDB")) {
        taxaNames <- as_chr(x[['id']])

      # Decontam omits the 'taxonomy' name, as in {"id":"Unc01pdq","metadata":[["Bacteria","__Fusobacteriota", ...
      } else if (is_null(names(x$metadata)) && length(x$metadata) == 1) {
        taxaNames <- as_chr(x$metadata[[1]])

      } else {
        return (as_chr(x$id))
      }
    }

    if (length(taxaNames) == 1) {
      if (is.na(taxaNames) || !nzchar(taxaNames)) return (as_chr(x$id))
      return (split_lineages(taxaNames)[1,])
    }

    return (taxaNames)
  })

  # Handle instances where some taxa strings have more levels than others.
  n   <- lengths(lineages)
  m   <- max(n)
  mtx <- matrix(NA_character_, nrow = length(lineages), ncol = m)
  mtx[cbind(rep(seq_along(lineages), n), sequence(n))] <- unlist(lineages, use.names = FALSE)

  finish_taxonomy(mtx, ids)
}

parse_hdf5_taxonomy <- function (h5) {

  ids <- h5$read("/observation/ids")

  if (!h5$exists("/observation/metadata/taxonomy")) return (NULL)

  mtx <- h5$read("/observation/metadata/taxonomy")
  if (!is.matrix(mtx)) mtx <- matrix(mtx, ncol = 1)
  if (nrow(mtx) != length(ids) && ncol(mtx) == length(ids)) mtx <- t(mtx)

  if (nrow(mtx) != length(ids))
    cli_abort("HDF5 taxonomy has {nrow(mtx)} rows but there are {length(ids)} OTUs.")

  finish_taxonomy(mtx, ids)
}




parse_json_sequences <- function (json) {

  # No sequence information
  if (!any(sapply(json$rows, function (x) 'sequence' %in% names(x$metadata) )))
    return (NULL)

  ids  <- vapply(json$rows, json_id, "")
  seqs <- vapply(json$rows, function (x) {
    s <- unlist(x$metadata$sequence)
    if (length(s) == 0 || !nzchar(s[[1]])) NA_character_ else as.character(s[[1]])
  }, "")

  setNames(seqs, ids)
}

parse_hdf5_sequences <- function (h5) {
  seqs <- NULL
  if (h5$exists("/observation/metadata/sequences")) {
    seqs        <- as.character(h5$read("/observation/metadata/sequences"))
    names(seqs) <- h5$read("/observation/ids")
  }
  return (seqs)
}



parse_json_tree <- function (json) {
  tree <- as.character(unlist(json[['phylogeny']]))
  if (isTRUE(startsWith(trimws(tree), '('))) tree else NULL
}

parse_hdf5_tree <- function (h5) {
  if (h5$exists("observation/group-metadata/phylogeny")) {
    tree <- as.character(h5$read("observation/group-metadata/phylogeny"))
    if (isTRUE(startsWith(trimws(tree), '('))) return (tree)
  }
  return (NULL)
}




default_taxa_ranks <- function (n) {
  if (n >= 6 && n <= 8)
    return(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:n])
  sprintf("Level.%i", seq_len(n))
}




