
#' Save an rbiom object to a file.
#'
#' Automatically creates directories and adds compression based on file name.
#' \describe{
#'   \item{`write_biom()` - }{ According to [BIOM format](http://biom-format.org/documentation/) specification. }
#'   \item{`write_xlsx()` - }{ Raw data and summary tables in Excel file format. See details. }
#'   \item{`write_fasta()` - }{ Sequences only in fasta format. `biom` may also be a named character vector. }
#'   \item{`write_tree()` - }{ Phylogenetic tree only in newick format. `biom` may also be a phylo object. }
#'   \item{`write_counts()`, `write_metadata()`, `write_taxonomy()` - }{ Tab-separated values. }
#' }
#'
#' @inherit documentation_default
#'
#' @param file  Path to the output file. File names ending in \code{.gz} or
#'        \code{.bz2} will be compressed accordingly.
#'        Setting `file=NULL` for `write_fasta()`, `write_tree()`, and
#'        `write_biom(format='json')` returns a string of the output which
#'        would have been written. For `write_biom(format='tab')`, `file=NULL`
#'        returns the tibble that would have been written.
#'
#' @param format  Options are \code{"tab"}, \code{"json"}, and \code{"hdf5"},
#'        corresponding to classic tabular format, BIOM format version 1.0 and
#'        biom version 2.1, respectively. NOTE: to write HDF5 formatted BIOM
#'        files, the \code{h5lite} R package must be installed. Default: `"json"`
#'
#' @param overwrite   Replace `file` if it already exists? All of the writers
#'        refuse to overwrite an existing file unless `overwrite = TRUE`.
#'        Default: `FALSE`
#'
#' @param depth   Passed on to [rarefy()]. For `write_xlsx()` only,
#'        `depth=0` disables rarefaction. Default: `NULL`
#'
#' @param seed   Random seed to use in rarefying. See [rarefy()] function
#'        for details. Must be a non-negative integer. Default: `0`
#'
#' @param quote,sep,...   Parameters passed on to [write.table()].
#'        Default: `quote=FALSE, sep="\t"`
#'
#' @return The normalized filepath that was written to (invisibly), unless
#'         `file=NULL` (see `file` argument above).
#'
#' @details
#'
#' Output is first written to a temporary file in the same directory and then
#' moved into place, so a failed write never leaves a partial file behind.
#'
#' In the tabular (`"tab"`) format, OTU IDs, sample names and taxonomy strings
#' that contain a tab, a double quote, or a line break are enclosed in double
#' quotes (with embedded quotes doubled), and missing taxonomic ranks are
#' written as empty strings. In the JSON format, numbers are written at full
#' precision and missing values are `null`. In the HDF5 format, the order of
#' metadata fields, logical fields, and `Date` fields (stored as ISO 8601
#' strings) are preserved through attributes on the datasets, which other
#' BIOM readers ignore.
#'
#' For `write_xlsx()`, `attributes(biom)` are saved as additional worksheets if
#' the attribute is a data frame, matrix, or dist -class object. An attribute
#' named 'Reads Per Step' is treated specially and merged with the usual 'Reads
#' Per Sample' tab.
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     write_tree(hmp50) %>% substr(1, 50)
#'
#'     outfile <- write_biom(hmp50, tempfile(fileext = ".biom.gz"))
#'     file.size(outfile)
#'     unlink(outfile)
#'
#'     if (requireNamespace("openxlsx", quietly = TRUE)) {
#'
#'       hmp10 <- rarefy(hmp50[1:10])
#'
#'       attr(hmp10, "Weighted UniFrac") <- bdiv_distmat(hmp10, 'weighted_unifrac')
#'       attr(hmp10, "Jaccard")          <- bdiv_distmat(hmp10, 'jaccard')
#'
#'       outfile <- write_xlsx(hmp10, tempfile(fileext = ".xlsx"))
#'       unlink(outfile)
#'     }
#'
write_biom <- function (biom, file, format = "json", overwrite = FALSE) {

  biom <- as_rbiom(biom)
  schema_check('write_biom', environment())

  format <- match.arg(tolower(format), CHOICES$biom.format)
  validate_bool('overwrite')

  #________________________________________________________
  # Allow NULL for json and tab
  #________________________________________________________
  if (format == "hdf5" || !is.null(file))
    file <- check_outfile(file, overwrite)


  switch (format,
    "hdf5" = write_biom_hdf5(biom, file),
    "json" = write_biom_json(biom, file),
    "tab"  = write_biom_tsv(biom, file) )
}



#________________________________________________________
# Validate an output path, create its directory, and refuse to
# clobber an existing file unless `overwrite` is TRUE.
#________________________________________________________
check_outfile <- function (file, overwrite = FALSE) {

  if (!is_scalar_character(file) || is_na(file) || !nzchar(file))
    cli_abort("`file` must be a file path, not {.type {file}}.")

  file <- normalizePath(file, winslash = "/", mustWork = FALSE)

  if (!dir.exists(dirname(file)))
    dir.create(dirname(file), recursive = TRUE)

  if (file.exists(file) && !isTRUE(overwrite))
    cli_abort(c(
      'x' = "Output file already exists: {.file {file}}",
      'i' = "Set `overwrite = TRUE` to replace it." ))

  return (file)
}



#________________________________________________________
# Move a finished temporary file into place.
#________________________________________________________
move_into_place <- function (tmp, file) {

  if (file.exists(file)) suppressWarnings(file.remove(file))

  if (!suppressWarnings(file.rename(tmp, file))) {
    ok <- suppressWarnings(file.copy(tmp, file, overwrite = TRUE))
    unlink(tmp)
    if (!ok) cli_abort("Can't write to {.file {file}}.")
  }

  return (invisible(file))
}



#________________________________________________________
# Quote strings for tab-separated output when they contain a tab, a
# double quote, or a line break. NA becomes an empty string.
#________________________________________________________
quote_tsv <- function (x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  i <- grepl('["\t\r\n]', x)
  x[i] <- paste0('"', gsub('"', '""', x[i], fixed = TRUE), '"')
  return (x)
}



#________________________________________________________
# "Bacteria; Firmicutes; ..." strings; missing ranks are empty.
#________________________________________________________
lineage_strings <- function (biom, sep = "; ", na = "") {
  tax <- as.matrix(biom$taxonomy[,-1,drop=FALSE])
  tax[is.na(tax)] <- na
  if (ncol(tax) == 0) return (rep("", nrow(tax)))
  unname(apply(tax, 1L, paste, collapse = sep))
}



#________________________________________________________
# BIOM Classic - Tabular
#________________________________________________________

write_biom_tsv <- function (biom, file) {


  # Sanity checks.
  #________________________________________________________
  if (any(biom$samples %in% c("#OTU ID", "Taxonomy")))
    cli_abort(c(x = "Sample names can't be '#OTU ID' or 'Taxonomy' when saving to tab format."))



  # Convert sparse matrix to tibble.
  #________________________________________________________
  tbl <- as.matrix(biom$counts) %>%
    tibble::as_tibble(rownames = "#OTU ID")



  # Add Taxonomy column when appropriate.
  #________________________________________________________
  if (ncol(biom$taxonomy) > 1)
    tbl[['Taxonomy']] <- lineage_strings(biom)


  if (is.null(file)) return (tbl)


  # Quote cells/headers that would otherwise break the format.
  #________________________________________________________
  out <- tbl
  colnames(out) <- quote_tsv(colnames(out))
  for (i in which(vapply(out, is.character, TRUE)))
    out[[i]] <- quote_tsv(out[[i]])


  write_wrapper(file, overwrite = TRUE, function (con) {
    write.table(out, con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  })


  return (invisible(file))
}



#________________________________________________________
# BIOM v1.0 - JSON
#________________________________________________________

write_biom_json <- function (biom, file) {

  metadata <- biom$metadata
  otus     <- biom$otus
  tax      <- as.matrix(biom$taxonomy[,-1,drop=FALSE])
  tax[is.na(tax)] <- ""
  seqs     <- biom$sequences
  dgT      <- as(biom$counts, 'TsparseMatrix')

  json <- toJSON(
    auto_unbox = TRUE,
    digits     = NA,
    na         = "null",
    null       = "null",
    x          = list(


      # Attributes
      #________________________________________________________
      id                  = biom$id,
      comment             = biom$comment,
      date                = biom$date,
      format              = "1.0.0",
      type                = "OTU table",
      format_url          = "http://biom-format.org",
      generated_by        = paste("rbiom", packageVersion("rbiom")),
      matrix_type         = "sparse",
      matrix_element_type = ifelse(all(biom$counts@x %% 1 == 0), "int", "float"),
      shape               = dim(biom$counts),


      # Phylogenetic Tree
      #________________________________________________________
      phylogeny = ifelse(is.null(biom$tree), "", write_tree(biom$tree)),


      # OTU IDs, Taxonomy, and Sequences
      #________________________________________________________
      rows = lapply(seq_along(otus), function (i) {
        dat <- list(taxonomy = I(unname(tax[i,])))
        if (!is.null(seqs)) dat[['sequence']] <- unname(seqs[[i]])
        list(id = otus[[i]], metadata = dat)
      }),


      # Sample IDs and Metadata
      #________________________________________________________
      columns = lapply(seq_len(nrow(metadata)), function (i) {
        row <- as.list(metadata[i,])
        list(
          id       = row[[1]],
          metadata = if (length(row) > 1) row[-1] else NULL )
      }),


      # Read counts
      #________________________________________________________
      data = mapply(
          FUN      = c,
          SIMPLIFY = FALSE,
          USE.NAMES = FALSE,
          dgT@i,
          dgT@j,
          dgT@x )

    ))


  if (is.null(file)) return (as.character(json))


  write_wrapper(file, overwrite = TRUE, function (con) {
    cat(enc2utf8(as.character(json)), file = con, sep = "")
  })

  return (invisible(file))
}




#________________________________________________________
# BIOM v2.1 - HDF5
#________________________________________________________

write_biom_hdf5 <- function (biom, file) {

  require_package('h5lite', 'to write HDF5 formatted BIOM files')


  # Write to a temporary file; move into place only on success.
  #________________________________________________________
  tmp <- tempfile(tmpdir = dirname(file), pattern = paste0(".", basename(file), "."), fileext = ".tmp")
  ok  <- FALSE
  h5  <- h5lite::h5_open(file = tmp)

  on.exit(add = TRUE, {
    try(h5$close(), silent = TRUE)
    if (!ok) unlink(tmp)
  })


  # Required top-level attributes:
  #________________________________________________________
  h5$write(attr = "creation-date",  as = "ascii", data = I(biom$date))
  h5$write(attr = "format-url",     as = "ascii", data = I("http://biom-format.org"))
  h5$write(attr = "format-version", as = "int64", data = c(2,1,0))
  h5$write(attr = "generated-by",   as = "ascii", data = I(paste("rbiom", packageVersion("rbiom"))))
  h5$write(attr = "id",                           data = I(biom$id))
  h5$write(attr = "nnz",            as = "int64", data = I(length(biom$counts@x)))
  h5$write(attr = "shape",          as = "int64", data = dim(biom$counts))
  h5$write(attr = "type",           as = "ascii", data = I("OTU table"))


  # Optional top-level attributes:
  #________________________________________________________
  if (nzchar(biom$comment))
    h5$write(attr = "comment", data = I(biom$comment))


  # Required groups:
  #________________________________________________________
  h5$create_group("observation/")
  h5$create_group("observation/matrix")
  h5$create_group("observation/metadata")
  h5$create_group("observation/group-metadata")
  h5$create_group("sample/")
  h5$create_group("sample/matrix")
  h5$create_group("sample/metadata")
  h5$create_group("sample/group-metadata")


  # Read counts by taxa (rows)
  #________________________________________________________
  dgT    <- as(biom$counts, 'TsparseMatrix')
  x      <- cbind(dgT@i, dgT@j, dgT@x)
  x      <- x[order(x[,1]),,drop=FALSE]
  indptr <- cumsum(unname(table(factor(x[,1]+1, 0:nrow(biom$counts)))))
  remove('dgT')

  h5$write(name = "observation/ids",                            data = rownames(biom$counts))
  h5$write(name = "observation/matrix/data",    as = "float64", data = x[,3])
  h5$write(name = "observation/matrix/indices", as = "int32",   data = x[,2])
  h5$write(name = "observation/matrix/indptr",  as = "int32",   data = indptr)


  # Read counts by sample (columns)
  #________________________________________________________
  x      <- x[order(x[,2]),,drop=FALSE]
  indptr <- cumsum(unname(table(factor(x[,2]+1, 0:ncol(biom$counts)))))

  h5$write(name = "sample/ids",                            data = colnames(biom$counts))
  h5$write(name = "sample/matrix/data",    as = "float64", data = x[,3])
  h5$write(name = "sample/matrix/indices", as = "int32",   data = x[,1])
  h5$write(name = "sample/matrix/indptr",  as = "int32",   data = indptr)


  # Sample Metadata
  #________________________________________________________
  if (ncol(biom$metadata) > 1) {

    fields <- names(biom$metadata)[-1]

    if (length(bad <- grep('/', fields, fixed = TRUE, value = TRUE)) > 0)
      cli_abort(c(
        'i' = "HDF5 BIOM spec prohibits slashes ({.var /}) in metadata field names.",
        'i' = "Either change the field name, or save to JSON format instead.",
        'x' = "Metadata field{?s} {.val {bad}} not encodable." ))

    # HDF5 groups are unordered; record the column order explicitly.
    h5$write(name = "sample/metadata", attr = "column_order", data = fields)

    for (field in fields) {

      name <- sprintf("sample/metadata/%s", field)
      data <- biom$metadata[[field]]
      cls  <- NULL

      if (is.factor(data)) {
        data <- as.character(data)

      } else if (is.logical(data)) {
        cls <- "logical" # stored as float64 0/1/NaN

      } else if (inherits(data, 'Date')) {
        cls  <- "Date"
        data <- format(data, "%Y-%m-%d")

      } else if (inherits(data, 'POSIXt')) {
        data <- format(data, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

      } else if (!is.numeric(data) && !is.character(data)) {
        data <- as.character(data)
      }

      h5$write(name = name, data = data)
      if (!is.null(cls)) h5$write(name = name, attr = "rbiom_class", data = cls)
    }
  }


  # Taxonomy
  #________________________________________________________
  if (ncol(biom$taxonomy) > 1) {
    name <- 'observation/metadata/taxonomy'
    data <- as.matrix(biom$taxonomy[,-1,drop=FALSE])
    data[is.na(data)] <- ""
    h5$write(name = name, data = data)
  }


  # Sequences
  #________________________________________________________
  if (!is.null(biom$sequences)) {
    name <- 'observation/metadata/sequences'
    h5$write(name = name, data = unname(biom$sequences))
  }


  # Phylogenetic Tree
  #________________________________________________________
  if (!is.null(biom$tree)) {
    name <- 'observation/group-metadata/phylogeny'
    h5$write(name = name, data = I(write_tree(biom$tree)))
    h5$write(name = name, attr = "data_type", data = "newick")
  }


  h5$close()
  move_into_place(tmp, file)
  ok <- TRUE

  return (invisible(file))
}




#' @rdname write_biom
#' @export
write_metadata <- function (biom, file, quote = FALSE, sep = "\t", overwrite = FALSE, ...) {

  write_wrapper(file, overwrite = overwrite, function (con) {
    as_rbiom(biom)$metadata %>%
      write.table(file = con, quote = quote, sep = sep, row.names = FALSE, ...)
  })
}


#' @rdname write_biom
#' @export
write_counts <- function (biom, file, quote = FALSE, sep = "\t", overwrite = FALSE, ...) {

  write_wrapper(file, overwrite = overwrite, function (con) {
    as_rbiom(biom)$counts %>%
      as.matrix() %>%
      write.table(file = con, sep = sep, quote = quote, ...)
  })
}


#' @rdname write_biom
#' @export
write_taxonomy <- function (biom, file, quote = FALSE, sep = "\t", overwrite = FALSE, ...) {

  write_wrapper(file, overwrite = overwrite, function (con) {
    as_rbiom(biom)$taxonomy %>%
      write.table(file = con, quote = quote, sep = sep, row.names = FALSE, ...)
  })
}




#' @rdname write_biom
#' @export
write_fasta <- function (biom, file = NULL, overwrite = FALSE) {

  if (is.character(biom) && !is.null(names(biom))) {
    seqs <- biom
  } else {
    seqs <- as_rbiom(biom)$sequences
  }

  if (is.null(seqs) || length(seqs) == 0)
    cli_abort(c('x' = "rbiom object does not have sequences."))

  # OTUs without a sequence (e.g. after biom_merge()) are skipped.
  seqs <- seqs[!is.na(seqs)]
  if (length(seqs) == 0)
    cli_abort(c('x' = "All sequences are NA."))

  lines <- sprintf(">%s\n%s", names(seqs), seqs)

  if (is.null(file)) return (paste0(paste(lines, collapse = "\n"), "\n"))

  write_wrapper(file, overwrite = overwrite, function (con) {
    writeLines(lines, con)
  })
}




#' @rdname write_biom
#' @export
write_tree <- function (biom, file = NULL, overwrite = FALSE) {

  tree <- if (inherits(biom, "phylo")) biom else as_rbiom(biom)$tree

  if (is.null(tree))
    cli_abort(c('x' = "rbiom object does not have a tree."))

  newick <- newick_string(tree)

  if (is_null(file)) return (newick)

  write_wrapper(file, overwrite = overwrite, function (con) {
    writeLines(text = newick, con = con)
  })
}



#________________________________________________________
# Quote Newick labels when they contain characters that are special in
# the format (parentheses, brackets, colon, semicolon, comma, quote),
# whitespace, or underscores (which unquoted readers turn into spaces).
# Embedded single quotes are doubled.
#________________________________________________________
quote_newick <- function (x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  i <- grepl("[][():;,'_[:space:]]", x)
  x[i] <- paste0("'", gsub("'", "''", x[i], fixed = TRUE), "'")
  return (x)
}



#________________________________________________________
# Serialize a phylo object to a Newick string. Iterative (no recursion)
# and O(n): tokens are emitted during a stack-based preorder traversal
# and pasted together once, so ladder-shaped trees with tens of
# thousands of tips are handled in a fraction of a second.
#________________________________________________________
newick_string <- function (tree) {

  if (!inherits(tree, 'phylo')) cli_abort("`tree` must be a phylo object, not {.type {tree}}.")

  edge   <- tree$edge
  parent <- edge[,1]
  child  <- edge[,2]
  n_tip  <- length(tree$tip.label)
  n_node <- max(edge)

  root <- setdiff(parent, child)
  if (length(root) != 1) cli_abort("Tree must have exactly one root node.")


  # Labels for tips and internal nodes.
  #________________________________________________________
  labels <- character(n_node)
  labels[seq_len(n_tip)] <- tree$tip.label
  if (!is.null(tree$node.label))
    labels[n_tip + seq_along(tree$node.label)] <- tree$node.label
  labels <- quote_newick(labels)


  # Branch length strings, indexed by child node.
  #________________________________________________________
  len <- character(n_node)
  if (!is.null(tree$edge.length)) {
    el <- tree$edge.length
    len[child] <- ifelse(is.na(el), "", paste0(":", as.character(el)))
  }


  # Per-node: child edges; pre-built tip and closing tokens.
  #________________________________________________________
  kids      <- tree_children(parent, child, n_node)
  tip_tok   <- paste0(labels, len)
  close_tok <- paste0(")", labels, len)
  close_tok[root] <- paste0(")", labels[root])


  # Stack-based preorder traversal emitting tokens.
  # Stack codes: v > 0 open node v; v < 0 close node -v; 0 comma.
  #________________________________________________________
  tokens <- character(3L * n_node)
  nt     <- 0L
  stack  <- integer(2L * n_node + 1L)
  sp     <- 1L
  stack[1L] <- root

  while (sp > 0L) {

    v  <- stack[sp]
    sp <- sp - 1L
    nt <- nt + 1L

    if (v == 0L) {
      tokens[nt] <- ","

    } else if (v < 0L) {
      tokens[nt] <- close_tok[-v]

    } else {

      ch <- kids[[v]]

      if (length(ch) == 0L) {
        tokens[nt] <- tip_tok[v]

      } else {
        tokens[nt] <- "("
        sp <- sp + 1L
        stack[sp] <- -v

        # Push children in reverse so the first child pops first,
        # with a comma marker between siblings.
        n  <- length(ch)
        pushed <- c(rbind(child[ch], c(rep(0L, n - 1L), NA)))
        pushed <- rev(pushed[!is.na(pushed)])
        stack[sp + seq_along(pushed)] <- pushed
        sp <- sp + length(pushed)
      }
    }
  }

  paste0(paste(tokens[seq_len(nt)], collapse = ""), ";")
}



#' Common code for writers.
#'
#' @noRd
#' @keywords internal
#'
#' @param outfile   The file to open a connection to.
#'
#' @param callback   A function that takes a connection as its only
#'        argument and writes content to it.
#'
#' @param overwrite   Allow replacing an existing file.
#'
#' @return The normalized filepath that was written to (invisibly).
#'
write_wrapper <- function (outfile, callback, overwrite = FALSE) {

  outfile <- check_outfile(outfile, overwrite)


  # Write to a temporary file in the same directory; move into place
  # only if the callback succeeds. The connection is always closed.
  #________________________________________________________
  tmp <- tempfile(tmpdir = dirname(outfile), pattern = paste0(".", basename(outfile), "."), fileext = ".tmp")
  con <- NULL
  ok  <- FALSE

  on.exit(add = TRUE, {
    if (!is.null(con)) try(close(con), silent = TRUE)
    if (!ok) unlink(tmp)
  })

  tryCatch(
    error = function (e) cli_abort(c(
      'x' = "Can't write to {.file {outfile}}.",
      'i' = conditionMessage(e) )),
    expr  = {

      # Binary mode: UTF-8 bytes as-is and LF line endings on every platform.
      lc  <- tolower(outfile)
      con <- if (endsWith(lc, ".gz"))  { gzfile(tmp, "wb")
      } else if (endsWith(lc, ".bz2")) { bzfile(tmp, "wb")
      } else                           { file(tmp, "wb") }

      callback(con)

      close(con)
      con <- NULL
    })

  move_into_place(tmp, outfile)
  ok <- TRUE

  return (invisible(outfile))
}




#' Export data to QIIME 2 or mothur.
#'
#' @description
#' Populates a directory with the following files, formatted according to
#' QIIME 2 or mothur's specifications.
#'
#' * `biom_counts.tsv`
#' * `biom_metadata.tsv`
#' * `biom_taxonomy.tsv`
#' * `biom_tree.nwk`
#' * `biom_seqs.fna`
#'
#' `biom_counts.tsv` will always be created. The others are dependent on
#' whether the content is present in the `biom` argument.
#'
#' For QIIME 2, metadata columns are typed `numeric` or `categorical`
#' (logical and `Date` fields are categorical; dates are written as
#' `YYYY-MM-DD`). For mothur, taxonomy rows use the OTU IDs from the counts
#' file, and missing ranks are written as `unclassified`.
#'
#'
#' @name export
#' @inherit documentation_default
#'
#' @param dir  Where to save the files. If the directory doesn't exist, it will
#'        be created. Default: `tempfile()`
#'
#' @param prefix  A string to prepend to each file name. Default: `'biom_'`
#'
#' @param overwrite   Replace files that already exist? Default: `FALSE`
#'
#' @return The normalized directory path that was written to (invisibly).
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     tdir <- tempfile()
#'
#'     write_qiime2(hmp50, tdir, 'qiime2_')
#'     write_mothur(hmp50, tdir, 'mothur_')
#'
#'     list.files(tdir)
#'
#'     readLines(file.path(tdir, 'qiime2_metadata.tsv'), n = 4)
#'
#'     readLines(file.path(tdir, 'mothur_taxonomy.tsv'), n = 3)
#'
#'     unlink(tdir, recursive = TRUE)
#'

write_mothur <- function (biom, dir = tempfile(), prefix = 'biom_', overwrite = FALSE) {

  biom <- as_rbiom(biom)

  dir <- normalizePath(dir, winslash = '/', mustWork = FALSE)
  if (!dir.exists(dir)) dir.create(path = dir, recursive = TRUE)

  validate_string('prefix')
  validate_bool('overwrite')


  # Counts
  write_wrapper(file.path(dir, paste0(prefix, 'counts.tsv')), overwrite = overwrite, function (con) {
    biom$counts %>%
      as.matrix() %>%
      tibble::as_tibble(rownames = "Representative_Sequence") %>%
      dplyr::mutate('total' = as.vector(rowSums(biom$counts)), .after = 1) %>%
      write.table(file=con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  })


  # Tree
  if (!is.null(biom$tree))
    write_tree(biom, file.path(dir, paste0(prefix, 'tree.nwk')), overwrite = overwrite)


  # Sequences
  if (!is.null(biom$sequences))
    write_fasta(biom, file.path(dir, paste0(prefix, 'seqs.fna')), overwrite = overwrite)


  # Taxonomy
  if (ncol(biom$taxonomy) > 1)
    write_wrapper(file.path(dir, paste0(prefix, 'taxonomy.tsv')), overwrite = overwrite, function (con) {
      tibble::tibble(
        'OTU'      = biom$otus,
        'Size'     = rowSums(biom$counts) %>% as.vector(),
        'Taxonomy' = lineage_strings(biom, sep = ';', na = 'unclassified') %>%
          gsub(' ', '_', ., fixed = TRUE) %>%
          paste0(';') ) %>%
        write.table(file=con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    })


  # Metadata
  if (ncol(biom$metadata) > 1)
    write_wrapper(file.path(dir, paste0(prefix, 'metadata.tsv')), overwrite = overwrite, function (con) {

      tbl <- biom$metadata

      colnames(tbl) %<>% gsub('[\\s\\t\\n]+', '_', ., perl = TRUE)
      colnames(tbl)[[1]] <- 'group'

      for (i in seq_len(ncol(tbl))) {
        x <- tbl[[i]]
        if (inherits(x, 'Date')) x <- format(x, "%Y-%m-%d")
        if (!is.numeric(x)) {
          x <- as.character(x)
          x[grep('[\\s\\t\\"]', x, perl = TRUE)] %<>% shQuote(type = 'cmd')
        }
        tbl[[i]] <- x
      }

      write.table(tbl, file = con, quote = FALSE, sep = "\t", row.names = FALSE, na = "")
    })


  return (invisible(dir))
}



#' @rdname export
#' @export

write_qiime2 <- function (biom, dir = tempfile(), prefix = 'biom_', overwrite = FALSE) {

  biom <- as_rbiom(biom)

  dir <- normalizePath(dir, winslash = '/', mustWork = FALSE)
  if (!dir.exists(dir)) dir.create(path = dir, recursive = TRUE)

  validate_string('prefix')
  validate_bool('overwrite')


  # Counts
  write_wrapper(file.path(dir, paste0(prefix, 'counts.tsv')), overwrite = overwrite, function (con) {
    biom$counts %>%
      as.matrix() %>%
      tibble::as_tibble(rownames = "#OTU ID") %>%
      write.table(file=con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  })


  # Tree
  if (!is.null(biom$tree))
    write_tree(biom, file.path(dir, paste0(prefix, 'tree.nwk')), overwrite = overwrite)


  # Sequences
  if (!is.null(biom$sequences))
    write_fasta(biom, file.path(dir, paste0(prefix, 'seqs.fna')), overwrite = overwrite)


  # Taxonomy
  if (ncol(biom$taxonomy) > 1)
    write_wrapper(file.path(dir, paste0(prefix, 'taxonomy.tsv')), overwrite = overwrite, function (con) {
      tibble::tibble(
        'Feature ID' = biom$otus,
        'Taxon'      = lineage_strings(biom, sep = '; ') %>% gsub(' ', '_', ., fixed = TRUE),
        'Confidence' = rep(1, biom$n_otus) ) %>%
        write.table(file=con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    })


  # Metadata
  if (ncol(biom$metadata) > 1)
    write_wrapper(file.path(dir, paste0(prefix, 'metadata.tsv')), overwrite = overwrite, function (con) {

      tbl     <- biom$metadata
      headers <- colnames(tbl)

      # As per https://docs.qiime2.org/2024.10/tutorials/metadata
      bad <- tolower(headers) %in% c('id', 'sampleid', 'sample id', 'sample-id', 'featureid', 'feature id', 'feature-id')
      bad <- bad | (headers %in% c('#SampleID', '#Sample ID', '#OTUID', '#OTU ID', 'sample_name'))
      headers <- ifelse(bad, paste0('_', headers), headers)

      headers[grep('[\\s\\t\\n\\"]', headers, perl = TRUE)] %<>% shQuote(type = 'cmd')
      headers %<>% gsub('\\"', '""', ., fixed = TRUE)

      headers[[1]] <- 'sample-id'
      cat(file = con, sep = '', paste(collapse = '\t', headers), '\n')
      cat(file = con, sep = '', '#q2:types')

      for (i in seq_len(ncol(tbl))[-1]) {

        x <- tbl[[i]]

        if (is.numeric(x)) {
          cat(file = con, sep = '', '\tnumeric')

        } else {
          cat(file = con, sep = '', '\tcategorical')
          if (inherits(x, 'Date')) x <- format(x, "%Y-%m-%d")
          x <- as.character(x)
          x[grep('[\\s\\t\\n\\"]', x, perl = TRUE)] %<>% shQuote(type = 'cmd')
          x %<>% gsub('\\"', '""', ., fixed = TRUE)
        }

        tbl[[i]] <- x
      }

      cat(file = con, sep = '', '\n')

      write.table(tbl, file = con, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, na = "")
    })


  return (invisible(dir))
}
