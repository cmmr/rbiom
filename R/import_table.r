
#' Read a data.frame from a file/URL/string. Auto-detects field separator.
#'
#' @noRd
#' @keywords internal
#'
#' @param src  A file name, URL, or literal data string.
#'        File types supported: csv, tsv, json, xls, and xlsx.
#'
#' @return A tibble data.frame. For delimited text, every column is character
#'         and the attribute `rbiom_untyped` is `TRUE`; the caller decides
#'         which columns to type-convert (see `convert_col_types()`), so that
#'         identifier columns such as `"1.10"` survive unchanged.
#'

import_table <- function (src) {

  #________________________________________________________
  # Get the data into a file.
  #________________________________________________________
  fpath  <- as_filepath(src)
  format <- fpath$format
  path   <- fpath$path
  on.exit(fpath$cleanup(), add = TRUE)


  #________________________________________________________
  # Read json/text/excel data into a data frame.
  #________________________________________________________
  untyped <- FALSE

  tbl <- if (format == 'excel') {
    require_package('readxl', 'to read Excel files')
    readxl::read_excel(path = path, .name_repair = "minimal")

  } else if (format == 'json') {
    as_tibble(jsonlite::read_json(path = path, simplifyVector = TRUE))

  } else if (format == 'text') {
    untyped <- TRUE
    df <- read_delim_text(path, src = fpath$id)
    for (i in seq_along(df)) {
      x <- df[[i]]
      x[x %in% c("", "NA")] <- NA
      df[[i]] <- x
    }
    as_tibble(df, .name_repair = "minimal")

  } else {
    cli_abort("`src` data type not recognized: {src}")
  }

  colnames(tbl) <- trimws(colnames(tbl))


  #________________________________________________________
  # Enforce unique column names.
  #________________________________________________________

  if (any(x <- duplicated(colnames(tbl)))) {
    x <- unique(colnames(tbl)[x])
    cli_abort("Duplicated column names in {.file {fpath$id}}: {.val {x}}")
  }

  attr(tbl, 'rbiom_untyped') <- untyped

  return (tbl)
}



#' Parse a delimited text file into a character data.frame.
#'
#' Handles tab/comma/semicolon separators (sniffed from the header line), a
#' UTF-8 byte-order mark, QIIME 2 `#q2:types` lines, `#` characters inside
#' values, apostrophes and backslashes, double-quoted fields (with `""` as
#' the embedded quote), a trailing delimiter, and a header with one fewer
#' field than the data (row names).
#'
#' @noRd
#' @keywords internal
#'
#' @param path   Path to a (possibly gzip/bzip2 compressed) text file.
#' @param src    How to refer to the file in messages.
#' @param header_pattern   Optional regex. The first line matching it is used
#'        as the header, and preceding lines are discarded.
#'
#' @return A data.frame of character columns. Empty cells are `""`.
#'
read_delim_text <- function (path, src = path, header_pattern = NULL) {

  lines <- tryCatch(
    expr  = readLines(path, warn = FALSE, encoding = "UTF-8"),
    error = function (e) cli_abort("Can't read {.file {src}}: {conditionMessage(e)}") )

  lines <- strip_bom(lines)
  lines <- sub("\r$", "", lines)
  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[!startsWith(lines, "#q2:")] # QIIME 2 column types

  if (length(lines) == 0)
    cli_abort("No data found in {.file {src}}.")


  #________________________________________________________
  # Locate the header line.
  #________________________________________________________
  if (!is.null(header_pattern)) {
    i <- grep(header_pattern, lines)
    if (length(i) == 0) i <- which(!startsWith(lines, "#"))
    if (length(i) > 0)  lines <- lines[i[[1]]:length(lines)]
  }

  sep   <- sniff_delim(lines[[1]])
  quote <- if (any(grepl(paste0('(^|', sep, ')\\s*"'), lines))) '"' else ''

  header <- scan(
    text = lines[[1]], what = "", sep = sep, quote = quote,
    strip.white = TRUE, quiet = TRUE, na.strings = character(0) )

  body <- lines[-1]


  #________________________________________________________
  # No data rows.
  #________________________________________________________
  if (length(body) == 0) {
    df <- as.data.frame(
      x                = rep(list(character(0)), length(header)),
      col.names        = header,
      check.names      = FALSE,
      stringsAsFactors = FALSE )
    return (df)
  }


  #________________________________________________________
  # Reconcile the number of fields in the header and the body.
  #________________________________________________________
  nf <- utils::count.fields(
    file = textConnection(body), sep = sep, quote = quote,
    comment.char = "", blank.lines.skip = TRUE )

  n_max <- max(nf, na.rm = TRUE)

  if (n_max == length(header) + 1) {
    header <- c("", header) # Row names

  } else if (n_max > length(header)) {
    i <- which(nf == n_max)[[1]]
    cli_abort(c(
      'x' = "Can't parse {.file {src}}.",
      'i' = "The header has {length(header)} field{?s}, but data line {i} has {n_max}." ))
  }


  #________________________________________________________
  # Parse the data lines.
  #________________________________________________________
  df <- tryCatch(
    error = function (e) cli_abort("Can't parse {.file {src}}: {conditionMessage(e)}"),
    expr  = utils::read.table(
      text             = body,
      sep              = sep,
      quote            = quote,
      header           = FALSE,
      col.names        = paste0("V", seq_along(header)),
      colClasses       = "character",
      check.names      = FALSE,
      comment.char     = "",
      allowEscapes     = FALSE,
      strip.white      = TRUE,
      fill             = TRUE,
      na.strings       = character(0),
      blank.lines.skip = TRUE,
      encoding         = "UTF-8",
      stringsAsFactors = FALSE ))

  names(df) <- header


  #________________________________________________________
  # A trailing delimiter creates an empty, unnamed column.
  #________________________________________________________
  while (ncol(df) > 1 && !nzchar(names(df)[[ncol(df)]]) && all(!nzchar(df[[ncol(df)]])))
    df[[ncol(df)]] <- NULL


  return (df)
}



#' Pick the field separator from a header line: the most frequent of tab,
#' comma, and semicolon (outside double quotes). Ties favor tab, then comma.
#'
#' @noRd
#' @keywords internal
sniff_delim <- function (line) {
  line <- gsub('"[^"]*"', '', line)
  n <- vapply(c("\t", ",", ";"), function (s) lengths(regmatches(line, gregexpr(s, line, fixed = TRUE))), 1L)
  if (all(n == 0)) return ("\t")
  names(n)[[which.max(n)]]
}



#' Remove a UTF-8 byte-order mark from the start of a string.
#'
#' @noRd
#' @keywords internal
strip_bom <- function (x) {
  if (length(x) > 0 && !is.na(x[[1]]))
    x[[1]] <- sub("^\uFEFF", "", x[[1]])
  return (x)
}



#' Convert character columns to numeric, logical, or Date where appropriate.
#'
#' Applied to tables read from text files, where every column starts as
#' character. Numbers are recognized by a regular expression that rejects
#' leading zeros (so that zero-padded identifiers stay character), logicals
#' must be spelled `TRUE`/`FALSE` (any case), and dates must be ISO 8601
#' `YYYY-MM-DD`. Columns named in `except` are left as character.
#'
#' @noRd
#' @keywords internal
convert_col_types <- function (df, except = character(0)) {

  num_re  <- "^[-+]?((0|[1-9][0-9]*)(\\.[0-9]*)?|\\.[0-9]+)([eE][-+]?[0-9]+)?$"
  bool_re <- "^(TRUE|FALSE|True|False|true|false)$"
  date_re <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}$"

  for (i in setdiff(names(df), except)) {

    x <- df[[i]]
    if (!is.character(x)) next

    v <- x[!is.na(x)]
    if (length(v) == 0) next

    if (all(grepl(bool_re, v))) {
      df[[i]] <- as.logical(toupper(x))

    } else if (all(grepl(num_re, v))) {
      df[[i]] <- as.numeric(x)

    } else if (all(grepl(date_re, v))) {
      d <- as.Date(x, format = "%Y-%m-%d")
      if (!anyNA(d[!is.na(x)])) df[[i]] <- d
    }
  }

  return (df)
}



#' Split lineage strings (`"Bacteria; Firmicutes; ..."`) into a character
#' matrix with one column per rank. Rows are padded with `NA`; empty ranks
#' become `NA`. Both `;` and `|` are accepted as separators.
#'
#' @noRd
#' @keywords internal
split_lineages <- function (x, sep = "[|;]") {

  x  <- as.character(x)
  ok <- !is.na(x)

  parts <- strsplit(paste0(ifelse(ok, x, ""), "\x01"), sep)
  lens  <- lengths(parts)
  vals  <- unlist(parts, use.names = FALSE)
  last  <- cumsum(lens)
  vals[last] <- sub("\x01$", "", vals[last])
  vals  <- trimws(vals)

  mat <- matrix(NA_character_, nrow = length(x), ncol = max(lens, 1L))
  mat[cbind(rep(seq_along(x), lens), sequence(lens))] <- vals
  mat[!nzchar(mat)] <- NA
  mat[!ok, ] <- NA

  return (mat)
}



#' Determine if a file path, URL, or literal data string.
#' is 'text', 'json', 'hdf5', or 'excel' format.
#'
#' @noRd
#' @keywords internal
#'
#' @param src  A file path, URL, or literal data string.
#' @param filename  The string to examine for file extensions.
#'
#' @return 'text', 'json', 'hdf5', or 'excel'.
#'         'AsIs' class when `src` is literal data rather than a file.
#'

file_type <- function (src, filename = src) {

  if (!is_nz_string(src)) cli_abort('`src` must be a string, not {.type {src}}.')
  if (is_url(src))        cli_abort('`src` cannot be a URL: {.url {src}}.')

  # Literal data: JSON (allowing a BOM and leading whitespace) or multi-line text.
  lead <- substr(trimws(strip_bom(src), which = "left"), 1, 1)
  if (lead %in% c('{', '['))          return (I('json'))
  if (grepl('\n', src, fixed = TRUE)) return (I('text'))

  if (!file.exists(src)) cli_abort('File not found: {filename}.')

  x <- sub('\\.(gz|bz2)$', '', tolower(filename))
  if (grepl('\\.(xls|xlsx)$', x))        return ('excel')
  if (grepl('\\.(txt|tsv|csv|tab)$', x)) return ('text')
  if (grepl('\\.json$', x))              return ('json')
  if (grepl('\\.(hdf|hdf5|h5)$', x))     return ('hdf5')

  con  <- base::file(src)
  cls  <- summary(con)$class
  peek <- suppressWarnings(readChar(con, 64L, useBytes = TRUE))
  close(con)

  if (cls %in% c('gzfile', 'bzfile', 'xzfile')) {
    lead <- substr(trimws(strip_bom(peek), which = "left"), 1, 1)
    return (ifelse(lead %in% c('{', '['), 'json', 'text'))
  }

  con  <- base::file(src, 'rb', raw = TRUE)
  peek <- readBin(con, 'raw', 64L)
  close(con)

  if (identical(peek[1:4], charToRaw('\x89HDF')))                          return ('hdf5')
  if (identical(peek[1:4], charToRaw('\x50\x4B\x03\x04')))                 return ('excel')
  if (identical(peek[1:8], charToRaw('\xD0\xCF\x11\xE0\xA1\xB1\x1A\xE1'))) return ('excel')

  # Skip a UTF-8 BOM and whitespace before checking for JSON.
  if (length(peek) >= 3 && identical(peek[1:3], as.raw(c(0xEF, 0xBB, 0xBF)))) peek <- peek[-(1:3)]
  peek <- peek[!peek %in% charToRaw(" \t\r\n")]
  if (length(peek) > 0 && peek[[1]] %in% charToRaw('{['))                  return ('json')

  return ('text')
}



#' Pull URLs to a local file. Save literal data to a file.
#' `$cleanup()` will remove any files created by this function.
#'
#' @noRd
#' @keywords internal
#'
#' @param src  A file path, URL, or literal data string.
#'
#' @return list(id, path, format, cleanup)
#'

as_filepath <- function (src) {

  if (is_url(src)) {

    id      <- if (nchar(src) > 100) basename(src) else src
    path    <- tempfile()
    cleanup <- function () unlink(path)
    x <- try(silent = TRUE, suppressWarnings(download.file(url = src, destfile = path, quiet = TRUE, mode = "wb")))
    if (!eq(0L, x)) {
      cleanup()
      cli_abort("Cannot retrieve URL {.url {src}}: {x}")
    }

    format <- file_type(src = path, filename = src)


  } else {

    format <- file_type(src = src)

    if (inherits(format, 'AsIs')) {

      # Save literal data string to a temp file.
      id      <- paste('Parsed from', format, 'string')
      path    <- tempfile()
      cleanup <- function () unlink(path)
      writeChar(object = enc2utf8(src), con = path, eos = NULL, useBytes = TRUE)
      class(format) <- NULL

    } else {
      path    <- normalizePath(src, winslash = '/')
      id      <- if (nchar(path) > 100) basename(path) else path
      cleanup <- function () NULL
    }
  }

  result <- list(id = id, path = path, format = format, cleanup = cleanup)


  return (result)
}



#' Find the column holding sample/OTU identifiers.
#'
#' Preference order: a column already named `dot`; row names; a column whose
#' (case-insensitive, `#`-stripped) name is one of `aliases`; and finally a
#' column with unique values that overlap `ids`. More than one overlapping
#' candidate is an error.
#'
#' @noRd
#' @keywords internal
#'
#' @return The table with the identifier column named `dot` and moved to the
#'         first position (row names are converted to that column).
#'
find_id_column <- function (tbl, ids, dot, aliases, what) {

  has_rn <- tibble::has_rownames(tbl)

  if (!hasName(tbl, dot) && !has_rn) {

    nm <- tolower(trimws(sub('^#\\s*', '', colnames(tbl))))
    i  <- which(nm %in% aliases)

    if (length(i) > 0) {
      i <- i[order(match(nm[i], aliases))][[1]]

    } else {
      cand <- which(vapply(tbl, function (x) {
        x <- trimws(as.character(x))
        !as.logical(anyDuplicated(x)) && any(x %in% ids)
      }, TRUE))

      if (length(cand) > 1)
        cli_abort(c(
          'x' = "Can't tell which column holds the {what} names.",
          'i' = "Candidate columns: {.val {colnames(tbl)[cand]}}.",
          'i' = "Rename the correct column to {.val {dot}}." ))

      if (length(cand) == 0)
        cli_abort(c(
          'x' = "{what} table must have row names or a {.val {dot}} column.",
          'i' = "No column matched the {what} names {.val {head(ids, 3)}}, ..." ))

      i <- cand
    }

    colnames(tbl)[[i]] <- dot
  }


  if (has_rn && hasName(tbl, dot)) {
    if (!all(rownames(tbl) == paste0('', seq_len(nrow(tbl)))))
      if (!all(rownames(tbl) == tbl[[dot]]))
        cli_abort(c(x = "Row names are not allowed when a {.val {dot}} column is present."))
    tbl %<>% tibble::remove_rownames()
  }

  if (tibble::has_rownames(tbl)) { tbl %<>% tibble::rownames_to_column(dot)
  } else                         { tbl %<>% relocate(dplyr::all_of(dot)) }

  return (tbl)
}



#' @noRd
#' @keywords internal
import_metadata <- function (value, sids) {

  #________________________________________________________
  # Erasing all the metadata.
  #________________________________________________________
  if (is.null(value)) return (tibble(.sample = sids))

  #________________________________________________________
  # Load metadata from a file.
  #________________________________________________________
  untyped <- FALSE
  if (is_scalar_character(value)) {
    value   <- import_table(src = value)
    untyped <- isTRUE(attr(value, 'rbiom_untyped'))
    attr(value, 'rbiom_untyped') <- NULL
  }

  #________________________________________________________
  # Convert data.frame/matrix to tibble.
  #________________________________________________________
  value <- tibble::as_tibble(value, rownames = NA, .name_repair = trimws)


  #________________________________________________________
  # Find the column with sample names; make it '.sample' and first.
  #________________________________________________________
  value <- find_id_column(
    tbl     = value,
    ids     = sids,
    dot     = '.sample',
    what    = "sample",
    aliases = c('sample', 'sample_id', 'sampleid', 'sample id', 'sample-id',
                'sample_name', 'samplename', 'sample name', 'id', '') )


  #________________________________________________________
  # Disallow column names starting with a "."
  #________________________________________________________
  if (length(i <- grep("^\\.", colnames(value)[-1], value = TRUE)) != 0)
    cli_abort(c(x = "Metadata field{?s} can't start with a '.': {.val {i}}."))


  #________________________________________________________
  # Text files: type the columns now that the IDs are known.
  #________________________________________________________
  if (untyped) value <- convert_col_types(value, except = '.sample')


  #________________________________________________________
  # Convert all character cols (except .sample) to factor.
  #________________________________________________________
  value[['.sample']] %<>% as.character() %>% trimws()
  for (i in colnames(value)[-1])
    if (is.character(value[[i]])) value[[i]] %<>% trimws() %>% as.factor()


  #________________________________________________________
  # Intersect with sample names from $counts.
  #________________________________________________________
  expected <- sids
  provided <- unique(value[['.sample']])

  if (length(intersect(expected, provided)) == 0)
    cli_abort(c(
      'i' = "No matching sample names.",
      '*' = "Expected: {.val {expected}}",
      '*' = "Provided: {.val {provided}}",
      'x' = "Can't subset to zero samples." ))

  if (length(i <- setdiff(expected, provided)) > 0)
    cli_warn(c('i' = paste(
      "Dropping {length(i)} sample{?s} from biom object",
      "since they are not in the new metadata: {.val {i}}." )))

  if (length(i <- setdiff(provided, expected)) > 0)
    cli_warn(c('i' = paste(
      "Ignoring metadata for {length(i)} sample{?s}",
      "not currently in biom object: {.val {i}}." )))

  value <- value[value[['.sample']] %in% expected,,drop=FALSE]
  remove('expected', 'provided', 'i')


  #________________________________________________________
  # Drop duplicate sample mappings. Warn if non-identical.
  #________________________________________________________
  if (as.logical(anyDuplicated(value[['.sample']]))) {

    not_equal <- function (x) (length(unique(as.character(x))) > 1)
    to_string <- function (i, x) toJSON(as.vector(x[i,,drop=FALSE]), auto_unbox = TRUE)

    for (sid in unique(value[['.sample']][duplicated(value[['.sample']])])) {

      x <- value[value[['.sample']] == sid,,drop=FALSE]
      x <- x[,vapply(x, not_equal, TRUE),drop=FALSE]
      if (ncol(x) == 0) next

      rlang::warn(
        message = c(i = paste0("Discarding subsequent metadata for sample ", sid)),
        body    = setNames(
          object = sapply(1:nrow(x), to_string, x),
          nm     = c('v', rep('x', nrow(x) - 1)) ))
    }

    value <- value[!duplicated(value[['.sample']]),,drop=FALSE]
    remove("not_equal", "to_string")
  }


  return (value)
}


#' @noRd
#' @keywords internal
import_taxonomy <- function (value, otus) {

  #________________________________________________________
  # Erasing all the taxonomy.
  #________________________________________________________
  if (is.null(value)) return (tibble(.otu = otus))


  #________________________________________________________
  # Load taxonomy from a file.
  #________________________________________________________
  if (is_scalar_character(value)) {
    value <- import_table(src = value)
    attr(value, 'rbiom_untyped') <- NULL
  }

  #________________________________________________________
  # Convert data.frame/matrix to tibble.
  #________________________________________________________
  value <- tibble::as_tibble(value, rownames = NA, .name_repair = trimws)


  #________________________________________________________
  # Find the column with OTU names; make it '.otu' and first.
  #________________________________________________________
  value <- find_id_column(
    tbl     = value,
    ids     = otus,
    dot     = '.otu',
    what    = "OTU",
    aliases = c('otu', 'otu_id', 'otuid', 'otu id', 'otu-id',
                'feature id', 'feature_id', 'featureid', 'feature-id',
                'asv', 'asv_id', 'asvid', 'representative_sequence', 'id', '') )


  #________________________________________________________
  # Disallow column names starting with a "."
  #________________________________________________________
  if (length(i <- grep("^\\.", colnames(value)[-1], value = TRUE)) != 0)
    cli_abort(c(x = "Taxonomy rank{?s} can't start with a '.': {.val {i}}."))


  #________________________________________________________
  # Convert all cols (except .otu) to factor.
  #________________________________________________________
  value[['.otu']] %<>% as.character() %>% trimws()
  for (i in seq_len(ncol(value))[-1])
    if (!is.factor(value[[i]])) value[[i]] %<>% as.character() %>% trimws() %>% as.factor()


  #________________________________________________________
  # Intersect with OTU names from $counts.
  #________________________________________________________
  expected <- otus
  provided <- unique(value[['.otu']])

  if (length(intersect(expected, provided)) == 0)
    cli_abort(c(
      'i' = "No matching OTU names between counts and taxonomy.",
      '*' = "Expected: {.val {expected}}",
      '*' = "Provided: {.val {provided}}",
      'x' = "Can't subset to zero OTUs." ))

  if (length(i <- setdiff(expected, provided)) > 0)
    cli_warn(c('i' = paste(
      "Dropping {length(i)} OTU{?s} from biom object",
      "since they are not in the new taxonomy: {.val {i}}." )))

  if (length(i <- setdiff(provided, expected)) > 0)
    cli_warn(c('i' = paste(
      "Ignoring taxonomy for {length(i)} OTU{?s}",
      "not currently in biom object: {.val {i}}." )))

  value <- value[value[['.otu']] %in% expected,,drop=FALSE]
  remove('expected', 'provided', 'i')


  #________________________________________________________
  # Drop duplicate OTU mappings. Warn if non-identical.
  #________________________________________________________
  if (as.logical(anyDuplicated(value[['.otu']]))) {

    if (ncol(value) == 1) {
      keep <- which(!duplicated(value[['.otu']]))

    } else {

      keep <- sapply(split(1:nrow(value), value[['.otu']]), function (i) {

        if (length(i) == 1) return (i)

        lineages <- apply(value[i,-1,drop=FALSE], 1L, paste, collapse = "; ")
        longest  <- which.max(nchar(lineages))

        if (length(unique(lineages)) > 1)
          rlang::warn(
            message = c(i = paste0("Discarding less verbose mapping for OTU ", value[i[[1]],1])),
            body    = setNames(
              object = lineages,
              nm     = ifelse(seq_along(lineages) == longest, 'v', 'x') ))

        return (i[[longest]])
      })
    }

    value <- value[sort(unname(keep)),,drop=FALSE]
    remove("keep")
  }


  #________________________________________________________
  # Drop 'Confidence' column (from qiime2).
  #________________________________________________________
  if (identical(tail(colnames(value), 1), 'Confidence'))
    value <- value[,seq_len(ncol(value) - 1)]


  #________________________________________________________
  # Split semicolon-separated lineages into one column per rank.
  #________________________________________________________
  if (ncol(value) == 2 && any(grepl("[|;]", as.character(value[[2]])))) {
    ssv <- split_lineages(value[[2]])
    colnames(ssv) <- default_taxa_ranks(ncol(ssv))
    ssv <- as_tibble(ssv, .name_repair = "minimal")
    for (i in seq_len(ncol(ssv))) ssv[[i]] <- as.factor(ssv[[i]])
    value <- dplyr::bind_cols(value[,'.otu'], ssv)
    remove('ssv')
  }


  #________________________________________________________
  # Drop unused factor levels.
  #________________________________________________________
  for (i in seq_len(ncol(value))[-1])
    value[[i]] %<>% { factor(., levels = intersect(levels(.), unique(.))) }


  return (value)
}

