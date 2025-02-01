#' Parse counts, metadata, taxonomy, and phylogeny from a BIOM file.
#' 
#' @inherit documentation_return.biom return
#'
#' @param src   Input data as either a file path, URL, or JSON string.
#'        BIOM files can be formatted according to 
#'        version 1.0 (JSON) or 2.1 (HDF5)
#'        \href{http://biom-format.org/documentation/}{specifications}, or as 
#'        classical tabular format. URLs must begin with \code{http://},
#'        \code{https://}, \code{ftp://}, or \code{ftps://}. JSON files must 
#'        have \code{\{} as their first character. Compressed (gzip or bzip2) 
#'        BIOM files are also supported. NOTE: to read HDF5 formatted BIOM 
#'        files, the BioConductor R package \code{rhdf5} must be installed.
#' 
#' @param ...   Properties to set in the new rbiom object, for example, 
#'        `metadata`, `id`, `comment`, or `tree`.
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
#'     # Phylogenetic tree
#'     biom$tree %>%
#'       tree_subset(1:10) %>%
#'       plot()
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
    
    require_package('rhdf5', 'to read HDF5 formatted BIOM files')
    if (!H5Fis_hdf5(fp))
      cli_abort("HDF5 file not recognized by rhdf5: {.file {src}}")
    
    h5        <- read_biom_hdf5(fp)
    counts    <- parse_hdf5_counts(h5)
    info      <- parse_hdf5_info(h5)
    sequences <- if (!hasName(dots, 'sequences')) parse_hdf5_sequences(h5)
    taxonomy  <- if (!hasName(dots, 'taxonomy'))  parse_hdf5_taxonomy(h5)
    metadata  <- if (!hasName(dots, 'metadata'))  parse_hdf5_metadata(h5)
    phylogeny <- if (!hasName(dots, 'tree'))      parse_hdf5_tree(h5)
    
    H5Fclose(h5)
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
    
    mtx       <- read_biom_tsv(fp)
    counts    <- parse_tsv_counts(mtx)
    info      <- list(id = fpath$id, type = 'OTU table')
    taxonomy  <- if (!hasName(dots, 'taxonomy')) parse_tsv_taxonomy(mtx)
    metadata  <- if (!hasName(dots, 'metadata')) data.frame(row.names=colnames(counts))
    sequences <- NULL
    phylogeny <- NULL
    
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
      'metadata'  = metadata, 
      'taxonomy'  = taxonomy, 
      'sequences' = sequences, 
      'tree'      = phylogeny, 
      'id'        = info$id, 
      'date'      = info$date, 
      'comment'   = info$comment ))
  
  args <- args[!duplicated(names(args))]
  biom <- do.call(rbiom$new, args)
  
  
  #________________________________________________________
  # Attach as_rbiom() call to provenance tracking
  #________________________________________________________
  # cl <- match.call()
  # cl[[1]] <- as.name("as_rbiom")
  # for (i in seq_along(cl)[-1]) {
  #   val <- eval.parent(cl[[i]])
  #   if (all(nchar(val) <= 200)) # Don't dump huge JSON strings
  #     cl[i] <- list(val)
  # }
  
  
  return (biom)
}



read_biom_tsv <- function (fp) {

  #________________________________________________________
  # Read in all the lines from the file
  #________________________________________________________

  lines <- tryCatch(
    error = function (e) stop("can't parse csv/tsv file - ", e),
    expr  = readLines(fp, warn=FALSE))
  
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) stop("Input file is empty.")


  #________________________________________________________
  # Find the header line; check if ',' or '\t' comes first.
  #________________________________________________________
  
  header <- lines[[c(grep("^#(Sample|OTU)", lines, ignore.case = TRUE), 1)[[1]]]]
  sep    <- head(c(intersect(strsplit(header, '',)[[1]], c(',', '\t')), ','), 1)


  #________________________________________________________
  # Read the file again, includes quote support.
  #________________________________________________________
  
  mat <- tryCatch(
    error = function (e) stop("can't parse file - ", e),
    expr  = as.matrix(read.table(
      file         = fp, 
      sep          = sep, 
      col.names    = strsplit(header, sep, fixed = TRUE)[[1]],
      row.names    = 1,
      colClasses   = 'character', 
      allowEscapes = TRUE,
      skipNul      = TRUE,
      as.is        = TRUE,
      check.names  = FALSE )))


  #________________________________________________________
  # Check for duplicate row/col names
  #________________________________________________________
  
  if (length(x <- rownames(mat)[duplicated(rownames(mat))])) cli_abort('Duplicated IDs: {.val {x}}')
  if (length(x <- colnames(mat)[duplicated(colnames(mat))])) cli_abort('Duplicated IDs: {.val {x}}')
  

  return (mat)
}

read_biom_json <- function (fp) {

  json <- tryCatch(
    expr  = read_json(path = fp),
    error = function (e) cli_abort("Can't read JSON file {.file {fp}}: {e}") )
  
  if (length(x <- setdiff(c('data', 'matrix_type', 'rows', 'columns'), names(json))))
    cli_abort("BIOM file is missing {.var {x}}")

  return (json)
}

read_biom_hdf5 <- function (fp) {
  
  h5 <- try(H5Fopen(fp, 'H5F_ACC_RDONLY', native = TRUE), silent=TRUE)
  if (inherits(h5, "try-error")) {
    try(H5Fclose(h5), silent=TRUE)
    stop(sprintf("Unable to parse HDF5 file. %s", as.character(h5)))
  }
  
  entries  <- with(h5ls(h5), paste(sep="/", group, name))
  expected <- c( "/observation/matrix/indptr", "/observation/matrix/indices", 
                 "/observation/ids", "/observation/matrix/data", "/sample/ids" )
  
  missing <- setdiff(expected, entries)
  if (length(missing) > 0) {
    try(H5Fclose(h5), silent=TRUE)
    stop(sprintf("BIOM file requires: %s", paste(collapse=",", missing)))
  }
  
  return (h5)
}



parse_json_metadata <- function (json) {
  
  if (as.logical(anyDuplicated(names(json$columns[[1]]$metadata))))
    cli_abort("Metadata column names are not unique.")
  
  df <- as.data.frame(
    check.names      = FALSE,
    stringsAsFactors = FALSE,
    vapply(names(json$columns[[1]]$metadata), function (i) {
      sapply(json$columns, function (x) {
        val <- as.character(x[['metadata']][[i]])
        ifelse(is_null(val), NA_character_, val)
      })
    }, character(length(json$columns)), USE.NAMES=TRUE)
  )
  
  rowIDs <- sapply(json$columns, function (x) unlist(x$id))
  if (length(rowIDs) == 1) df <- data.frame(t(df))
  rownames(df) <- rowIDs
  

  # Convert strings to numbers
  for (i in seq_len(ncol(df))) {
    df[[i]][which(df[[i]] == "NA")] <- NA
    if (all(grepl("^-?(0|[1-9]\\d*)(\\.\\d*[1-9]|)(e[\\+\\-]{0,1}[0-9]+|)$", na.omit(df[[i]]))))
      df[[i]] <- as.numeric(df[[i]])
  }

  return (df)
}

parse_hdf5_metadata <- function (h5) {

  keys <- names(h5$sample$metadata)
  vals <- lapply(keys, function (k) {
    
    obj_class <- typeof(h5$sample$metadata[[k]])
    
    if (eq(obj_class, "character")) {
      as(h5$sample$metadata[[k]], "character")
    } else {
      as(h5$sample$metadata[[k]], "numeric")
    }
  })

  as.data.frame(
    row.names        = h5$sample$ids,
    check.names      = FALSE,
    stringsAsFactors = FALSE,
    setNames(vals, keys)
  )
}



parse_json_info <- function (json) {
  info <- list()
  for (i in setdiff(names(json), c("rows", "columns", "data", "phylogeny")))
    info[[i]] <- if (is_null(unlist(json[[i]]))) NA else unlist(json[[i]])
  return (info)
}

parse_hdf5_info <- function (h5) {
  attrs <- h5readAttributes(h5, "/")
  for (i in names(attrs)) {
    attrs[[i]] <- as(attrs[[i]], typeof(attrs[[i]]))
  }
  return (attrs)
}



parse_tsv_counts <- function (mtx) {

  # Only keep columns that are all numbers
  allNumbers <- function (x) all(grepl("^\\d+(\\.\\d+|)([Ee][\\+\\-]{0,1}[0-9]+|)$", x))
  mtx <- mtx[, apply(mtx, 2L, allNumbers), drop=FALSE]
  mtx <- matrix(
    data     = as.numeric(mtx),
    nrow     = nrow(mtx),
    dimnames = list(rownames(mtx), colnames(mtx)) )

  if (length(mtx) == 0 || sum(mtx) == 0)
    stop("No abundance counts.")

  as.simple_triplet_matrix(mtx)
}

parse_json_counts <- function (json) {

  if (!any(c('sparse', 'dense') %in% json$matrix_type))
    stop("BIOM file's matrix_type must be either 'sparse' or 'dense'")

  if (length(json$data) == 0)
    stop("BIOM file does not have any count data.")

  TaxaIDs   <- sapply(json$rows,    function (x) unlist(x$id))
  SampleIDs <- sapply(json$columns, function (x) unlist(x$id))

  if (json$matrix_type == "sparse")
    counts <- simple_triplet_matrix(
                    i        = sapply(json$data, function (x) x[[1]]) + 1,
                    j        = sapply(json$data, function (x) x[[2]]) + 1,
                    v        = sapply(json$data, function (x) x[[3]]),
                    nrow     = length(TaxaIDs),
                    ncol     = length(SampleIDs),
                    dimnames = list(TaxaIDs, SampleIDs))

  if (json$matrix_type == "dense")
    counts <- as.simple_triplet_matrix(
                    matrix(
                      data     = unlist(json$data),
                      byrow    = TRUE,
                      nrow     = length(TaxaIDs),
                      ncol     = length(SampleIDs),
                      dimnames = list(TaxaIDs, SampleIDs)))

  return (counts)
}

parse_hdf5_counts <- function (h5) {

  indptr <- as.numeric(h5$observation$matrix$indptr)

  simple_triplet_matrix(
        i        = unlist(sapply(1:(length(indptr)-1), function (i) rep(i, diff(indptr[c(i,i+1)])))),
        j        = as.numeric(h5$observation$matrix$indices) + 1,
        v        = as.numeric(h5$observation$matrix$data),
        nrow     = length(h5$observation$ids),
        ncol     = length(h5$sample$ids),
        dimnames = list(
          as.character(h5$observation$ids),
          as.character(h5$sample$ids)
      ))
}



parse_tsv_taxonomy <- function (mtx) {

  # Discard columns that are all numbers
  allNumbers <- function (x) all(grepl("^\\d+(\\.\\d+|)(e[\\+\\-]{0,1}[0-9]+|)$", x))
  mtx <- mtx[, !apply(mtx, 2L, allNumbers), drop=FALSE]

  # Look for a taxonomy column, otherwise try to parse the taxa IDs
  txCol      <- head(grep("taxonomy", colnames(mtx), ignore.case=TRUE), 1)
  taxaNames  <- if (length(txCol)) mtx[,txCol] else rownames(mtx)
  taxa_table <- strsplit(taxaNames, "[\\|;]\\ *")

  # Handle instances where some taxa strings have more levels than others.
  n <- sapply(taxa_table, length)
  m <- max(n)
  i <- which(n < m)
  taxa_table[i] <- lapply(taxa_table[i], function (x) c(x, rep(NA, m - length(x))))
  taxa_table    <- matrix(unlist(taxa_table), nrow=length(taxa_table), byrow=TRUE)

  rownames(taxa_table) <- rownames(mtx)
  colnames(taxa_table) <- default_taxa_ranks(ncol(taxa_table))

  # Better to return a completely empty table than one with just the taxa IDs
  if (eq(unname(taxa_table[,1]), rownames(taxa_table)))
    taxa_table <- taxa_table[,-1,drop=FALSE]

  return (taxa_table)
}

parse_json_taxonomy <- function (json) {

  taxa_table <- sapply(json$rows, simplify = "array", function (x) {
    
    taxaNames <- unlist(x$metadata$taxonomy)
    
    # The taxa names aren't where they're supposed to be
    if (is_null(taxaNames)) {
      
      # MicrobiomeDB puts the taxa string in the 'ID' field, e.g, {"metadata":null,"id":"Archaea;Euryarchaeota;...
      if (eq(json[['generated_by']], "MicrobiomeDB")) {
        taxaNames <- x[['id']]
      
      # Decontam omits the 'taxonomy' name, as in {"id":"Unc01pdq","metadata":[["Bacteria","__Fusobacteriota", ...
      } else if (is_null(names(x$metadata)) && length(x$metadata) == 1) {
        taxaNames <- unlist(x$metadata[[1]])
        
      } else {
        return (unlist(x$id))
      }
    }
    
    
    if (length(taxaNames) == 1) {
      if (nchar(taxaNames) == 0) return (unlist(x$id))
      return (strsplit(taxaNames, "[\\|;]\\ *")[[1]])
    }
    
    return (taxaNames)
  })

  if (inherits(taxa_table, "matrix")) {
    # 'Normal' situation of same number of ranks (>1) per taxa string.
    
    taxa_table <- t(taxa_table)
    
  } else if (inherits(taxa_table, "list")) {
    # Handle instances where some taxa strings have more levels than others.
    
    n <- sapply(taxa_table, length)
    m <- max(n)
    i <- which(n < m)
    taxa_table[i] <- lapply(taxa_table[i], function (x) c(x, rep(NA, m - length(x))))
    taxa_table    <- t(simplify2array(taxa_table))
    
  } else {
    # There are no taxa strings, just the ID.
    
    taxa_table <- as.matrix(taxa_table, ncol=1)
  }

  rownames(taxa_table) <- sapply(json$rows, function (x) unlist(x$id))
  colnames(taxa_table) <- default_taxa_ranks(ncol(taxa_table))

  #taxa_table <- PB.SanitizeTaxonomy(taxa_table, env=parent.frame())

  return (taxa_table)
}

parse_hdf5_taxonomy <- function (h5) {
  
  if ("taxonomy" %in% names(h5$observation$metadata)) {
    taxa_table           <- t(h5$observation$metadata$taxonomy)
    rownames(taxa_table) <- as.character(h5$observation$ids)
    colnames(taxa_table) <- default_taxa_ranks(ncol(taxa_table))
    
  } else {
    ids        <- as.character(h5$observation$ids)
    taxa_table <- matrix(
      data     = character(0), 
      nrow     = length(ids), 
      ncol     = 0, 
      dimnames = list(ids, character(0)) )
  }

  #taxa_table <- PB.SanitizeTaxonomy(taxa_table, env=parent.frame())
    
  return (taxa_table)
}




parse_json_sequences <- function (json) {
  
  # No sequence information
  if (!any(sapply(json$rows, function (x) 'sequence' %in% names(x$metadata) )))
    return (NULL)
  
  ids  <- sapply(json$rows, simplify = "array", function (x) { unlist(x$id) })
  seqs <- sapply(json$rows, simplify = "array", function (x) {
    if (is_null(x$metadata$sequence)) NA else unlist(x$metadata$sequence)
  })
  
  res <- setNames(seqs, ids)
  
  return (res)
}

parse_hdf5_sequences <- function (h5) {
  
  # No sequence information
  if (!"sequences" %in% names(h5$observation$metadata))
    return (NULL)
  
  ids  <- as.character(h5$observation$ids)
  seqs <- as.character(h5$observation$metadata$sequences)
  
  res <- setNames(seqs, ids)
    
  return (res)
}



parse_json_tree <- function (json) {
  tree <- as.character(unlist(json[['phylogeny']]))
  if (isTRUE(startsWith(tree, '('))) tree else NULL
}

parse_hdf5_tree <- function (h5) {
  tree <- as.character(h5$observation$`group-metadata`$phylogeny)
  if (isTRUE(startsWith(tree, '('))) tree else NULL
}




default_taxa_ranks <- function (n) {
  if (n >= 6 && n <= 8)
    return(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:n])
  sprintf("Level.%i", 1:n)
}




