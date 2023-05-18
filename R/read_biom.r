#' Extracts counts, metadata, taxonomy, and phylogeny from a biom file.
#'
#' @param src  Input data as either a file path, URL, or JSON string.
#'        \code{read_biom} can read BIOM files formatted according to both the
#'        version 1.0 (JSON) and 2.1 (HDF5)
#'        \href{http://biom-format.org/documentation/}{specifications} as well
#'        as classical tabular format. URLs must begin with \kbd{http://},
#'        \kbd{https://}, \kbd{ftp://}, or \kbd{ftps://}. JSON files must have
#'        \code{\{} as their first non-whitespace character. Compressed (gzip
#'        or bzip2) BIOM files are also supported. NOTE: to read HDF5 formatted
#'        BIOM files, the BioConductor R package \code{rhdf5} must be 
#'        installed.
#'     
#' @param tree  The default value of \code{auto} will read the tree from the 
#'        BIOM file specified in \code{src}, if present. The value \code{TRUE} 
#'        will do the same, but will generate an error message if a tree is not 
#'        present. Setting \code{tree=FALSE} will return a \code{BIOM} object 
#'        without any tree data. You may also provide a file path, URL, or
#'        Newick string to load that tree data into the final \code{BIOM} 
#'        object.
#'     
#' @param prune  Should samples and taxa with zero observations be discarded?
#'        (Default: \code{cleanup})
#' 
#' @param cleanup  Renames ambiguous taxons and removes leading underscores. 
#'        Also converts character metadata into factors and dates based on 
#'        heuristics. (Default: FALSE)
#' 
#' @return A \code{BIOM} class object containing the parsed data. This object
#'     can be treated as a list with the following named elements:
#'     \describe{
#'         \item{counts}{A numeric \code{slam} sparse matrix of observation
#'         counts. Taxa (OTUs) as rows and samples as columns.}
#'         \item{metadata}{A data frame containing any embedded metadata.
#'         Row names are sample IDs.}
#'         \item{taxonomy}{Character matrix of taxonomic names, if given.
#'         Row names are taxa (OTU) IDs. Column rows are named Kingdom,
#'         Phylum, Class, Order, Family, Genus, Species, and Strain, or
#'         TaxLvl.1, TaxLvl.2, ... , TaxLvl.N when more than 8 levels of
#'         taxonomy are encoded in the biom file.}
#'         \item{phylogeny}{An object of class \code{phylo} defining the
#'         phylogenetic relationships between the taxa. Although the
#'         official specification for BIOM only includes phylogenetic trees
#'         in BIOM version 2.1, if a BIOM version 1.0 file includes a
#'         \code{phylogeny} entry with newick data, then it will be loaded
#'         here as well. The \pkg{ape} package has additional functions
#'         for working with \code{phylo} objects.}
#'         \item{sequences}{A named character vector, where the names are
#'         taxonomic identifiers and the values are the sequences they
#'         represent. These values are not part of the official BIOM
#'         specification, but will be read and written when defined.}
#'         \item{info}{A list of other attributes defined in the BIOM file,
#'         such as \code{id}, \code{type}, \code{format}, \code{format_url},
#'         \code{generated_by}, \code{date}, \code{matrix_type},
#'         \code{matrix_element_type}, \code{Comment}, and \code{shape}}
#'        }
#'
#'     \code{metadata}, \code{taxonomy}, and \code{phylogeny} are optional
#'     components of the BIOM file specification and therefore will be empty
#'     in the returned object when they are not provided by the BIOM file.
#' @export
#' @examples
#'     library(rbiom)
#'
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read_biom(infile)
#'
#'     summary(biom)
#'
#'     # Taxa Abundances
#'     as.matrix(biom$counts[1:4,1:4])
#'
#'     top5 <- names(head(rev(sort(slam::row_sums(biom$counts))), 5))
#'     biom$taxonomy[top5,c('Family', 'Genus')]
#'     as.matrix(biom$counts[top5, 1:6])
#'
#'     # Metadata
#'     table(biom$metadata$Sex, biom$metadata$`Body Site`)
#'     sprintf("Mean age: %.1f", mean(biom$metadata$Age))
#'
#'     # Phylogenetic tree
#'     tree <- biom$phylogeny
#'     top5.tree <- rbiom::subtree(tree, top5)
#'     ape::plot.phylo(top5.tree)
#'


read_biom <- function (src, tree='auto', prune=cleanup, cleanup=FALSE) {

  #________________________________________________________
  # Sanity check input values
  #________________________________________________________

  if (length(src) != 1 || !is(src, "character"))
    stop(simpleError("Data source for read_biom() must be a single string."))


  #________________________________________________________
  # Get the url or text for the src/BIOM data into a file
  #________________________________________________________

  if (length(grep("^(ht|f)tps{0,1}://.+", src)) == 1) {

    fp <- tempfile()
    on.exit(unlink(fp), add=TRUE)
    
    # To do: switch to curl::curl_download
    if (!identical(0L, try(download.file(src, fp, quiet=TRUE), silent=TRUE)))
        stop(simpleError(sprintf("Cannot retrieve URL %s", src)))

  } else if (length(grep("^[ \t\n]*\\{", src)) == 1) {

    fp <- tempfile(fileext=".biom")
    on.exit(unlink(fp), add=TRUE)
    if (!is_null(try(writeChar(src, fp, eos=NULL), silent=TRUE)))
        stop(simpleError(sprintf("Cannot write text to file %s", fp)))

  } else {

    fp <- normalizePath(src)
  }

  if (!file.exists(fp))
    stop(simpleError(sprintf("Cannot locate file %s", fp)))
  
  # If date field is missing from the biom file, use file creation time.
  fp_date <- strftime(file.info(fp)[['ctime']], "%Y-%m-%dT%H:%M:%SZ", tz="UTC")
  

  #________________________________________________________
  # Decompress files that are in gzip or bzip2 format
  #________________________________________________________
  
  format <- guess_format(fp)
  
  if (endsWith(format, ".gz")) {
    fp <- R.utils::gunzip(fp, destname=tempfile(), remove=FALSE)
    on.exit(unlink(fp), add=TRUE)
    
  } else if (endsWith(format, ".bz2")) {
    fp <- R.utils::bunzip2(fp, destname=tempfile(), remove=FALSE)
    on.exit(unlink(fp), add=TRUE)
  }
  
  
  
  #________________________________________________________
  # Process the file according to its internal format
  #________________________________________________________
  
  if (startsWith(format, "hdf5")) {
    
    #___________________#
    # HDF5 file format  #
    #___________________#
    
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
      stop(simpleError(paste0(
        "\n",
        "Error: rbiom requires the R package 'rhdf5' to be installed\n",
        "in order to read and write HDF5 formatted BIOM files.\n\n",
        "Please run the following commands to install 'rhdf5':\n",
        "   install.packages('BiocManager')\n",
        "   BiocManager::install('rhdf5')\n\n" )))
    }
    
    if (!rhdf5::H5Fis_hdf5(fp)) {
      stop(simpleError("HDF5 file not recognized by rhdf5."))
    }
    
    
    hdf5      <- read_biom_hdf5(fp)
    counts    <- parse_hdf5_counts(hdf5)
    sequences <- parse_hdf5_sequences(hdf5)
    taxonomy  <- parse_hdf5_taxonomy(hdf5)
    metadata  <- parse_hdf5_metadata(hdf5)
    info      <- parse_hdf5_info(hdf5)
    phylogeny <- parse_hdf5_tree(hdf5, tree)
    
    rhdf5::H5Fclose(hdf5)
    remove("hdf5")
    
  } else if (startsWith(format, "json")) {
    
    #___________________#
    # JSON file format  #
    #___________________#
    
    json      <- read_biom_json(fp)
    counts    <- parse_json_counts(json)
    sequences <- parse_json_sequences(json)
    taxonomy  <- parse_json_taxonomy(json)
    metadata  <- parse_json_metadata(json)
    info      <- parse_json_info(json)
    phylogeny <- parse_json_tree(json, tree)
    
    remove("json")
    
  } else {
    
    #___________________#
    # TSV file format   #
    #___________________#
    
    if (identical(tree, TRUE))
      stop(simpleError("It is impossible to load a phylogenetic tree from a BIOM file in tab-separated format."))
    
    mtx       <- read_biom_tsv(fp)
    counts    <- parse_tsv_counts(mtx)
    sequences <- NULL
    taxonomy  <- parse_tsv_taxonomy(mtx)
    metadata  <- data.frame(row.names=colnames(counts))
    info      <- list(id=tools::md5sum(fp)[[1]], type="OTU table")
    phylogeny <- NULL
  }
  
  
  #________________________________________________________
  # Return everything we've computed as a BIOM class object.
  #________________________________________________________
  
  if (!is.character(info[['comment']]))
    info[['comment']] <- NULL
  
  if (is_null(info[['date']])) {
    if (!is_null(info[['creation-date']])) {
      info[['date']] <- info[['creation-date']]
    } else {
      info[['date']] <- fp_date
    }
  }
  
  biom <- structure(
    class   = c("BIOM", "list"),
    comment = info[['comment']] %||% "",
    list( 'counts'    = counts,
          'metadata'  = metadata,
          'taxonomy'  = taxonomy,
          'phylogeny' = phylogeny,
          'sequences' = sequences,
          'info'      = info[names(info) != "comment"]
    )
  )
  
  
  #________________________________________________________
  # Clean up taxa names and metadata column classes.
  #________________________________________________________
  if (isTRUE(cleanup)) {
    biom[['taxonomy']] <- rbiom::taxonomy(biom, unc     = "grouped")
    biom[['metadata']] <- rbiom::metadata(biom, cleanup = TRUE)
  }
  
  #________________________________________________________
  # Discard samples/taxa with zero observations
  #________________________________________________________
  if (isTRUE(prune)) {
    biom <- repair(biom, prune=TRUE)
  }
  
  
  #________________________________________________________
  # Attach read_biom() call to provenance tracking
  #________________________________________________________
  cl <- match.call()
  cl[[1]] <- as.name("read_biom")
  for (i in seq_along(cl)[-1]) {
    val <- eval.parent(cl[[i]])
    if (all(nchar(val) <= 200)) # Don't dump huge JSON strings
      cl[i] <- list(val)
  }
  attr(biom, 'history') <- paste("biom <-", deparse1(cl))
  
  
  #________________________________________________________
  # Determine if these counts are pre-rarefied
  #________________________________________________________
  if (isTRUE(length(d <- depth(biom)) == 1))
    attr(biom, 'rarefaction') <- d
  
  
  return (biom)
}



#' Inspects a file to see if it's most likely hdf5, json, or tab-delimited.
#' Also reports gzip or bzip2 compression.
#' 
#' @name guess_format
#' 
#' @param fp  The path to a BIOM file.
#' 
#' @return One of \code{c("tsv", "tsv.gz", "tsv.bz2", "json", "json.gz", 
#'         "json.bz2", "hdf5")}.
#' @export
#' @examples
#'     
#'     fp <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     guess_format(fp)
#'
guess_format <- function (fp) {
  
  stopifnot(file.exists(fp))
  
  con <- file(fp)
  on.exit(close(con))
  
  format <- local({
    if (identical(readChar(con, 1L), '{'))       return ("json")
    if (identical(readChar(con, 4L), "\x89HDF")) return ("hdf5")
    return ("tsv")
  })
  
  format <- switch(
    EXPR = summary(con)$class,
    'gzfile' = paste0(format, ".gz"),
    'bzfile' = paste0(format, ".bz2"),
    format)
  
  return (format)
}



read_biom_tsv <- function (fp) {

  #________________________________________________________
  # Read in all the lines from the file
  #________________________________________________________

  lines <- try(readLines(fp, warn=FALSE), silent=TRUE)
  if (is(lines, "try-error"))
    stop(simpleError(sprintf("Unable to parse tab delimited file. %s", as.character(lines))))


  #________________________________________________________
  # Write to a temp file all the lines that have content and don't begin with '#'
  #________________________________________________________

  lines <- trimws(lines)
  lines <- c(
    grep("^#SampleID", lines, value=TRUE),             # Keep
    grep("^#OTU",      lines, value=TRUE),             # Keep
    grep("^#",         lines, value=TRUE, invert=TRUE) # Discard
  )
  fp <- tempfile()
  writeLines(lines, fp, "\n")


  #________________________________________________________
  # Comma or a tab first? To infer csv vs tsv format.
  #________________________________________________________

  csv <- min(gregexpr(",",   lines[[1]])[[1]], fixed=TRUE)
  tsv <- min(gregexpr("\\t", lines[[1]])[[1]], fixed=TRUE)

  if (csv > -1 && tsv > -1) { importFn <- if (csv < tsv) read.csv else read.delim
  } else if (csv > -1)      { importFn <- read.csv
  } else if (tsv > -1)      { importFn <- read.delim
  } else                    { stop(simpleError("Cannot find comma or tab delimiters in file.")) }

  mat <- try(importFn(fp, header=FALSE, colClasses="character"), silent=TRUE)
  if (is(mat, "try-error"))
    stop(simpleError(sprintf("Error parsing file: %s", as.character(mat))))

  mat <- try(as.matrix(mat), silent=TRUE)
  if (is(mat, "try-error"))
    stop(simpleError(sprintf("Error converting to matrix: %s", as.character(mat))))


  #________________________________________________________
  # Ensure we have at least one taxa (row headers) and one sample (column headers)
  #________________________________________________________

  if (nrow(mat) < 2 || ncol(mat) < 2) {
    msg <- "Unable to parse provided file: found %i rows and %i columns"
    stop(simpleError(sprintf(msg, nrow(mat), ncol(mat))))
  }
  
  
  #________________________________________________________
  # Check for duplicate taxa or sample names
  #________________________________________________________
  
  dupTaxaNames   <- unique(mat[duplicated(mat[,1]),1])
  dupSampleNames <- unique(mat[1,duplicated(mat[1,])])
  if (length(dupTaxaNames)   > 5) dupTaxaNames[[5]]   <- sprintf("... +%i more", length(dupTaxaNames)   - 4)
  if (length(dupSampleNames) > 5) dupSampleNames[[5]] <- sprintf("... +%i more", length(dupSampleNames) - 4)
  if (length(dupTaxaNames)   > 0) stop(simpleError(sprintf("Duplicate taxa names: %s", paste(collapse=",", dupTaxaNames))))
  if (length(dupSampleNames) > 0) stop(simpleError(sprintf("Duplicate sample names: %s", paste(collapse=",", dupSampleNames))))
  
  
  #________________________________________________________
  # Move the Taxa and Sample names into the matrix headers
  #________________________________________________________
  
  dimnames(mat) <- list(mat[,1], mat[1,])
  mat <- mat[-1,-1, drop=FALSE]

  return (mat)
}

read_biom_json <- function (fp) {

  json <- try(jsonlite::read_json(path=fp), silent=TRUE)
  if (is(json, "try-error"))
    stop(simpleError(sprintf("Unable to parse JSON file. %s", as.character(json))))

  if (!all(c('data', 'matrix_type', 'rows', 'columns') %in% names(json)))
      stop(simpleError("BIOM file requires: data, matrix_type, shape, rows, columns"))

  return (json)
}

read_biom_hdf5 <- function (fp) {
  
  hdf5 <- try(rhdf5::H5Fopen(fp), silent=TRUE)
  if (is(hdf5, "try-error")) {
    try(rhdf5::H5Fclose(hdf5), silent=TRUE)
    stop(simpleError(sprintf("Unable to parse HDF5 file. %s", as.character(hdf5))))
  }
  
  entries  <- with(rhdf5::h5ls(hdf5), paste(sep="/", group, name))
  expected <- c( "/observation/matrix/indptr", "/observation/matrix/indices", 
                 "/observation/ids", "/observation/matrix/data", "/sample/ids" )
  
  missing  <- setdiff(expected, entries)
  if (length(missing) > 0) {
    try(rhdf5::H5Fclose(hdf5), silent=TRUE)
    stop(simpleError(sprintf("BIOM file requires: %s", paste(collapse=",", missing))))
  }
  
  return (hdf5)
}



parse_json_metadata <- function (json) {
  
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

parse_hdf5_metadata <- function (hdf5) {

  keys <- names(hdf5$sample$metadata)
  vals <- lapply(keys, function (k) {
    
    obj_class <- typeof(hdf5$sample$metadata[[k]])
    
    if (identical(obj_class, "character")) {
      as(hdf5$sample$metadata[[k]], "character")
    } else {
      as(hdf5$sample$metadata[[k]], "numeric")
    }
  })

  as.data.frame(
    row.names        = hdf5$sample$ids,
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

parse_hdf5_info <- function (hdf5) {
  attrs <- rhdf5::h5readAttributes(hdf5, "/")
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
    stop(simpleError("No abundance counts."))

  slam::as.simple_triplet_matrix(mtx)
}

parse_json_counts <- function (json) {

  if (!any(c('sparse', 'dense') %in% json$matrix_type))
    stop(simpleError("BIOM file's matrix_type must be either 'sparse' or 'dense'"))

  if (length(json$data) == 0)
    stop(simpleError("BIOM file does not have any count data."))

  TaxaIDs   <- sapply(json$rows,    function (x) unlist(x$id))
  SampleIDs <- sapply(json$columns, function (x) unlist(x$id))

  if (json$matrix_type == "sparse")
    counts <- slam::simple_triplet_matrix(
                    i        = sapply(json$data, function (x) x[[1]]) + 1,
                    j        = sapply(json$data, function (x) x[[2]]) + 1,
                    v        = sapply(json$data, function (x) x[[3]]),
                    nrow     = length(TaxaIDs),
                    ncol     = length(SampleIDs),
                    dimnames = list(TaxaIDs, SampleIDs))

  if (json$matrix_type == "dense")
    counts <- slam::as.simple_triplet_matrix(
                    matrix(
                      data     = unlist(json$data),
                      byrow    = TRUE,
                      nrow     = length(TaxaIDs),
                      ncol     = length(SampleIDs),
                      dimnames = list(TaxaIDs, SampleIDs)))

  return (counts)
}

parse_hdf5_counts <- function (hdf5) {

  indptr <- as.numeric(hdf5$observation$matrix$indptr)

  slam::simple_triplet_matrix(
        i        = unlist(sapply(1:(length(indptr)-1), function (i) rep(i, diff(indptr[c(i,i+1)])))),
        j        = as.numeric(hdf5$observation$matrix$indices) + 1,
        v        = as.numeric(hdf5$observation$matrix$data),
        nrow     = length(hdf5$observation$ids),
        ncol     = length(hdf5$sample$ids),
        dimnames = list(
          as.character(hdf5$observation$ids),
          as.character(hdf5$sample$ids)
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
  if (identical(unname(taxa_table[,1]), rownames(taxa_table)))
    taxa_table <- taxa_table[,-1,drop=FALSE]

  return (taxa_table)
}

parse_json_taxonomy <- function (json) {

  taxa_table <- sapply(json$rows, simplify = "array", function (x) {
    
    taxaNames <- unlist(x$metadata$taxonomy)
    
    # The taxa names aren't where they're supposed to be
    if (is_null(taxaNames)) {
      
      # MicrobiomeDB puts the taxa string in the 'ID' field, e.g, {"metadata":null,"id":"Archaea;Euryarchaeota;...
      if (identical(json[['generated_by']], "MicrobiomeDB")) {
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

  if (is(taxa_table, "matrix")) {
    # 'Normal' situation of same number of ranks (>1) per taxa string.
    
    taxa_table <- t(taxa_table)
    
  } else if (is(taxa_table, "list")) {
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

parse_hdf5_taxonomy <- function (hdf5) {
  
  if ("taxonomy" %in% names(hdf5$observation$metadata)) {
    taxa_table           <- t(hdf5$observation$metadata$taxonomy)
    rownames(taxa_table) <- as.character(hdf5$observation$ids)
    colnames(taxa_table) <- default_taxa_ranks(ncol(taxa_table))
    
  } else {
    ids        <- as.character(hdf5$observation$ids)
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

parse_hdf5_sequences <- function (hdf5) {
  
  # No sequence information
  if (!"sequences" %in% names(hdf5$observation$metadata))
    return (NULL)
  
  ids  <- as.character(hdf5$observation$ids)
  seqs <- as.character(hdf5$observation$metadata$sequences)
  
  res <- setNames(seqs, ids)
    
  return (res)
}



parse_json_tree <- function (json, tree_mode) {
  
  # Obey the tree argument
  #________________________________________________________
  if (identical(tree_mode, TRUE)) {
    tree <- unlist(json[['phylogeny']])
    
  } else if (identical(tree_mode, FALSE)) {
    return (NULL)
    
  } else if (identical(tree_mode, 'auto')) {
    tree <- unlist(json[['phylogeny']])
    if (is_null(tree))       return (NULL)
    if (is.na(tree))         return (NULL)
    if (!is.character(tree)) return (NULL)
    if (!length(tree) == 1)  return (NULL)
    if (!nchar(tree) >= 1)   return (NULL)
    
  } else {
    tree <- tree_mode
  }
  

  # Try to read it, assuming newick format
  #________________________________________________________
  tree <- try(rbiom::read_tree(tree), silent=TRUE)
  if (is(tree, "try-error")) {
    errmsg <- sprintf("Unable to read embedded phylogeny. %s", as.character(tree))
    if (!identical(tree_mode, 'auto')) stop(errmsg)
    cat(file=stderr(), errmsg)
    return (NULL)
  }


  # Make sure it has all the OTUs from the table
  #________________________________________________________
  TaxaIDs <- sapply(json$rows, function (x) unlist(x$id))
  missing <- setdiff(TaxaIDs, tree$tip.label)
  if (length(missing) > 0) {
    if (length(missing) > 6)
      missing <- c(missing[1:5], sprintf("+ %i more", length(missing) - 5))
    errmsg <- sprintf("OTUs missing from tree: %s\n", paste(collapse=",", missing))
    if (!identical(tree_mode, 'auto')) stop(errmsg)
    cat(file=stderr(), errmsg)
    return (NULL)
  }
  
  
  # Drop any extra taxa found in the tree
  #________________________________________________________
  if (length(tree$tip.label) > length(TaxaIDs))
    tree <- rbiom::subtree(tree, TaxaIDs)


  return (tree)
}

parse_hdf5_tree <- function (hdf5, tree_mode) {
  
  # Obey the tree argument
  #________________________________________________________
  if (identical(tree_mode, FALSE)) {
    return (NULL)
  
  } else if (identical(tree_mode, TRUE) || identical(tree_mode, 'auto')) {
    
    # See if a tree is included in the BIOM file
    #________________________________________________________
    tree <- as.character(hdf5$observation$`group-metadata`$phylogeny)
    
    errmsg <- NULL
    
    if (is_null(tree) || identical(tree, character(0))) {
      errmsg <- "There is no tree in this BIOM file."
      
    } else if (!length(tree)) {
      errmsg <- "There is no tree in this BIOM file."
      
    } else {
    
      # Assume it's newick format unless otherwise indicated
      #________________________________________________________
      attrs <- rhdf5::h5readAttributes(hdf5, "observation/group-metadata/phylogeny")
      if ("data_type" %in% names(attrs)) {
        data_type <- tolower(as.character(attrs[['data_type']]))
        if (!identical(data_type, "newick"))
          errmsg <- sprintf("Phylogeny is not Newick format, is '%s'.", data_type)
      }
    }
    
    if (!is_null(errmsg)) {
      if (identical(tree_mode, TRUE)) stop(errmsg)
      return (NULL)
    }
    
  } else {
    tree <- tree_mode
  }
  
  
  # Try to read the Newick-formatted tree
  #________________________________________________________
  tree <- try(rbiom::read_tree(tree), silent=TRUE)
  if (is(tree, "try-error")) {
    errmsg <- sprintf("Unable to read embedded phylogeny. %s", as.character(tree))
    if (!identical(tree_mode, 'auto')) stop(errmsg)
    cat(file=stderr(), errmsg)
    return (NULL)
  }
  
  
  # Make sure it has all the OTUs from the table
  #________________________________________________________
  TaxaIDs <- as.character(hdf5$observation$ids)
  missing <- setdiff(TaxaIDs, tree$tip.label)
  if (length(missing) > 0) {
    if (length(missing) > 6)
      missing <- c(missing[1:5], sprintf("+ %i more", length(missing) - 5))
    errmsg <- sprintf("OTUs missing from tree: %s\n", paste(collapse=",", missing))
    if (!identical(tree_mode, 'auto')) stop(errmsg)
    cat(file=stderr(), errmsg)
    return (NULL)
  }
  
  
  # Drop any extra taxa found in the tree
  #________________________________________________________
  if (length(tree$tip.label) > length(TaxaIDs))
    tree <- rbiom::subtree(tree, TaxaIDs)
  
  
  return (tree)
}




default_taxa_ranks <- function (n) {
  if (n >= 6 && n <= 8)
    return(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:n])
  sprintf("Level.%i", 1:n)
}




