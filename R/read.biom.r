#' Extracts counts, metadata, taxonomy, and phylogeny from a biom file.
#'
#' @param src  Input data as either a file path, URL, or JSON string.
#'     \code{read.biom} can read BIOM files formatted according to both the
#'     version 1.0 (JSON) and 2.1 (HDF5)
#'     \href{http://biom-format.org/documentation/}{specifications} as well as
#'     classical tabular format. URLs must begin with \kbd{http://},
#'     \kbd{https://}, \kbd{ftp://}, or \kbd{ftps://}. JSON files must have
#'     \code{\{} as their first non-whitespace character. Compressed (gzip or
#'     bzip2) biom files are also supported.
#' @param tree  The default value of \code{auto} will read the tree from the 
#'     biom file specified in \code{src}, if present. The value \code{TRUE} 
#'     will do the same, but will generate an error message if a tree is not 
#'     present. Setting \code{tree=FALSE} will return a \code{BIOM} object 
#'     without any tree data. You may also provide a file path, URL, or Newick
#'     string to load that tree data into the final \code{BIOM} object.
#' @param prune  Should samples and taxa with zero observations be discarded?
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
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
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


read.biom <- function (src, tree='auto', prune=FALSE) {

  #--------------------------------------------------------------
  # Sanity check input values
  #--------------------------------------------------------------

  if (length(src) != 1 | !is(src, "character"))
    stop(simpleError("Data source for read.biom() must be a single string."))


  #--------------------------------------------------------------
  # Get the url or text for the src/BIOM data into a file
  #--------------------------------------------------------------

  if (length(grep("^(ht|f)tps{0,1}://.+", src)) == 1) {

    fp <- tempfile(fileext=basename(src))
    on.exit(unlink(fp), add=TRUE)
    if (!identical(0L, try(download.file(src, fp, quiet=TRUE), silent=TRUE)))
        stop(simpleError(sprintf("Cannot retrieve URL %s", src)))

  } else if (length(grep("^[ \t\n]*\\{", src)) == 1) {

    fp <- tempfile(fileext=".biom")
    on.exit(unlink(fp), add=TRUE)
    if (!is.null(try(writeChar(src, fp), silent=TRUE)))
        stop(simpleError(sprintf("Cannot write text to file %s", fp)))

  } else {

    fp <- normalizePath(src)
  }

  if (!file.exists(fp))
    stop(simpleError(sprintf("Cannot locate file %s", fp)))


  #--------------------------------------------------------------
  # Uncompress files that are in gzip or bzip2 format
  #--------------------------------------------------------------

  file_con   <- file(fp)
  file_class <- summary(file_con)$class
  close.connection(file_con)
  
  if (file_class %in% c("gzfile", "bzfile")) {
    
    if (identical(file_class, "gzfile"))
      fp <- R.utils::gunzip(fp, destname=tempfile(), remove=FALSE)
    
    if (identical(file_class, "bzfile"))
      fp <- R.utils::bunzip2(fp, destname=tempfile(), remove=FALSE)
    
    on.exit(unlink(fp), add=TRUE)
  }
  
  remove("file_con", "file_class")
  
  
  #--------------------------------------------------------------
  # Process the file according to its internal format
  #--------------------------------------------------------------
  
  if (rhdf5::H5Fis_hdf5(fp)) {
    
    #-=-=-=-=-=-=-=-=-=-#
    # HDF5 file format  #
    #-=-=-=-=-=-=-=-=-=-#
    
    hdf5      <- PB.HDF5.ReadHDF5(fp)
    counts    <- PB.HDF5.Counts(hdf5)
    sequences <- PB.HDF5.Sequences(hdf5)
    taxonomy  <- PB.HDF5.Taxonomy(hdf5)
    metadata  <- PB.HDF5.Metadata(hdf5)
    info      <- PB.HDF5.Info(hdf5)
    phylogeny <- PB.HDF5.Tree(hdf5, tree)
    
    rhdf5::H5Fclose(hdf5)
    remove("hdf5")
    
  } else if (identical("{", readChar(fp, 1))) {
    
    #-=-=-=-=-=-=-=-=-=-#
    # JSON file format  #
    #-=-=-=-=-=-=-=-=-=-#
    
    json      <- PB.JSON.ReadJSON(fp)
    counts    <- PB.JSON.Counts(json)
    sequences <- PB.JSON.Sequences(json)
    taxonomy  <- PB.JSON.Taxonomy(json)
    metadata  <- PB.JSON.Metadata(json)
    info      <- PB.JSON.Info(json)
    phylogeny <- PB.JSON.Tree(json, tree)
    
    remove("json")
    
  } else {
    
    #-=-=-=-=-=-=-=-=-=-#
    # TSV file format   #
    #-=-=-=-=-=-=-=-=-=-#
    
    if (identical(tree, TRUE))
      stop(simpleError("It is impossible to load a phylogenetic tree from a BIOM file in tab-separated format."))
    
    mtx       <- PB.TSV.ReadTSV(fp)
    counts    <- PB.TSV.Counts(mtx)
    sequences <- NULL
    taxonomy  <- PB.TSV.Taxonomy(mtx)
    metadata  <- data.frame(row.names=colnames(counts))
    info      <- list(id=tools::md5sum(fp)[[1]], type="OTU table")
    phylogeny <- NULL
  }
  
  
  #--------------------------------------------------------------
  # Discard samples/taxa with zero observations
  #--------------------------------------------------------------
  
  if (identical(prune, TRUE)) {
    counts   <- counts[slam::row_sums(counts) > 0, slam::col_sums(counts) > 0]
    taxonomy <- taxonomy[rownames(counts),,drop=FALSE]
    metadata <- metadata[colnames(counts),,drop=FALSE]
    if (!is.null(sequences))
      sequences <- sequences[rownames(counts)]
    if (!is.null(phylogeny))
      phylogeny <- subtree(phylogeny, rownames(counts))
  }
  
  
  #--------------------------------------------------------------
  # Return everything we've computed as a BIOM class object
  #--------------------------------------------------------------
  
  structure(
    class = c("BIOM", "list"),
    list( 'counts'    = counts,
          'metadata'  = metadata,
          'taxonomy'  = taxonomy,
          'phylogeny' = phylogeny,
          'sequences' = sequences,
          'info'      = info
    )
  )


}



PB.TSV.ReadTSV <- function (fp) {

  #--------------------------------------------------------------
  # Read in all the lines from the file
  #--------------------------------------------------------------

  lines <- try(readLines(fp, warn=FALSE), silent=TRUE)
  if (is(lines, "try-error"))
    stop(simpleError(sprintf("Unable to parse tab delimited file. %s", as.character(lines))))


  #--------------------------------------------------------------
  # Write to a temp file all the lines that have content and don't begin with '#'
  #--------------------------------------------------------------

  lines <- trimws(lines)
  lines <- c(
    grep("^#SampleID", lines, value=TRUE),             # Keep
    grep("^#OTU",      lines, value=TRUE),             # Keep
    grep("^#",         lines, value=TRUE, invert=TRUE) # Discard
  )
  fp <- tempfile()
  writeLines(lines, fp, "\n")


  #--------------------------------------------------------------
  # See if we see a comma or a tab first, thereby deducing csv or tsv format
  #--------------------------------------------------------------

  csv <- min(gregexpr(",",   lines[[1]])[[1]], fixed=TRUE)
  tsv <- min(gregexpr("\\t", lines[[1]])[[1]], fixed=TRUE)

  if (csv > -1 & tsv > -1) { importFn <- if (csv < tsv) read.csv else read.delim
  } else if (csv > -1)     { importFn <- read.csv
  } else if (tsv > -1)     { importFn <- read.delim
  } else                   { stop(simpleError("Cannot find comma or tab delimiters in file.")) }

  mat <- try(importFn(fp, header=FALSE, colClasses="character"), silent=TRUE)
  if (is(mat, "try-error"))
    stop(simpleError(sprintf("Error parsing file: %s", as.character(mat))))

  mat <- try(as.matrix(mat), silent=TRUE)
  if (is(mat, "try-error"))
    stop(simpleError(sprintf("Error converting to matrix: %s", as.character(mat))))


  #--------------------------------------------------------------
  # Ensure we have at least one taxa (row headers) and one sample (column headers)
  #--------------------------------------------------------------

  if (nrow(mat) < 2 | ncol(mat) < 2) {
    msg <- "Unable to parse provided file: found %i rows and %i columns"
    stop(simpleError(sprintf(msg, nrow(mat), ncol(mat))))
  }
  
  
  #--------------------------------------------------------------
  # Check for duplicate taxa or sample names
  #--------------------------------------------------------------
  
  dupTaxaNames   <- unique(mat[duplicated(mat[,1]),1])
  dupSampleNames <- unique(mat[1,duplicated(mat[1,])])
  if (length(dupTaxaNames)   > 5) dupTaxaNames[[5]]   <- sprintf("... +%i more", length(dupTaxaNames)   - 4)
  if (length(dupSampleNames) > 5) dupSampleNames[[5]] <- sprintf("... +%i more", length(dupSampleNames) - 4)
  if (length(dupTaxaNames)   > 0) stop(simpleError(sprintf("Duplicate taxa names: %s", paste(collapse=",", dupTaxaNames))))
  if (length(dupSampleNames) > 0) stop(simpleError(sprintf("Duplicate sample names: %s", paste(collapse=",", dupSampleNames))))
  
  
  #--------------------------------------------------------------
  # Move the Taxa and Sample names into the matrix headers
  #--------------------------------------------------------------
  
  dimnames(mat) <- list(mat[,1], mat[1,])
  mat <- mat[-1,-1, drop=FALSE]

  return (mat)
}

PB.JSON.ReadJSON <- function (fp) {

  json <- try(rjson::fromJSON(file=fp), silent=TRUE)
  if (is(json, "try-error"))
    stop(simpleError(sprintf("Unable to parse JSON file. %s", as.character(json))))

  if (!all(c('data', 'matrix_type', 'rows', 'columns') %in% names(json)))
      stop(simpleError("BIOM file requires: data, matrix_type, shape, rows, columns"))

  return (json)
}

PB.HDF5.ReadHDF5 <- function (fp) {
  
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



PB.JSON.Metadata <- function (json) {

  df <- as.data.frame(
    row.names        = sapply(json$columns, function (x) x[['id']]),
    check.names      = FALSE,
    stringsAsFactors = FALSE,
    vapply(names(json$columns[[1]]$metadata), function (i) {
      sapply(json$columns, function (x) as.character(x[['metadata']][[i]]))
    }, character(length(json$columns)), USE.NAMES=TRUE)
  )

  # Convert strings to numbers
  for (i in seq_len(ncol(df))) {
    df[[i]][which(df[[i]] == "NA")] <- NA
    if (all(grepl("^-?(0|[1-9]\\d*)(\\.\\d*[1-9]|)(e[\\+\\-]{0,1}[0-9]+|)$", na.omit(df[[i]]))))
      df[[i]] <- as.numeric(df[[i]])
  }

  return (df)
}

PB.HDF5.Metadata <- function (hdf5) {

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



PB.JSON.Info <- function (json) {
  fields <- sort(names(json))
  fields <- fields[!fields %in% c("rows", "columns", "data", "phylogeny")]
  json[fields]
}

PB.HDF5.Info <- function (hdf5) {
  attrs <- rhdf5::h5readAttributes(hdf5, "/")
  for (i in names(attrs)) {
    attrs[[i]] <- as(attrs[[i]], typeof(attrs[[i]]))
  }
  return (attrs)
}



PB.TSV.Counts <- function (mtx) {

  # Only keep columns that are all numbers
  allNumbers <- function (x) all(grepl("^\\d+(\\.\\d+|)([Ee][\\+\\-]{0,1}[0-9]+|)$", x))
  mtx <- mtx[, apply(mtx, 2L, allNumbers), drop=FALSE]
  mtx <- matrix(
    data     = as.numeric(mtx),
    nrow     = nrow(mtx),
    dimnames = list(rownames(mtx), colnames(mtx)) )

  if (length(mtx) == 0 | sum(mtx) == 0)
    stop(simpleError("No abundance counts."))

  slam::as.simple_triplet_matrix(mtx)
}

PB.JSON.Counts <- function (json) {

  if (!any(c('sparse', 'dense') %in% json$matrix_type))
    stop(simpleError("BIOM file's matrix_type must be either 'sparse' or 'dense'"))

  if (length(json$data) == 0)
    stop(simpleError("BIOM file does not have any count data."))

  TaxaIDs   <- sapply(json$rows,    function (x) x$id)
  SampleIDs <- sapply(json$columns, function (x) x$id)

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

PB.HDF5.Counts <- function (hdf5) {

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



PB.TSV.Taxonomy <- function (mtx) {

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
  colnames(taxa_table) <- PB.TaxaLevelNames(ncol(taxa_table))

  # Better to return a completely empty table than one with just the taxa IDs
  if (identical(unname(taxa_table[,1]), rownames(taxa_table)))
    taxa_table <- taxa_table[,-1,drop=FALSE]

  return (taxa_table)
}

PB.JSON.Taxonomy <- function (json) {

  # No taxonomic information
  if (!any(sapply(json$rows, function (x) 'taxonomy' %in% names(x$metadata) )))
    return (NULL)


  taxa_table <- sapply(json$rows, simplify = "array", function (x) {
    taxaNames <- x$metadata$taxonomy
    
    if (is.null(taxaNames))
      return (x$id)
    
    if (length(taxaNames) == 1) {
      if (nchar(taxaNames) == 0) return (x$id)
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

  rownames(taxa_table) <- sapply(json$rows, function (x) x$id)
  colnames(taxa_table) <- PB.TaxaLevelNames(ncol(taxa_table))

  #taxa_table <- PB.SanitizeTaxonomy(taxa_table, env=parent.frame())

  return (taxa_table)
}

PB.HDF5.Taxonomy <- function (hdf5) {
  
  if ("taxonomy" %in% names(hdf5$observation$metadata)) {
    taxa_table           <- t(hdf5$observation$metadata$taxonomy)
    rownames(taxa_table) <- as.character(hdf5$observation$ids)
    colnames(taxa_table) <- PB.TaxaLevelNames(ncol(taxa_table))
    
  } else {
    ids        <- as.character(hdf5$observation$ids)
    taxa_table <- matrix(nrow=length(ids), ncol=0, dimnames=list(ids, character(0)))
  }

  #taxa_table <- PB.SanitizeTaxonomy(taxa_table, env=parent.frame())
    
  return (taxa_table)
}




PB.JSON.Sequences <- function (json) {
  
  # No sequence information
  if (!any(sapply(json$rows, function (x) 'sequence' %in% names(x$metadata) )))
    return (NULL)
  
  ids  <- sapply(json$rows, simplify = "array", function (x) { x$id })
  seqs <- sapply(json$rows, simplify = "array", function (x) {
    if (is.null(x$metadata$sequence)) NA else x$metadata$sequence
  })
  
  res <- setNames(seqs, ids)
  
  return (res)
}

PB.HDF5.Sequences <- function (hdf5) {
  
  # No sequence information
  if (!"sequences" %in% names(hdf5$observation$metadata))
    return (NULL)
  
  ids  <- as.character(hdf5$observation$ids)
  seqs <- as.character(hdf5$observation$metadata$sequences)
  
  res <- setNames(seqs, ids)
    
  return (res)
}



PB.JSON.Tree <- function (json, tree_mode) {
  
  # Obey the tree argument
  #------------------------------------------------------
  if (identical(tree_mode, TRUE)) {
    tree <- json[['phylogeny']]
    
  } else if (identical(tree_mode, FALSE)) {
    return (NULL)
    
  } else if (identical(tree_mode, 'auto')) {
    if (!'phylogeny' %in% names(json))      return (NULL)
    if (is.null(json[['phylogeny']]))       return (NULL)
    if (is.na(json[['phylogeny']]))         return (NULL)
    if (!is.character(json[['phylogeny']])) return (NULL)
    if (!length(json[['phylogeny']]) == 1)  return (NULL)
    if (!nchar(json[['phylogeny']]) >= 1)   return (NULL)
    tree <- json[['phylogeny']]
    
  } else {
    tree <- tree_mode
  }
  

  # Try to read it, assuming newick format
  #------------------------------------------------------
  tree <- try(rbiom::read.tree(tree), silent=TRUE)
  if (is(tree, "try-error")) {
    errmsg <- sprintf("Unable to read embedded phylogeny. %s", as.character(tree))
    if (!identical(tree_mode, 'auto')) stop(errmsg)
    cat(file=stderr(), errmsg)
    return (NULL)
  }


  # Make sure it has all the OTUs from the table
  #------------------------------------------------------
  TaxaIDs <- sapply(json$rows, function (x) x$id)
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
  #------------------------------------------------------
  if (length(tree$tip.label) > length(TaxaIDs))
    tree <- rbiom::subtree(tree, TaxaIDs)


  return (tree)
}

PB.HDF5.Tree <- function (hdf5, tree_mode) {
  
  # Obey the tree argument
  #------------------------------------------------------
  if (identical(tree_mode, FALSE)) {
    return (NULL)
  
  } else if (identical(tree_mode, TRUE) | identical(tree_mode, 'auto')) {
    
    # See if a tree is included in the BIOM file
    #------------------------------------------------------
    tree <- as.character(hdf5$observation$`group-metadata`$phylogeny)
    
    errmsg <- NULL
    
    if (is.null(tree) || identical(tree, character(0))) {
      errmsg <- "There is no tree in this BIOM file."
      
    } else if (!length(tree)) {
      errmsg <- "There is no tree in this BIOM file."
      
    } else {
    
      # Assume it's newick format unless otherwise indicated
      #------------------------------------------------------
      attrs <- rhdf5::h5readAttributes(hdf5, "observation/group-metadata/phylogeny")
      if ("data_type" %in% names(attrs)) {
        data_type <- tolower(as.character(attrs[['data_type']]))
        if (!identical(data_type, "newick"))
          errmsg <- sprintf("Phylogeny is not Newick format, is '%s'.", data_type)
      }
    }
    
    if (!is.null(errmsg) & identical(tree_mode, TRUE)) {
      stop(errmsg)
    }
    
  } else {
    tree <- tree_mode
  }
  
  
  # Try to read the Newick-formatted tree
  #------------------------------------------------------
  tree <- try(rbiom::read.tree(tree), silent=TRUE)
  if (is(tree, "try-error")) {
    errmsg <- sprintf("Unable to read embedded phylogeny. %s", as.character(tree))
    if (!identical(tree_mode, 'auto')) stop(errmsg)
    cat(file=stderr(), errmsg)
    return (NULL)
  }
  
  
  # Make sure it has all the OTUs from the table
  #------------------------------------------------------
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
  #------------------------------------------------------
  if (length(tree$tip.label) > length(TaxaIDs))
    tree <- rbiom::subtree(tree, TaxaIDs)
  
  
  return (tree)
}




PB.TaxaLevelNames <- function (n) {
  if (n >= 6 && n <= 8)
    return(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:n])
  sprintf("Level.%i", 1:n)
}




# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Rename poorly described taxa
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Definition of poorly named taxa:
# #  - ambiguious keywords
# #  - numbers outside of parentheses
# #  - length of less than 2 characters
# #  - two capital letters in a row
# #------------------------------------------------------------------
#
# PB.SanitizeTaxonomy <- function (tt, env=NULL) {
#
#   if (is.null(env)) env <- parent.frame()
#
#   tt <- sub("^.*?__\ *", "", tt)
#   tt <- sub("^[a-z]_",   "", tt)
#   tt <- sub("^[\ _]+",   "", tt)
#   tt <- sub("[\ ]+$",    "", tt)
#   tt[which(tt == "")] <- NA
#
#   invalid_case_sensi <- "([A-Z]{2})"
#   invalid_case_insen <- paste(sep="|",
#     "unknown",    "uncultured",      "unclassified",    "unidentified",
#     "group",      "subsection",      "family",          "lineage",
#     "candidate",  "incertae.*sedis", "(^.$)" )
#   exception_rewrites <- c(
#     "TM7"                  = "Saccharibacteria (TM7)",
#     "RF9"                  = "Mollicutes RF9",
#     "Escherichia_Shigella" = "Escherichia/Shigella")
#
#   # Search and replace bad names with NA
#   taxaNames <- unique(as.vector(tt))
#   invalids1 <- grep(invalid_case_insen, taxaNames, value=TRUE, ignore.case=TRUE)
#   invalids2 <- grep(invalid_case_sensi, taxaNames, value=TRUE, ignore.case=FALSE)
#   invalids  <- unique(c(invalids1, invalids2))
#   invalids  <- setNames(rep(NA, length(invalids)), invalids)
#   rewrites  <- c(exception_rewrites, invalids)
#   rewrites  <- rewrites[!duplicated(names(rewrites))]
#   tt[]      <- plyr::revalue(as.vector(as.matrix(tt)), rewrites, FALSE)
#
#   tt <- cbind(tt, rownames(tt))
#
#   # Replace NA with "<Last prior Non-NA value> (<Next Non-NA value or OTU Name>)"
#   env$tt <- tt
#
#   tt <- withCluster(nTasks=nrow(tt), env=env, {
#
#     foreach::`%dopar%`(
#       foreach::foreach(set=sets, .combine='cbind', .options.snow=opts),
#       apply(tt[set,,drop=FALSE], 1, function (row) {
#         cols <- which(is.na(row))
#         if (length(cols) > 0)
#           for (i in split(seq_along(cols), cumsum(c(0, diff(cols) > 1))))
#             row[cols[i]] <- paste0(row[cols[min(i)] - 1], " (", row[cols[max(i)] + 1], ")")
#         return (row)
#       })
#     )
#   })
#
#   return (t(tt[1:(nrow(tt) - 1),,drop=FALSE]))
# }
