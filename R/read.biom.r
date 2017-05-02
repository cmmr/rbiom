#' Extracts counts, metadata, taxonomy, and phylogeny from a biom file.
#'
#' @param src  Input data as either a file path, URL, or JSON string.
#'     \code{read.biom} can read BIOM files formatted according to both the
#'     version 1.0 (JSON) and 2.1 (HDF5)
#'     \href{http://biom-format.org/documentation/}{specifications} as well as
#'     classical tabular format. URLs must begin with \kbd{http://}, 
#'     \kbd{https://}, \kbd{ftp://}, or \kbd{ftps://}. JSON files must have
#'     \code{\{} as their first non-whitespace character.
#' @param progressbar  Whether to display a progress bar and status messages
#'     (logical). Will automatically tie in with \pkg{shiny} if run within a
#'     \pkg{shiny} session.
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
#'         \item{info}{A list of other attributes defined in the BIOM file,
#'         such as \code{id}, \code{type}, \code{format}, \code{format_url},
#'         \code{generated_by}, \code{date}, \code{matrix_type},
#'         \code{matrix_element_type}, and \code{shape}}
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
#'     table(biom$metadata$Sex, biom$metadata$Body.Site)
#'     sprintf("Mean age: %.1f", mean(biom$metadata$Age))
#'
#'     # Phylogenetic tree
#'     tree <- biom$phylogeny
#'     top5.tree <- ape::drop.tip(tree, setdiff(tree$tip.label, top5))
#'     plot(top5.tree)
#'


read.biom <- function (src, progressbar=FALSE) {

  #--------------------------------------------------------------
  # Sanity check input values
  #--------------------------------------------------------------

  if (length(src) != 1 | !is(src, "character"))
    stop(simpleError("Data source for read.biom() must be a single string."))


  pb <- progressBar(progressbar=progressbar)
  on.exit(pb$close())


  #--------------------------------------------------------------
  # Get the url or text data into a file
  #--------------------------------------------------------------

  if (length(grep("^(ht|f)tps{0,1}://.+", src)) == 1) {

    pb$set(0, 'Downloading BIOM file')

    file <- tempfile(fileext=basename(src))
    if (!identical(0, try(download.file(src, file, quiet=TRUE), silent=TRUE)))
        stop(simpleError(sprintf("Cannot retrieve URL %s", src)))

  } else if (length(grep("^[ \t\n]*\\{", src)) == 1) {

    file <- tempfile(fileext=".biom")
    if (!is.null(try(writeChar(src, file), silent=TRUE)))
        stop(simpleError(sprintf("Cannot write text to file %s", file)))

  } else {

    file <- normalizePath(src)
  }

  if (!file.exists(file))
    stop(simpleError(sprintf("Cannot locate file %s", file)))



  #--------------------------------------------------------------
  # Process the file according to it's internal format
  #--------------------------------------------------------------

  pb$set(0, 'Determining file type')

  if (h5::is.h5file(file)) {

    #-=-=-=-=-=-=-=-=-=-#
    # HDF5 file format  #
    #-=-=-=-=-=-=-=-=-=-#

    pb$set(0.1, 'Reading HDF5 BIOM file');        hdf5      <- PB.HDF5.ReadHDF5(file)
    pb$set(0.2, 'Assembling OTU table');          counts    <- PB.HDF5.Counts(hdf5)
    pb$set(0.7, 'Processing taxonomic lineages'); taxonomy  <- PB.HDF5.Taxonomy(hdf5)
    pb$set(0.8, 'Extracting metadata');           metadata  <- PB.HDF5.Metadata(hdf5)
    pb$set(0.9, 'Extracting attributes');         info      <- PB.HDF5.Info(hdf5)
    pb$set(1.0, 'Extracting phylogeny');          phylogeny <- PB.HDF5.Tree(hdf5)

    h5::h5close(hdf5)
    remove("hdf5")

  } else if (identical("{", readChar(file, 1))) {

    #-=-=-=-=-=-=-=-=-=-#
    # JSON file format  #
    #-=-=-=-=-=-=-=-=-=-#

    pb$set(0.1, 'Reading JSON BIOM file');        json      <- PB.JSON.ReadJSON(file)
    pb$set(0.2, 'Assembling OTU table');          counts    <- PB.JSON.Counts(json)
    pb$set(0.7, 'Processing taxonomic lineages'); taxonomy  <- PB.JSON.Taxonomy(json)
    pb$set(0.8, 'Extracting metadata');           metadata  <- PB.JSON.Metadata(json)
    pb$set(0.9, 'Extracting attributes');         info      <- PB.JSON.Info(json)
    pb$set(1.0, 'Extracting phylogeny');          phylogeny <- PB.JSON.Tree(json)

    remove("json")
    
  } else {
    
    #-=-=-=-=-=-=-=-=-=-#
    # TSV file format   #
    #-=-=-=-=-=-=-=-=-=-#

    pb$set(0.1, 'Reading tabular data file');     mtx       <- PB.TSV.ReadTSV(file)
    pb$set(0.3, 'Assembling OTU table');          counts    <- PB.TSV.Counts(mtx)
    pb$set(0.8, 'Processing taxonomic lineages'); taxonomy  <- PB.TSV.Taxonomy(mtx)
    pb$set(0.9, 'Extracting metadata');           metadata  <- data.frame(row.names=colnames(counts))
    pb$set(1.0, 'Extracting attributes');         info      <- list()
    pb$set(1.0, 'Extracting phylogeny');          phylogeny <- NULL
    
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
          'info'      = info
    )
  )


}



PB.TSV.ReadTSV <- function (fp) {

  lines <- try(readLines(fp, warn=FALSE), silent=TRUE)
  if (is(lines, "try-error"))
    stop(simpleError(sprintf("Unable to parse JSON file. %s", as.character(lines))))
  
  # Split on tabs
  lines <- sapply(lines, strsplit, split="\t", fixed=TRUE, USE.NAMES=FALSE)
  
  # Use the last line to deduce the number of columns that aren't headers
  lines <- lines[which(sapply(lines, length) == length(lines[[length(lines)]]))]
  
  # We should now have at a minimum:
  #   - First row = Sample Names
  #   - Second row = counts for the first taxa/otu
  #   - First column = Taxa Names
  #   - Second column = counts for the first sample
  
  if (length(lines) < 2 | length(lines[[1]]) < 2)
    stop(simpleError("Unable to parse provided file: unknown type."))
  
  # Cram the tsv lines into a character matrix
  mtx <- matrix(unlist(lines), nrow=length(lines), byrow=TRUE)
  
  if (any(duplicated(mtx[,1]))) stop(simpleError("Duplicate taxa names."))
  if (any(duplicated(mtx[1,]))) stop(simpleError("Duplicate sample names."))
  
  dimnames(mtx) <- list(mtx[,1], mtx[1,])
  mtx <- mtx[-1,-1, drop=FALSE]
  
  return (mtx)
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

  hdf5 <- try(h5::h5file(name=fp, mode="r"), silent=TRUE)
  if (is(hdf5, "try-error"))
    stop(simpleError(sprintf("Unable to parse HDF5 file. %s", as.character(hdf5))))

  expected <- c( "/observation/matrix/indptr", "/observation/matrix/indices", "/observation/ids",
                 "/observation/matrix/data", "/observation/metadata/taxonomy", "/sample/ids" )

  missing  <- setdiff(expected, h5::list.datasets(hdf5))
  if (length(missing) > 0)
      stop(simpleError(sprintf("BIOM file requires: %s", paste(collapse=",", missing))))

  return (hdf5)
}



PB.JSON.Metadata <- function (json) {

  df <- as.data.frame(
    row.names        = sapply(json$columns, function (x) x[['id']]),
    check.names      = FALSE,
    stringsAsFactors = FALSE,
    vapply(names(json$columns[[1]]$metadata), function (i) {
      sapply(json$columns, function (x) x[['metadata']][[i]])
    }, character(length(json$columns)), USE.NAMES=TRUE)
  )
  
  # Convert strings to numbers
  for (i in seq_len(ncol(df))) {
    df[[i]][which(df[[i]] == "NA")] <- NA
    if (all(grepl("^-?(0|[1-9]\\d*)(\\.\\d*[1-9]|)$", na.omit(df[[i]]))))
      df[[i]] <- as.numeric(df[[i]])
  }
  
  return (df)
}

PB.HDF5.Metadata <- function (hdf5) {

  fields <- grep("^/sample/metadata/", h5::list.datasets(hdf5), value=TRUE)

  as.data.frame(
    row.names        = h5::readDataSet(hdf5["/sample/ids"]),
    check.names      = FALSE,
    stringsAsFactors = FALSE,
    setNames(
      lapply(fields, function (f) h5::readDataSet(hdf5[f])),
      sub("/sample/metadata/", "", fields))
  )
}



PB.JSON.Info <- function (json) {
  fields <- sort(names(json))
  fields <- fields[!fields %in% c("rows", "columns", "data", "phylogeny")]
  json[fields]
}

PB.HDF5.Info <- function (hdf5) {
  fields <- sort(h5::list.attributes(hdf5))
  sapply(fields, h5::h5attr, .Object=hdf5, USE.NAMES=TRUE)
}



PB.TSV.Counts <- function (mtx) {
  
  # Only keep column that are all numbers
  allNumbers <- function (x) all(grepl("^\\d+(\\.\\d+|)$", x))
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

  indptr <- h5::readDataSet(hdf5["observation/matrix/indptr"])

  slam::simple_triplet_matrix(
        i        = unlist(sapply(1:(length(indptr)-1), function (i) rep(i, diff(indptr[c(i,i+1)])))),
        j        = h5::readDataSet(hdf5["/observation/matrix/indices"]) + 1,
        v        = h5::readDataSet(hdf5["/observation/matrix/data"]),
        nrow     = length(h5::readDataSet(hdf5["/observation/ids"])),
        ncol     = length(h5::readDataSet(hdf5["/sample/ids"])),
        dimnames = list(
          h5::readDataSet(hdf5["/observation/ids"]),
          h5::readDataSet(hdf5["/sample/ids"])
      ))
}



PB.TSV.Taxonomy <- function (mtx) {
  
  # Discard columns that are all numbers
  allNumbers <- function (x) all(grepl("^\\d+(\\.\\d+|)$", x))
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
    taxaNames <- if (is.null(x$metadata$taxonomy)) x$id else x$metadata$taxonomy
    unname(sapply(taxaNames, function (y) { strsplit(y, "[\\|;]\\ *")[[1]] }))
  })

  # Handle instances where some taxa strings have more levels than others.
  if (class(taxa_table) == "list") {
    n <- sapply(taxa_table, length)
    m <- max(n)
    i <- which(n < m)
    taxa_table[i] <- lapply(taxa_table[i], function (x) c(x, rep(NA, m - length(x))))
    taxa_table    <- t(simplify2array(taxa_table))
  } else {
    taxa_table <- t(taxa_table)
  }

  rownames(taxa_table) <- sapply(json$rows, function (x) x$id)
  colnames(taxa_table) <- PB.TaxaLevelNames(ncol(taxa_table))

  #taxa_table <- PB.SanitizeTaxonomy(taxa_table, env=parent.frame())

  return (taxa_table)
}

PB.HDF5.Taxonomy <- function (hdf5) {

  taxa_table           <- h5::readDataSet(hdf5["/observation/metadata/taxonomy"])
  rownames(taxa_table) <- h5::readDataSet(hdf5["/observation/ids"])
  colnames(taxa_table) <- PB.TaxaLevelNames(ncol(taxa_table))

  #taxa_table <- PB.SanitizeTaxonomy(taxa_table, env=parent.frame())

  return (taxa_table)
}



PB.JSON.Tree <- function (json) {

  # See if a tree is included in the BIOM file
  #------------------------------------------------------
  newick <- json[['phylogeny']]
  if (is.null(newick)) return (NULL)
  if (is.na(newick))   return (NULL)
  if (!nchar(newick))  return (NULL)


  # Try to read it with ape
  #------------------------------------------------------
  tree <- try(ape::read.tree(text=newick), silent=TRUE)
  if (is(tree, "try-error")) {
    errmsg <- sprintf("Unable to read embedded phylogeny. %s", as.character(tree))
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
    cat(file=stderr(), errmsg)
    return (NULL)
  }


  return (tree)
}

PB.HDF5.Tree <- function (hdf5) {

  # See if a tree is included in the BIOM file
  #------------------------------------------------------
  hdf5Path <- "/observation/group-metadata/phylogeny"
  if (!hdf5Path %in% h5::list.datasets(hdf5))
    return (NULL)


  # Assume it's newick format unless otherwise indicated
  #------------------------------------------------------
  data_type <- "newick"
  if ("data_type" %in% h5::list.attributes(hdf5[hdf5Path]))
    data_type <- tolower(h5::h5attr(hdf5[hdf5Path], 'data_type'))

  if (data_type != "newick")
    return (NULL)


  # Try to read it with ape
  #------------------------------------------------------
  tree <- try(ape::read.tree(text=h5::readDataSet(hdf5[hdf5Path])), silent=TRUE)
  if (is(tree, "try-error"))
    stop(simpleError(sprintf("Unable to read phylogeny. %s", as.character(tree))))


  # Make sure it has all the OTUs from the table
  #------------------------------------------------------
  TaxaIDs <- h5::readDataSet(hdf5["/observation/ids"])
  missing <- setdiff(TaxaIDs, tree$tip.label)
  if (length(missing) > 0) {
    if (length(missing) > 6) missing <- c(missing[1:5], sprintf("+ %i more", length(missing) - 5))
    stop(simpleError(sprintf("OTUs missing from tree: %s", paste(collapse=",", missing))))
  }


  return (tree)
}




PB.TaxaLevelNames <- function (n) {
  if (n <= 8)
    return(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:n])
  sprintf("TaxLvl.%i", 1:n)
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