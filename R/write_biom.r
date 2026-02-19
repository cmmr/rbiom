
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
#'        `write_biom(format='json')`, and returns a string of the output which 
#'        would have been written. For `write_biom(format='tab')`, `file=NULL`
#'        returns the tibble that would have been written.
#' 
#' @param format  Options are \code{"tab"}, \code{"json"}, and \code{"hdf5"}, 
#'        corresponding to classic tabular format, BIOM format version 1.0 and 
#'        biom version 2.1, respectively. NOTE: to write HDF5 formatted BIOM 
#'        files, the \code{h5lite} R package must be installed. Default: `"json"`
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
#'     if (FALSE) {
#'     
#'       hmp10        <- hmp50$clone()
#'       hmp10$counts <- hmp10$counts[,1:10] %>% rarefy_cols()
#'       
#'       attr(hmp10, "Weighted UniFrac") <- bdiv_distmat(hmp10, 'wunifrac')
#'       attr(hmp10, "Jaccard")          <- bdiv_distmat(hmp10, 'jaccard')
#'       
#'       outfile <- write_xlsx(hmp10, tempfile(fileext = ".xlsx"))
#'     }
#' 
write_biom <- function (biom, file, format = "json") {
  
  biom <- as_rbiom(biom)
  
  format <- match.arg(tolower(format), c("tab", "json", "hdf5"))
  
  #________________________________________________________
  # Allow NULL for json and tab
  #________________________________________________________
  if (format == "hdf5" || !is.null(file)) {
    
    stopifnot(is_scalar_character(file) && !is_na(file))
  
  
    #________________________________________________________
    # Refuse to overwrite existing files.
    #________________________________________________________
    
    file <- normalizePath(file, winslash = "/", mustWork = FALSE)
    
    if (!dir.exists(dirname(file)))
      dir.create(dirname(file), recursive = TRUE)
    
    if (file.exists(file))
        stop("Output file already exists: ", file)
  }
  
  
  switch (format, 
    "hdf5" = write_biom_hdf5(biom, file),
    "json" = write_biom_json(biom, file),
    "tab"  = write_biom_tsv(biom, file) )
}



#________________________________________________________
# BIOM Classic - Tabular
#________________________________________________________

write_biom_tsv <- function (biom, file) {
  
  
  # Sanity checks.
  #________________________________________________________
  if (any(biom$samples %in% c("#OTU ID", "Taxonomy")))
    cli_abort(c(x = "Sample names can't be '#OTU ID' or 'Taxonomy' when saving to tab format."))
  
  
  
  # Convert slam matrix to tibble.
  #________________________________________________________
  tbl <- as.matrix(biom$counts) %>%
    tibble::as_tibble(rownames = "#OTU ID")
  
  
  
  # Add Taxonomy column when appropriate.
  #________________________________________________________
  if (ncol(biom$taxonomy) > 1)
    tbl[['Taxonomy']] <- apply(biom$taxonomy[,-1], 1L, paste, collapse="; ")
  
  
  if (is.null(file)) return (tbl)
  
  
  # Compress when writing according to file extension.
  #________________________________________________________
  if        (grepl("\\.gz$",  tolower(file))) { con <- gzfile(file, "w")
  } else if (grepl("\\.bz2$", tolower(file))) { con <- bzfile(file, "w")
  } else                                      { con <- base::file(file, "w") }
  
  write.table(tbl, con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  close(con)
  
  
  return (invisible(file))
}



#________________________________________________________
# BIOM v1.0 - JSON
#________________________________________________________

write_biom_json <- function (biom, file) {
  
  metadata <- biom$metadata
  dgT      <- as(biom$counts, 'TsparseMatrix')
  
  json <- toJSON(
    auto_unbox = TRUE, 
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
      rows = apply(biom$taxonomy, 1L, function (x) {
        
        otu <- x[[1]]
        dat <- list(taxonomy=unname(x[-1]))
        
        if (!is.null(biom$sequences))
          dat[['sequence']] <- biom$sequences[[otu]]
        
        list(id=otu, metadata=dat)
      }) |> setNames(NULL),
      
      
      # Sample IDs and Metadata
      #________________________________________________________
      columns = lapply(seq_len(nrow(metadata)), function (i) {
        row <- as.list(metadata[i,])
        list(
          id       = row[[1]], 
          metadata = if (length(row) > 1) row[-1] else NULL )
      }) |> setNames(NULL),
      
      
      # Read counts
      #________________________________________________________
      data = mapply(
          FUN      = c,
          SIMPLIFY = FALSE,
          dgT@i,
          dgT@j,
          dgT@x
        ) |> setNames(NULL)
      
    ))
  
  
  if (is.null(file)) return (as.character(json))
  
  if        (grepl("\\.gz$",  tolower(file))) { con <- gzfile(file, "w")
  } else if (grepl("\\.bz2$", tolower(file))) { con <- bzfile(file, "w")
  } else                                      { con <- base::file(file, "w") }
  
  res <- try(writeChar(json, con, eos=NULL), silent = TRUE)
  if (inherits(res, "try-error"))
    stop(sprintf("Can't save to '%s': %s", file, res))
  
  close(con)
  
  return (invisible(file))
}




#________________________________________________________
# BIOM v2.1 - HDF5
#________________________________________________________

write_biom_hdf5 <- function (biom, file) {
  
  require_package('h5lite', 'to write HDF5 formatted BIOM files')
  
  h5 <- h5lite::h5_open(file = file)


  # Required top-level attributes:
  #________________________________________________________
  h5$write(attr = "creation-date",  as = "ascii", data = I(biom$date))
  h5$write(attr = "format-url",     as = "ascii", data = I("http://biom-format.org"))
  h5$write(attr = "format-version", as = "int64", data = c(2,1,0))
  h5$write(attr = "generated-by",   as = "ascii", data = I(paste("rbiom", packageVersion("rbiom"))))
  h5$write(attr = "id",             as = "ascii", data = I(biom$id))
  h5$write(attr = "nnz",            as = "int64", data = I(length(biom$counts@x)))
  h5$write(attr = "shape",          as = "int64", data = dim(biom$counts))
  h5$write(attr = "type",           as = "ascii", data = I("OTU table"))


  # Optional top-level attributes:
  #________________________________________________________
  if (nzchar(biom$comment))
    h5$write(attr = "comment", as = "ascii", data = I(biom$comment))


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
  
  h5$write(name = "observation/ids",            as = "ascii",   data = rownames(biom$counts))
  h5$write(name = "observation/matrix/data",    as = "float64", data = x[,3])
  h5$write(name = "observation/matrix/indices", as = "int32",   data = x[,2])
  h5$write(name = "observation/matrix/indptr",  as = "int32",   data = indptr)
  
  
  # Read counts by sample (columns)
  #________________________________________________________
  x      <- x[order(x[,2]),,drop=FALSE]
  indptr <- cumsum(unname(table(factor(x[,2]+1, 0:ncol(biom$counts)))))
  
  h5$write(name = "sample/ids",            as = "ascii",   data = colnames(biom$counts))
  h5$write(name = "sample/matrix/data",    as = "float64", data = x[,3])
  h5$write(name = "sample/matrix/indices", as = "int32",   data = x[,1])
  h5$write(name = "sample/matrix/indptr",  as = "int32",   data = indptr)
  
  
  # Sample Metadata
  #________________________________________________________
  if (ncol(biom$metadata) > 1) {
    plyr::l_ply(names(biom$metadata)[-1], function (field) {
      
      if (grepl('/', field, fixed = TRUE))
        cli_abort(c(
          'i' = "HDF5 BIOM spec prohibits slashes ({.var /}) in metadata field names.", 
          'i' = "Either change the field name, or save to JSON format instead.", 
          'x' = "Metadata field {.val {field}} not encodable." ))
      
      name <- sprintf("sample/metadata/%s", field)
      data <- biom$metadata[[field]]
      
      if (is.factor(data))
        data <- as.character(data)
      
      h5$write(name = name, data = data)
    })
  }
  
  
  # Taxonomy
  #________________________________________________________
  if (ncol(biom$taxonomy) > 1) {
    name <- 'observation/metadata/taxonomy'
    data <- as.matrix(biom$taxonomy[,-1])
    h5$write(name = name, as = "ascii", data = data)
  }
  
  
  # Sequences
  #________________________________________________________
  if (!is.null(biom$sequences)) {
    name <- 'observation/metadata/sequences'
    h5$write(name = name, data = biom$sequences)
  }
  
  
  # Phylogenetic Tree
  #________________________________________________________
  if (!is.null(biom$tree)) {
    name <- 'observation/group-metadata/phylogeny'
    h5$write(name = name, data = I(write_tree(biom$tree)))
    h5$write(name = name, attr = "data_type", data = "newick")
  }
  
  return (invisible(file))
}




#' @rdname write_biom
#' @export
write_metadata <- function (biom, file, quote = FALSE, sep = "\t", ...) {
  
  write_wrapper(file, function (con) {
    as_rbiom(biom)$metadata %>%
      write.table(file = con, quote = quote, sep = sep, row.names = FALSE, ...)
  })
}


#' @rdname write_biom
#' @export
write_counts <- function (biom, file, quote = FALSE, sep = "\t", ...) {
  
  write_wrapper(file, function (con) {
    as_rbiom(biom)$counts %>% 
      as.matrix() %>% 
      write.table(file = con, sep = sep, quote = quote, ...)
  })
}


#' @rdname write_biom
#' @export
write_taxonomy <- function (biom, file, quote = FALSE, sep = "\t", ...) {
  
  write_wrapper(file, function (con) {
    as_rbiom(biom)$taxonomy %>%
      write.table(file = con, quote = quote, sep = sep, row.names = FALSE, ...)
  })
}




#' @rdname write_biom
#' @export
write_fasta <- function (biom, file = NULL) {
  
  if (is.character(biom) && !is.null(names(biom))) {
    seqs <- biom
  } else {
    seqs <- as_rbiom(biom)$sequences
  }
  
  if (is.null(seqs) || length(seqs) == 0)
    cli_abort(c('x' = "rbiom object does not have sequences."))
  
  
  write_wrapper(file, function (con) {
    cat(file=con, sep="\n", sprintf(">%s\n%s", names(seqs), seqs))
  })
}




#' @rdname write_biom
#' @export
write_tree <- function (biom, file = NULL) {
  
  tree <- if (inherits(biom, "phylo")) biom else as_rbiom(biom)$tree
  
  if (is.null(tree))
    cli_abort(c('x' = "rbiom object does not have a tree."))
  
  
  rootNode <- setdiff(tree$edge[,1], tree$edge[,2])
  parentAt <- stats::aggregate(1:nrow(tree$edge), by=list(tree$edge[,1]), c, simplify=FALSE)
  parentAt <- setNames(lapply(parentAt[,2], unlist), parentAt[,1])
  
  fx <- function (root=NULL) ({
    
    nodes <- parentAt[[as.character(root)]]
    
    if (length(nodes) == 0) {
      
      nodeLabel <- tree$tip.label[root]
      
      if (any(grepl(" ", nodeLabel, fixed=TRUE))) {
        if (any(grepl("_", nodeLabel, fixed=TRUE))) {
          nodeLabel <- paste0("'", nodeLabel, "'")
        } else {
          nodeLabel <- gsub(" ", "_", nodeLabel)
        }
      }
      return (nodeLabel)
    }
    
    children <- tree$edge[nodes, 2]
    children <- sapply(children, fx)
    
    if (!is_null(tree$edge.length))
      children <- paste(sep=":", children, tree$edge.length[nodes])
    
    sprintf("(%s)", paste(collapse=",", children))
  })
  
  newick <- paste0(fx(rootNode), ";")
  
  
  if (is_null(file)) return (newick)
  
  
  write_wrapper(file, function (con) {
    writeLines(text=newick, con=con, sep="")
  })
  
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
#' @return The normalized filepath that was written to (invisibly).
#' 
write_wrapper <- function (outfile, callback) {
  
  stopifnot(is_scalar_character(outfile) && !is_na(outfile))
  stopifnot(isTRUE(nchar(outfile) > 0))
  
  
  outfile <- normalizePath(outfile, winslash = "/", mustWork = FALSE)
  
  tryCatch(
    error = function (e) stop("Can't write to '", outfile, "'.\n", e),
    expr  = local({
      
      if (!dir.exists(dirname(outfile)))
        dir.create(dirname(outfile), recursive = TRUE)
      
      con <- if (endsWith(outfile, ".gz"))  { gzfile(outfile, "w")
      } else if (endsWith(outfile, ".bz2")) { bzfile(outfile, "w")
      } else                                { file(outfile, "w") }
      
      callback(con)
      
      close(con)
    }))
  
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
#' 
#' @name export
#' @inherit documentation_default
#' 
#' @param dir  Where to save the files. If the directory doesn't exist, it will 
#'        be created. Default: `tempfile()`
#' 
#' @param prefix  A string to prepend to each file name. Default: `'biom_'`
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

write_mothur <- function (biom, dir = tempfile(), prefix = 'biom_') {
  
  biom <- as_rbiom(biom)
  
  dir <- normalizePath(dir, winslash = '/', mustWork = FALSE)
  if (!dir.exists(dir)) dir.create(path = dir, recursive = TRUE)
  
  validate_string('prefix')
  
  
  # Counts
  write_wrapper(file.path(dir, paste0(prefix, 'counts.tsv')), function (con) {
    as_rbiom(biom)$counts %>% 
      as.matrix() %>%
      tibble::as_tibble(rownames = "Representative_Sequence") %>%
      dplyr::mutate('total' = as.vector(rowSums(biom$counts)), .after = 1) %>% 
      write.table(file=con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  })
  
  
  # Tree
  if (!is.null(biom$tree))
    write_tree(biom, file.path(dir, paste0(prefix, 'tree.nwk')))
  
  
  # Sequences
  if (!is.null(biom$sequences))
    write_fasta(biom, file.path(dir, paste0(prefix, 'seqs.fna')))
  
  
  # Taxonomy
  if (ncol(biom$taxonomy) > 1)
    write_wrapper(file.path(dir, paste0(prefix, 'taxonomy.tsv')), function (con) {
      tibble::tibble(
        'OTU'      = seq_len(biom$n_otus),
        'Size'     = rowSums(biom$counts) %>% as.vector(),
        'Taxonomy' = apply(as.matrix(biom$taxonomy[,-1]), 1L, paste, collapse=';') %>% gsub(' ', '_', ., fixed = TRUE) ) %>%
        write.table(file=con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    })
  
  
  # Metadata
  if (ncol(biom$metadata) > 1)
    write_wrapper(file.path(dir, paste0(prefix, 'metadata.tsv')), function (con) {
      
      tbl <- as_rbiom(biom)$metadata
      
      colnames(tbl) %<>% gsub('[\\s\\t\\n]+', '_', ., perl = TRUE)
      colnames(tbl)[[1]] <- 'group'
      tbl[[1]] %<>% as.character()
      
      for (i in seq_len(ncol(tbl)))
        if (is.factor(tbl[[i]])) {
          x <- levels(tbl[[i]])
          x[grep('[\\s\\t\\"]', x, perl = TRUE)] %<>% shQuote(type = 'cmd')
          levels(tbl[[i]]) <- x
        }
      
      write.table(tbl, file = con, quote = FALSE, sep = "\t", row.names = FALSE)
    })
  
  
  return (invisible(dir))
}



#' @rdname export
#' @export

write_qiime2 <- function (biom, dir = tempfile(), prefix = 'biom_') {
  
  biom <- as_rbiom(biom)
  
  dir <- normalizePath(dir, winslash = '/', mustWork = FALSE)
  if (!dir.exists(dir)) dir.create(path = dir, recursive = TRUE)
  
  validate_string('prefix')
  
  
  # Counts
  write_wrapper(file.path(dir, paste0(prefix, 'counts.tsv')), function (con) {
    as_rbiom(biom)$counts %>% 
      as.matrix() %>%
      tibble::as_tibble(rownames = "#OTU ID") %>%
      write.table(file=con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  })
  
  
  # Tree
  if (!is.null(biom$tree))
    write_tree(biom, file.path(dir, paste0(prefix, 'tree.nwk')))
  
  
  # Sequences
  if (!is.null(biom$sequences))
    write_fasta(biom, file.path(dir, paste0(prefix, 'seqs.fna')))
  
  
  # Taxonomy
  if (ncol(biom$taxonomy) > 1)
    write_wrapper(file.path(dir, paste0(prefix, 'taxonomy.tsv')), function (con) {
      tibble::tibble(
        'Feature ID' = biom$otus,
        'Taxon'      = apply(as.matrix(biom$taxonomy[,-1]), 1L, paste, collapse='; ') %>% gsub(' ', '_', ., fixed = TRUE),
        'Confidence' = rep(1, biom$n_otus) ) %>%
        write.table(file=con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    })
  
  
  # Metadata
  if (ncol(biom$metadata) > 1)
    write_wrapper(file.path(dir, paste0(prefix, 'metadata.tsv')), function (con) {
      
      tbl     <- as_rbiom(biom)$metadata
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
        if (is.factor(tbl[[i]])) {
          cat(file = con, sep = '', '\tcategorical')
          x <- levels(tbl[[i]])
          x[grep('[\\s\\t\\n\\"]', x, perl = TRUE)] %<>% shQuote(type = 'cmd')
          x %<>% gsub('\\"', '""', ., fixed = TRUE)
          levels(tbl[[i]]) <- x
        } else {
          cat(file = con, sep = '', '\tnumeric')
        }
      }
      
      cat(file = con, sep = '', '\n')
      
      write.table(tbl, file = con, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    })
  
  
  return (invisible(dir))
}
