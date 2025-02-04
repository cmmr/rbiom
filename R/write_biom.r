
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
#'        files, the BioConductor R package \code{rhdf5} must be installed.
#'        Default: `"json"`
#' 
#' @param depth,n   Passed on to [rarefy_cols()]. For `write_xlsx()` only, 
#'        `depth=0` disables rarefaction. Default: `depth=0.1, n=NULL`
#'                   
#' @param seed   Random seed to use in rarefying. See [rarefy_cols()] function
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
#'       attr(hmp10, "Weighted UniFrac")   <- bdiv_distmat(hmp10, 'unifrac')
#'       attr(hmp10, "Unweighted Jaccard") <- bdiv_distmat(hmp10, 'jaccard', weighted=FALSE)
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
      matrix_element_type = ifelse(all(biom$counts[['v']] %% 1 == 0), "int", "float"),
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
        biom$counts$i - 1, 
        biom$counts$j - 1, 
        biom$counts$v )
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
  
  require_package('rhdf5', 'to write HDF5 formatted BIOM files')
  
  res <- try(h5createFile(file), silent = TRUE)
  if (!isTRUE(res))
    cli_abort("Can't create file {.file {file}}: {res}")
  
  invisible(h5createGroup(file, 'observation'))
  invisible(h5createGroup(file, 'observation/matrix'))
  invisible(h5createGroup(file, 'observation/metadata'))
  invisible(h5createGroup(file, 'observation/group-metadata'))
  invisible(h5createGroup(file, 'sample'))
  invisible(h5createGroup(file, 'sample/matrix'))
  invisible(h5createGroup(file, 'sample/metadata'))
  invisible(h5createGroup(file, 'sample/group-metadata'))
  
  h5 <- try(H5Fopen(file, 'H5F_ACC_RDWR', native = TRUE), silent = TRUE)
  if (!inherits(h5, "H5IdComponent"))
    cli_abort("Can't open HDF5 file {.file {file}}: {h5}")
  
  
  
  # Attributes
  #________________________________________________________
  h5writeAttribute(biom$id,                    h5, 'id')
  h5writeAttribute("OTU table",                h5, 'type')
  h5writeAttribute(biom$comment,               h5, 'comment')
  h5writeAttribute("http://biom-format.org",   h5, 'format-url')
  h5writeAttribute(as.integer(c(2,1,0)),       h5, 'format-version') #, 3)
  h5writeAttribute(biom$date,                  h5, 'creation-date')
  h5writeAttribute(dim(biom$counts),           h5, 'shape') #, 2)
  h5writeAttribute(length(biom$counts[['v']]), h5, 'nnz')
  h5writeAttribute(paste("rbiom", packageVersion("rbiom")), h5, 'generated-by')
  
  
  
  # Read counts by taxa (rows)
  #________________________________________________________
  x <- matrix(c(biom$counts$i - 1, biom$counts$j - 1, biom$counts$v), byrow=FALSE, ncol=3)
  
  x      <- x[order(x[,1]),,drop=FALSE]
  indptr <- cumsum(unname(table(factor(x[,1]+1, 0:nrow(biom$counts)))))
  
  h5writeDataset(rownames(biom$counts), h5, 'observation/ids')
  h5writeDataset(as.numeric(x[,3]),     h5, 'observation/matrix/data')
  h5writeDataset(as.integer(x[,2]),     h5, 'observation/matrix/indices')
  h5writeDataset(as.integer(indptr),    h5, 'observation/matrix/indptr')
  
  
  
  # Read counts by sample (columns)
  #________________________________________________________
  x      <- x[order(x[,2]),,drop=FALSE]
  indptr <- cumsum(unname(table(factor(x[,2]+1, 0:ncol(biom$counts)))))
  
  h5writeDataset(colnames(biom$counts), h5, 'sample/ids')
  h5writeDataset(as.numeric(x[,3]),     h5, 'sample/matrix/data')
  h5writeDataset(as.integer(x[,1]),     h5, 'sample/matrix/indices')
  h5writeDataset(as.integer(indptr),    h5, 'sample/matrix/indptr')
  
  
  
  # Sample Metadata
  #________________________________________________________
  if (ncol(biom$metadata) > 1) {
    plyr::l_ply(names(biom$metadata)[-1], function (field) {
      
      if (grepl('/', field, fixed = TRUE))
        cli_abort(c(
          'i' = "HDF5 BIOM spec prohibits slashes ({.var /}) in metadata field names.", 
          'i' = "Either change the field name, or save to JSON format instead.", 
          'x' = "Metadata field {.val {field}} not encodable." ))
      
      h5path <- sprintf("sample/metadata/%s", field)
      values <- biom$metadata[[field]]
      
      if (is.numeric(values)) {
        if (all(values %% 1 == 0, na.rm=TRUE))
          values <- as.integer(values)
        
      } else {
        values <- as.character(values)
      }
      
      h5writeDataset(values, h5, h5path)
      
    })
  }
  
  
  
  # Taxonomy
  #________________________________________________________
  if (ncol(biom$taxonomy) > 1) {
    h5path      <- 'observation/metadata/taxonomy'
    x           <- t(as.matrix(biom$taxonomy[,-1]))
    dimnames(x) <- list(NULL, NULL)
    h5writeDataset(x, h5, h5path)
  }
  
  
  
  # Sequences
  #________________________________________________________
  if (!is.null(biom$sequences)) {
    h5path <- 'observation/metadata/sequences'
    x      <- unname(biom$sequences)
    h5writeDataset(x, h5, h5path)
  }
  
  
  
  # Phylogenetic Tree
  #________________________________________________________
  if (!is.null(biom$tree)) {
    
    h5path <- 'observation/group-metadata/phylogeny'
    x      <- write_tree(biom$tree)
    h5writeDataset(x, h5, h5path)
    
    h5d <- H5Dopen(h5, h5path)
    h5writeAttribute("newick", h5d, 'data_type')
    H5Dclose(h5d)
  }
  
  
  H5Fclose(h5)
  
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
      dplyr::mutate('total' = as.vector(row_sums(biom$counts)), .after = 1) %>% 
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
        'Size'     = row_sums(biom$counts) %>% as.vector(),
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



