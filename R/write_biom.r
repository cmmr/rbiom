
#' Write counts, metadata, taxonomy, and phylogeny to a biom file.
#' 
#' @param biom  The BIOM object to save to the file. If another class of
#'        object is given, it will be coerced to matrix and output in
#'        tabular format, provided it is numeric with rownames and colnames.
#' 
#' @param file  Path to the output file. If the file name ends in \code{".gz"} 
#'        or \code{".bz2"}, the file contents will be compressed accordingly.
#' 
#' @param format  Options are \code{"tab"}, \code{"json"}, and \code{"hdf5"}, 
#'        corresponding to classic tabular format, biom format version 1.0 and 
#'        biom version 2.1, respectively. 
#'        See \url{http://biom-format.org/documentation/} for details.
#'         NOTE: to write HDF5 formatted BIOM files, the BioConductor R package 
#'        \code{rhdf5} must be installed.
#' 
#' @return The normalized filepath that was written to (invisibly).
#' @export


write_biom <- function (biom, file, format="json") {
  
  stopifnot(is_scalar_character(file) && !is_na(file))
  stopifnot(is_string(format, c("tab", "json", "hdf5")))
  
  
  #________________________________________________________
  # Refuse to overwrite existing files.
  #________________________________________________________
  
  file <- normalizePath(file, winslash = "/", mustWork = FALSE)
  
  if (!dir.exists(dirname(file)))
    dir.create(dirname(file), recursive = TRUE)
  
  if (file.exists(file))
      stop("Output file already exists: ", file)
  
  
  
  #________________________________________________________
  # Try to convert non-BIOM objects to a matrix for tsv output.
  #________________________________________________________
  if (!is(biom, "BIOM")) {
    
    mtx <- try(as.matrix(biom), silent = TRUE)
    if (is.matrix(mtx) && is.numeric(mtx))
      if (!is_null(rownames(mtx)) && !is_null(colnames(mtx)))
        return (write_biom.tsv(mtx, file))
    
    stop("Invalid BIOM object.")
  }
  
  
  #________________________________________________________
  # Default values for required fields
  #________________________________________________________
  
  if (is_null(biom$info$type))     biom$info$type <- "OTU table"
  if (is.na(biom$info$type))       biom$info$type <- "OTU table"
  if (length(biom$info$type) != 1) biom$info$id   <- "OTU table"
  if (nchar(biom$info$type) == 0)  biom$info$type <- "OTU table"
  if (is_null(biom$info$id))       biom$info$id   <- "NA"
  if (is.na(biom$info$id))         biom$info$id   <- "NA"
  if (length(biom$info$id) != 1)   biom$info$id   <- "NA"
  if (nchar(biom$info$id) == 0)    biom$info$id   <- "NA"
  
  switch (format, 
    "hdf5" = write_biom.2.1(biom, file),
    "json" = write_biom.1.0(biom, file),
    "tab"  = write_biom.tsv(biom, file),
    stop("unknown output format: ", format)
  )
}



#________________________________________________________
# BIOM Classic - Tabular
#________________________________________________________

write_biom.tsv <- function (biom, file) {
  
  # biom could be either a BIOM or matrix object
  counts   <- if (is(biom, "BIOM")) biom$counts   else biom
  taxonomy <- if (is(biom, "BIOM")) biom$taxonomy else NULL
  
  mtx <- {
    rbind(
      matrix(
        data = c("#OTU ID", colnames(counts)),
        nrow = 1
      ),
      cbind(
        matrix(
          data = rownames(counts),
          ncol = 1
        ),
        matrix(
          data = as.character(as.matrix(counts)), 
          nrow = nrow(counts)
        )
      )
    )
  }
  
  if (!is_null(taxonomy)) {
    if (ncol(taxonomy) > 0) {
      mtx <- {
        cbind(
          mtx,
          matrix(
            data = c("Taxonomy", apply(taxonomy, 1L, paste, collapse="; ")),
            ncol = 1
          )
        )
      }
    }
  }
  
  if        (grepl("\\.gz$",  tolower(file))) { con <- gzfile(file, "w")
  } else if (grepl("\\.bz2$", tolower(file))) { con <- bzfile(file, "w")
  } else                                      { con <- base::file(file, "w") }
  
  write.table(mtx, con, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  close(con)
  
  return (invisible(file))
}



#________________________________________________________
# BIOM v1.0 - JSON
#________________________________________________________

write_biom.1.0 <- function (biom, file) {
  
  for (i in names(biom$metadata))
    if (is(biom$metadata[[i]], "factor"))
      biom$metadata[[i]] <- as.character(biom$metadata[[i]])
  
  
  json <- jsonlite::toJSON(
    auto_unbox = TRUE, 
    x          = list(
    
    
    # Attributes
    #________________________________________________________
    id                  = biom$info$id   %||% "",
    type                = biom$info$type %||% "",
    format              = "1.0.0", 
    format_url          = "http://biom-format.org",
    generated_by        = paste("rbiom", utils::packageDescription('rbiom')$Version),
    date                = strftime(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC"),
    matrix_type         = "sparse",
    matrix_element_type = ifelse(sum(biom$counts$v %% 1) == 0, "int", "float"),
    shape               = c(biom$counts$nrow, biom$counts$ncol),
    comment             = biom$info$comment %||% "",
    
    
    # Phylogenetic Tree
    #________________________________________________________
    phylogeny = ifelse(is_null(biom$phylogeny), "", write_tree(biom$phylogeny)),
    
    
    # Taxonomy
    #________________________________________________________
    rows = lapply(1:biom$counts$nrow, function (i) {
      TaxaID   <- biom$counts$dimnames[[1]][i]
      Metadata <- list(taxonomy=unname(biom$taxonomy[i,]))
      
      if (is(biom[['sequences']], "character"))
        Metadata[['sequence']] <- biom$sequences[[TaxaID]]
      
      list(id=TaxaID, metadata=Metadata)
    }),
    
    
    # Sample IDs and Metadata
    #________________________________________________________
    columns = lapply(1:biom$counts$ncol, function (j) list(
      id = biom$counts$dimnames[[2]][j],
      metadata = {
        md_fields <- names(biom$metadata)
        if (length(md_fields) > 0) {
          as.list(biom$metadata[j,,drop=FALSE])
        } else {
          NULL
        }
      }
    )),
    
    
    # Read counts
    #________________________________________________________
    data = lapply(seq_along(biom$counts$v), function (k) 
      c(biom$counts$i[k] - 1, biom$counts$j[k] - 1, biom$counts$v[k])
    )
  ))
  
  
  if        (grepl("\\.gz$",  tolower(file))) { con <- gzfile(file, "w")
  } else if (grepl("\\.bz2$", tolower(file))) { con <- bzfile(file, "w")
  } else                                      { con <- base::file(file, "w") }
  
  res <- try(writeChar(json, con, eos=NULL), silent = TRUE)
  if (is(res, "try-error"))
    stop(sprintf("Can't save to '%s': %s", file, res))
  
  close(con)
  
  return (invisible(file))
}




#________________________________________________________
# BIOM v2.1 - HDF5
#________________________________________________________

write_biom.2.1 <- function (biom, file) {
  
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop(paste0(
      "\n",
      "Error: rbiom requires the R package 'rhdf5' to be installed\n",
      "in order to read and write HDF5 formatted BIOM files.\n\n",
      "Please run the following commands to install 'rhdf5':\n",
      "   install.packages('BiocManager')\n",
      "   BiocManager::install('rhdf5')\n\n" ))
  }
  
  res <- try(rhdf5::h5createFile(file), silent = TRUE)
  if (!identical(res, TRUE))
    stop(sprintf("Can't create file '%s': %s", file, as.character(res)))
  
  invisible(rhdf5::h5createGroup(file, '/observation'))
  invisible(rhdf5::h5createGroup(file, '/observation/matrix'))
  invisible(rhdf5::h5createGroup(file, '/observation/metadata'))
  invisible(rhdf5::h5createGroup(file, '/observation/group-metadata'))
  invisible(rhdf5::h5createGroup(file, '/sample'))
  invisible(rhdf5::h5createGroup(file, '/sample/matrix'))
  invisible(rhdf5::h5createGroup(file, '/sample/metadata'))
  invisible(rhdf5::h5createGroup(file, '/sample/group-metadata'))
  
  h5 <- try(rhdf5::H5Fopen(file), silent = TRUE)
  if (!is(h5, "H5IdComponent"))
    stop(sprintf("Can't open file '%s': %s", file, as.character(h5)))
  
  
  
  # Attributes
  #________________________________________________________
  rhdf5::h5writeAttribute(as.character(biom$info$id      %||% ""),                    h5, 'id')
  rhdf5::h5writeAttribute(as.character(biom$info$type    %||% ""),                    h5, 'type')
  rhdf5::h5writeAttribute(as.character(biom$info$comment %||% ""),                    h5, 'comment')
  rhdf5::h5writeAttribute("http://biom-format.org",                                   h5, 'format-url')
  rhdf5::h5writeAttribute(as.integer(c(2,1,0)),                                       h5, 'format-version', 3)
  rhdf5::h5writeAttribute(paste("rbiom", utils::packageDescription('rbiom')$Version), h5, 'generated-by')
  rhdf5::h5writeAttribute(strftime(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC"),       h5, 'creation-date')
  rhdf5::h5writeAttribute(as.integer(c(biom$counts$nrow, biom$counts$ncol)),          h5, 'shape', 2)
  rhdf5::h5writeAttribute(as.integer(length(biom$counts$v)),                          h5, 'nnz')
  
  
  
  
  # Read counts by taxa (rows)
  #________________________________________________________
  x <- matrix(c(biom$counts$i - 1, biom$counts$j - 1, biom$counts$v), byrow=FALSE, ncol=3)
  x <- x[order(x[,1]),,drop=FALSE]
  rhdf5::h5writeDataset(as.character(biom$counts$dimnames[[1]]),                                h5, 'observation/ids')
  rhdf5::h5writeDataset(as.numeric(x[,3]),                                                      h5, 'observation/matrix/data')
  rhdf5::h5writeDataset(as.integer(x[,2]),                                                      h5, 'observation/matrix/indices')
  rhdf5::h5writeDataset(as.integer(cumsum(unname(table(factor(x[,1]+1, 0:biom$counts$nrow))))), h5, 'observation/matrix/indptr')
  
  
  
  # Read counts by sample (columns)
  #________________________________________________________
  x <- x[order(x[,2]),,drop=FALSE]
  rhdf5::h5writeDataset(as.character(biom$counts$dimnames[[2]]),                                h5, 'sample/ids')
  rhdf5::h5writeDataset(as.numeric(x[,3]),                                                      h5, 'sample/matrix/data')
  rhdf5::h5writeDataset(as.integer(x[,1]),                                                      h5, 'sample/matrix/indices')
  rhdf5::h5writeDataset(as.integer(cumsum(unname(table(factor(x[,2]+1, 0:biom$counts$ncol))))), h5, 'sample/matrix/indptr')
  
  
  
  # Sample Metadata
  #________________________________________________________
  if (is(biom[['metadata']], "data.frame")) {
    plyr::l_ply(setdiff(names(biom$metadata), 'SampleID'), function (field) {
      
      h5path <- sprintf("/sample/metadata/%s", field)
      values <- biom$metadata[biom$counts$dimnames[[2]], field]
      
      if (is.numeric(values)) {
        if (all(values %% 1 == 0, na.rm=TRUE))
          values <- as.integer(values)
        
      } else if (!is.logical(values)) {
        values <- as.character(values)
      }
      
      rhdf5::h5writeDataset(values, h5, h5path)
      
    })
  }
  
  
  
  # Taxonomy
  #________________________________________________________
  if (is(biom[['taxonomy']], "matrix")) {
    h5path      <- 'observation/metadata/taxonomy'
    x           <- t(biom$taxonomy[biom$counts$dimnames[[1]],,drop=FALSE])
    dimnames(x) <- list(NULL, NULL)
    rhdf5::h5writeDataset(x, h5, h5path)
  }
  
  
  
  # Sequences
  #________________________________________________________
  if (is(biom[['sequences']], "character")) {
    h5path <- 'observation/metadata/sequences'
    x      <- unname(biom$sequences[biom$counts$dimnames[[1]]])
    rhdf5::h5writeDataset(x, h5, h5path)
  }
  
  
  
  # Phylogenetic Tree
  #________________________________________________________
  if (is(biom[['phylogeny']], "phylo")) {
    
    h5path <- '/observation/group-metadata/phylogeny'
    x      <- write_tree(biom[['phylogeny']])
    rhdf5::h5writeDataset(x, h5, h5path)
    
    h5path <- h5&'observation'&'group-metadata'&'phylogeny'
    rhdf5::h5writeAttribute("newick", h5path, 'data_type')
  }
  
  rhdf5::H5Fflush(h5)
  rhdf5::H5Fclose(h5)
  
  return (invisible(file))
}
