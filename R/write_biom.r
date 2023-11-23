
#' Write counts, metadata, taxonomy, and phylogeny to a biom file.
#' 
#' @param biom  The \code{rbiom} object to save to the file.
#' 
#' @param file  Path to the output file. If the file name ends in \code{".gz"} 
#'        or \code{".bz2"}, the file contents will be compressed accordingly.
#' 
#' @param format  Options are \code{"tab"}, \code{"json"}, and \code{"hdf5"}, 
#'        corresponding to classic tabular format, BIOM format version 1.0 and 
#'        biom version 2.1, respectively. 
#'        See \url{http://biom-format.org/documentation/} for details.
#'         NOTE: to write HDF5 formatted BIOM files, the BioConductor R package 
#'        \code{rhdf5} must be installed.
#' 
#' @return The normalized filepath that was written to (invisibly).
#' @export


write_biom <- function (biom, file, format="json") {
  
  validate_biom(clone = TRUE)
  
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
  # BIOM fields pertaining to BIOM generator.
  #________________________________________________________
  
  biom[['info']][['type']]         <- if.null("OTU table")
  biom[['info']][['format_url']]   <- "http://biom-format.org"
  biom[['info']][['generated_by']] <- if.null(paste("rbiom", packageVersion("rbiom")))
  biom[['info']][['date']]         <- if.null(strftime(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC"))
  biom[['info']][['shape']]        <- dim(biom[['counts']])
  biom[['info']][['matrix_type']]  <- "sparse"
  biom[['info']][['nnz']]          <- length(biom[['counts']][['v']])
  biom[['info']][['matrix_element_type']] <- ifelse(
    test = all(biom[['counts']][['v']] %% 1 == 0), 
    yes  = "int", 
    no   = "float" )
  
  
  
  #________________________________________________________
  # Convert factors to character, tibbles to base types.
  #________________________________________________________
  metadata <- metadata %>%
    mutate(across(where(as.factor), as.character)) %>%
    tibble::column_to_rownames('.sample') %>%
    as.data.frame()
  
  taxonomy <- taxonomy %>%
    tibble::column_to_rownames('.otu') %>%
    as.matrix()
  
  
    
  switch (format, 
    "hdf5" = write_biom_hdf5(biom, file),
    "json" = write_biom_json(biom, file),
    "tab"  = write_biom_tsv(biom, file),
    stop("unknown output format: ", format)
  )
}



#________________________________________________________
# BIOM Classic - Tabular
#________________________________________________________

write_biom_tsv <- function (biom, file) {
  
  counts   <- biom[['counts']]
  taxonomy <- biom[['taxonomy']]
  
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

write_biom_json <- function (biom, file) {
  
  
  json <- with(biom, {
    
    
    
    json <- jsonlite::toJSON(
      auto_unbox = TRUE, 
      x          = list(
        
        
        # Attributes
        #________________________________________________________
        format              = "1.0.0", 
        id                  = info[['id']],
        type                = info[['type']],
        format_url          = info[['format_url']],
        generated_by        = info[['generated_by']],
        date                = info[['date']],
        matrix_type         = info[['matrix_type']],
        matrix_element_type = info[['matrix_element_type']],
        shape               = info[['shape']],
        comment             = info[['comment']],
        
        
        # Phylogenetic Tree
        #________________________________________________________
        phylogeny = ifelse(is_null(phylogeny), "", write_tree(phylogeny)),
        
        
        # Taxonomy
        #________________________________________________________
        rows = lapply(1:counts$nrow, function (i) {
          TaxaID   <- counts$dimnames[[1]][i]
          Metadata <- list(taxonomy=unname(taxonomy[i,]))
          
          if (is(biom[['sequences']], "character"))
            Metadata[['sequence']] <- sequences[[TaxaID]]
          
          list(id=TaxaID, metadata=Metadata)
        }),
        
        
        # Sample IDs and Metadata
        #________________________________________________________
        columns = lapply(1:counts$ncol, function (j) list(
          id = counts$dimnames[[2]][j],
          metadata = {
            md_fields <- names(metadata)
            if (length(md_fields) > 0) {
              as.list(metadata[j,,drop=FALSE])
            } else {
              NULL
            }
          }
        )),
        
        
        # Read counts
        #________________________________________________________
        data = lapply(seq_along(counts$v), function (k) 
          c(counts$i[k] - 1, counts$j[k] - 1, counts$v[k])
        )
      ))
    
    return (json)
  })
  
  
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

write_biom_hdf5 <- function (biom, file) {
  
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
  if (!eq(res, TRUE))
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
  
  
  
  with(biom, {
    
    
    # Attributes
    #________________________________________________________
    rhdf5::h5writeAttribute(as.character(info[['id']]),           h5, 'id')
    rhdf5::h5writeAttribute(as.character(info[['type']]),         h5, 'type')
    rhdf5::h5writeAttribute(as.character(info[['comment']]),      h5, 'comment')
    rhdf5::h5writeAttribute(as.character(info[['format_url']]),   h5, 'format-url')
    rhdf5::h5writeAttribute(as.integer(c(2,1,0)),                 h5, 'format-version', 3)
    rhdf5::h5writeAttribute(as.character(info[['generated_by']]), h5, 'generated-by')
    rhdf5::h5writeAttribute(as.character(info[['date']]),         h5, 'creation-date')
    rhdf5::h5writeAttribute(as.integer(info[['shape']]),          h5, 'shape', 2)
    rhdf5::h5writeAttribute(as.integer(info[['nnz']]),            h5, 'nnz')
    
    
    
    
    # Read counts by taxa (rows)
    #________________________________________________________
    x <- matrix(c(counts$i - 1, counts$j - 1, counts$v), byrow=FALSE, ncol=3)
    x <- x[order(x[,1]),,drop=FALSE]
    rhdf5::h5writeDataset(as.character(counts$dimnames[[1]]),                                h5, 'observation/ids')
    rhdf5::h5writeDataset(as.numeric(x[,3]),                                                 h5, 'observation/matrix/data')
    rhdf5::h5writeDataset(as.integer(x[,2]),                                                 h5, 'observation/matrix/indices')
    rhdf5::h5writeDataset(as.integer(cumsum(unname(table(factor(x[,1]+1, 0:counts$nrow))))), h5, 'observation/matrix/indptr')
    
    
    
    # Read counts by sample (columns)
    #________________________________________________________
    x <- x[order(x[,2]),,drop=FALSE]
    rhdf5::h5writeDataset(as.character(counts$dimnames[[2]]),                                h5, 'sample/ids')
    rhdf5::h5writeDataset(as.numeric(x[,3]),                                                 h5, 'sample/matrix/data')
    rhdf5::h5writeDataset(as.integer(x[,1]),                                                 h5, 'sample/matrix/indices')
    rhdf5::h5writeDataset(as.integer(cumsum(unname(table(factor(x[,2]+1, 0:counts$ncol))))), h5, 'sample/matrix/indptr')
    
    
    
    # Sample Metadata
    #________________________________________________________
    if (is(biom[['metadata']], "data.frame")) {
      plyr::l_ply(setdiff(names(metadata), 'SampleID'), function (field) {
        
        h5path <- sprintf("/sample/metadata/%s", field)
        values <- metadata[counts$dimnames[[2]], field]
        
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
      x           <- t(taxonomy[counts$dimnames[[1]],,drop=FALSE])
      dimnames(x) <- list(NULL, NULL)
      rhdf5::h5writeDataset(x, h5, h5path)
    }
    
    
    
    # Sequences
    #________________________________________________________
    if (is(biom[['sequences']], "character")) {
      h5path <- 'observation/metadata/sequences'
      x      <- unname(sequences[counts$dimnames[[1]]])
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
    
  })
  
  
  rhdf5::H5Fflush(h5)
  rhdf5::H5Fclose(h5)
  
  return (invisible(file))
}
