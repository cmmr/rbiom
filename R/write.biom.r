#' Write counts, metadata, taxonomy, and phylogeny to a biom file.
#'
#' @param biom  The BIOM object to save to the file.
#' @param file  Path to the output file.
#' @param format  Options are \bold{\dQuote{tab}}, 
#'     \bold{\dQuote{json}}, and \bold{\dQuote{hdfs}}, 
#'     corresponding to classic tabular format, biom format version 1.0 and 
#'     biom version 2.1, respectively. Abbreviations are also accepted. See
#'     \url{http://biom-format.org/documentation/} for details.
#' @return On success, returns \code{NULL} invisibly.
#' @export


write.biom <- function (biom, file, format="hdf5") {
  
  #--------------------------------------------------------------
  # Sanity Checks
  #--------------------------------------------------------------
  
  if (!is(biom, "BIOM"))
      stop(simpleError("Invalid BIOM object."))
  
  
  #--------------------------------------------------------------
  # Select the approriate format engine
  #--------------------------------------------------------------
  
  opts   <- c("tab", "json", "hdf5")
  format <- tolower(head(format, 1))
  format <- c(opts, format)[pmatch(format, opts, nomatch=4)]
  
  switch (format, 
    "hdf5" = write.biom.2.1(biom, file),
    "json" = write.biom.1.0(biom, file),
    "tab"  = write.biom.tsv(biom, file),
    stop(simpleError(sprintf("unknown output format: '%s'", format)))
  )
  
}



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# BIOM Classic - Tabular
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

write.biom.tsv <- function (biom, file) {
  
  mtx <- {
    rbind(
      matrix(
        data = c("#SampleID", colnames(biom$counts)),
        nrow = 1
      ),
      cbind(
        matrix(
          data = rownames(biom$counts),
          ncol = 1
        ),
        matrix(
          data = as.character(as.matrix(biom$counts)), 
          nrow = nrow(biom$counts)
        )
      )
    )
  }
  
  if (ncol(biom$taxonomy) > 0) {
    mtx <- {
      cbind(
        mtx,
        matrix(
          data = c("Taxonomy", apply(biom$taxonomy, 1L, paste, collapse="; ")),
          ncol = 1
        )
      )
    }
  }
  
  write.table(mtx, file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  return (invisible(NULL))
}



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# BIOM v1.0 - JSON
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


write.biom.1.0 <- function (biom, file) {
  
  for (i in names(biom$metadata))
    if (is(biom$metadata[[i]], "factor"))
      biom$metadata[[i]] <- as.character(biom$metadata[[i]])
  
  
  json <- rjson::toJSON(list(
    
    
    # Attributes
    #------------------------------------------------------
    id                  = biom$info$id,
    type                = biom$info$type,
    format              = "1.0.0", 
    format_url          = "http://biom-format.org",
    generated_by        = paste("rbiom", utils::packageDescription('rbiom')$Version),
    date                = strftime(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC"),
    matrix_type         = "sparse",
    matrix_element_type = ifelse(sum(biom$counts$v %% 1) == 0, "int", "float"),
    shape               = c(biom$counts$nrow, biom$counts$ncol),
    
    
    # Phylogenetic Tree
    #------------------------------------------------------
    phylogeny = ifelse(is.null(biom$phylogeny), "", ape::write.tree(biom$phylogeny)),
    
    
    # Taxonomy
    #------------------------------------------------------
    rows = lapply(1:biom$counts$nrow, function (i) list(
      id = biom$counts$dimnames[[1]][i],
      metadata = list(
        taxonomy = unname(biom$taxonomy[i,])
      )
    )),
    
    
    # Sample IDs and Metadata
    #------------------------------------------------------
    columns = lapply(1:biom$counts$ncol, function (j) list(
      id = biom$counts$dimnames[[2]][j],
      metadata = {
        md_fields <- names(biom$metadata)
        if (length(md_fields) > 0) {
          sapply(md_fields, function (k) biom$metadata[j,k])
        } else {
          NULL
        }
      }
    )),
    
    
    # Read counts
    #------------------------------------------------------
    data = lapply(seq_along(biom$counts$v), function (k) 
      c(biom$counts$i[k] - 1, biom$counts$j[k] - 1, biom$counts$v[k])
    )
  ))
  
  
  res <- try(writeChar(json, file, eos=NULL), silent = TRUE)
  if (is(res, "try-error"))
    stop(simpleError(sprintf("Can't save to '%s': %s", file, res)))
  
  return (invisible(NULL))
}




#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# BIOM v2.1 - HDF5
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

write.biom.2.1 <- function (biom, file) {
  
  hdf5 <- try(h5::h5file(name = file, mode = "w"), silent = TRUE)
  if (is(hdf5, "try-error"))
    stop(simpleError(sprintf("Can't save to '%s': %s", file, hdf5)))
  
  
  # Attributes
  #------------------------------------------------------
  h5::h5attr(hdf5, 'id')             <- biom$info$id
  h5::h5attr(hdf5, 'type')           <- biom$info$type
  h5::h5attr(hdf5, 'format-url')     <- "http://biom-format.org"
  h5::h5attr(hdf5, 'format-version') <- c(2,1,0)
  h5::h5attr(hdf5, 'generated-by')   <- paste("rbiom", utils::packageDescription('rbiom')$Version)
  h5::h5attr(hdf5, 'creation-date')  <- strftime(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC")
  h5::h5attr(hdf5, 'shape')          <- c(biom$counts$nrow, biom$counts$ncol)
  h5::h5attr(hdf5, 'nnz')            <- length(biom$counts$v)
  
  
  
  # Read counts by taxa (rows)
  #------------------------------------------------------
  x <- matrix(c(biom$counts$i - 1, biom$counts$j - 1, biom$counts$v), byrow=FALSE, ncol=3)
  x <- x[order(x[,1]),,drop=FALSE]
  hdf5['observation/ids']            <- biom$counts$dimnames[[1]]
  hdf5['observation/matrix/data']    <- x[,3]
  hdf5['observation/matrix/indices'] <- x[,2]
  hdf5['observation/matrix/indptr']  <- cumsum(unname(table(factor(x[,1]+1, 0:biom$counts$nrow))))
  
  
  
  # Read counts by sample (columns)
  #------------------------------------------------------
  x <- x[order(x[,2]),,drop=FALSE]
  hdf5['sample/ids']                 <- biom$counts$dimnames[[2]]
  hdf5['sample/matrix/data']         <- x[,3]
  hdf5['sample/matrix/indices']      <- x[,1]
  hdf5['sample/matrix/indptr']       <- cumsum(unname(table(factor(x[,2]+1, 0:biom$counts$ncol))))
  
  
  
  # Sample Metadata
  #------------------------------------------------------
  if (is(biom[['metadata']], "data.frame")) {
    plyr::l_ply(setdiff(names(biom$metadata), 'SampleID'), function (field) {
      h5path       <- sprintf("/sample/metadata/%s", field)
      values       <- biom$metadata[biom$counts$dimnames[[2]], field]
      hdf5[h5path] <- if (is.factor(values)) as.character(values) else values
    })
  }
  
  
  
  # Taxonomy
  #------------------------------------------------------
  if (is(biom[['taxonomy']], "matrix")) {
    x           <- biom$taxonomy[biom$counts$dimnames[[1]],,drop=FALSE]
    dimnames(x) <- list(NULL, NULL)
    hdf5['observation/metadata/taxonomy'] <- x
  }
  
  
  
  # Phylogenetic Tree
  #------------------------------------------------------
  if (is(biom[['phylogeny']], "phylo")) {
    h5path       <- '/observation/group-metadata/phylogeny'
    hdf5[h5path] <- ape::write.tree(biom[['phylogeny']])
    h5::h5attr(hdf5[h5path], 'data_type') <- "newick"
  }
  
  
  h5::h5close(hdf5)
  
  return (invisible(NULL))
}
