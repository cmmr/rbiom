#' Set the sample names.
#' 
#' @param x   A \code{BIOM} object, as returned from \link{read.biom}.
#' @param value   A character vector. A named character vector can be used to
#'        only change some of the sample names.
#' @family setters
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     
#'     head(sample.names(biom))
#'     
#'     sample.names(biom) <- sub('HMP', 'Sample_', sample.names(biom))
#'     head(sample.names(biom))
#'     
#'     sample.names(biom) <- c('Sample_02' = 'One', 'Sample_03' = 'Two')
#'     head(sample.names(biom))
#'

`sample.names<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.character(value))
  
  if (is.null(names(value))) {
    new_ids <- value
    
  } else {
    stopifnot(!any(duplicated(names(value))))
    stopifnot(all(names(value) %in% sample.names(x)))
    
    new_ids <- sample.names(x)
    new_ids[names(value)] <- unname(value)
  }
  
  stopifnot(!any(duplicated(new_ids)))
  
  rownames(x[['metadata']]) <- new_ids
  colnames(x[['counts']])   <- new_ids
  
  return (x)
}


#' Set the taxa names.
#' 
#' @param x   A \code{BIOM} object, as returned from \link{read.biom}.
#' @param value   A character vector. A named character vector can be used to
#'        only change some of the taxa names.
#' @family setters
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     
#'     head(taxa.names(biom))
#'     
#'     taxa.names(biom) <- sub('Unc', 'Uncultured_', taxa.names(biom))
#'     head(taxa.names(biom))
#'     
#'     taxa.names(biom) <- c('AnmMass2' = 'One', 'PreBivi6' = 'Two')
#'     head(taxa.names(biom))
#'

`taxa.names<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.character(value))
  
  if (is.null(names(value))) {
    new_ids <- value
    
  } else {
    stopifnot(!any(duplicated(names(value))))
    stopifnot(all(names(value) %in% taxa.names(x)))
    
    new_ids <- taxa.names(x)
    new_ids[names(value)] <- unname(value)
  }
  
  stopifnot(!any(duplicated(new_ids)))
  
  if (has.phylogeny(x)) { # Likely ordered differently
    idx <- match(taxa.names(x), x[['phylogeny']][['tip.label']])
    x[['phylogeny']][['tip.label']][idx] <- new_ids
  }
  
  if (has.sequences(x))
    names(x[['sequences']]) <- new_ids
  
  rownames(x[['counts']])   <- new_ids
  rownames(x[['taxonomy']]) <- new_ids
  
  return (x)
}


#' Set the taxa rank names.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read.biom}.
#' @param value  A character vector. A named character vector can be used to
#'        only change some of the taxa rank names.
#' @family setters
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     
#'     taxa.ranks(biom)
#'     
#'     taxa.ranks(biom) <- c("OTU" = "ASV")
#'     taxa.ranks(biom)
#'     
#'     taxa.ranks(biom) <- paste0("Level", 1:7)
#'     taxa.ranks(biom)
#'

`taxa.ranks<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.character(value))
  
  if (is.null(names(value))) {
    new_ids <- value
    
  } else {
    stopifnot(!any(duplicated(names(value))))
    stopifnot(all(names(value) %in% taxa.ranks(x)))
    
    new_ids <- taxa.ranks(x)
    new_ids[names(value)] <- unname(value)
  }
  
  stopifnot(!any(duplicated(new_ids)))
  
  colnames(x[['taxonomy']]) <- new_ids
  
  return (x)
}


#' Set the abundance counts.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read.biom}.
#' @param value A numeric matrix. Rownames and colnames must be identical to
#'        the current \code{counts()} value.
#' @family setters
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     taxa <- c('Unc53100', 'Unc00c7g', 'Unc25731')
#'     
#'     counts(biom)[taxa,1:5]
#'     
#'     biom <- rarefy(biom, 100)
#'     counts(biom)[taxa,1:5]
#'     
#'     counts(biom) <- counts(biom) / 100
#'     counts(biom)[taxa,1:5]
#'

`counts<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.character(value))
  stopifnot(identical(dimnames(counts(x)), dimnames(value)))
  
  x[['counts']] <- slam::as.simple_triplet_matrix(value)
  
  return (x)
}


#' Set the phylogenetic tree.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read.biom}.
#' @param value A \code{phylo} class object with tip.labels matching 
#'        \code{taxa.names(x)}. If there are more tip.labels than
#'        taxa names, then the tree will be subset.
#' @family setters
#' @export
#' @examples
#' \dontrun{
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     phylogeny(biom) <- read.tree('path/to/newick.tre')
#' }
#' 

`phylogeny<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is(value, 'phylo'))
  stopifnot(all(taxa.names(x) %in% value$tip.labels))
  
  if (length(value$tip.labels) > length(taxa.names(x)))
    value <- subtree(value, taxa.names(x))
  
  x[['phylogeny']] <- value
  
  return (x)
}


#' Set nucleotide sequences associated with each taxonomic identifier.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read.biom}.
#' @param value A named character vector. Names must match \code{taxa.names(x)}.
#' @family setters
#' @export
#' @examples
#' \dontrun{
#'     library(rbiom)
#'     
#'     sequences(biom) <- read.fasta('path/to/sequences.fa')
#'     sequences(biom) <- c(OTU1 = 'ATCGGGTA', OTU2 = 'GGCATTAGC')
#' }
#' 

`sequences<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.character(value))
  stopifnot(all(taxa.names(x) %in% names(value)))
  
  x[['sequences']] <- value[taxa.names(x)]
  
  return (x)
}


#' Set a \code{BIOM} object's id or comment.
#' 
#' The BIOM specification includes \code{id} and \code{comment} fields
#' for free-form text.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read.biom}.
#' @param value The identifier to add (character vector of length 1).
#' @family setters
#' @export
#' @examples
#'     library(rbiom)
#'     biom <- hmp50
#'     
#'     info(biom)$id
#'     
#'     id(biom) <- "My new title/id"
#'     info(biom)$id
#'     
#'     comment(biom) <- "A description of this study"
#'     comment(biom)
#'

`id<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.character(value))
  stopifnot(length(value) == 1)
  
  x[['info']][['id']] <- value
  
  return (x)
}
