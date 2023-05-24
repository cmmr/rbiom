#' Set the sample names.
#' 
#' @param x   A \code{BIOM} object, as returned from \link{read_biom}.
#' @param value   A character vector. A named character vector can be used to
#'        only change some of the sample names.
#' @family setters
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     
#'     head(sample_names(biom))
#'     
#'     sample_names(biom) <- sub('HMP', 'Sample_', sample_names(biom))
#'     head(sample_names(biom))
#'     
#'     sample_names(biom) <- c('Sample_02' = 'One', 'Sample_03' = 'Two')
#'     head(sample_names(biom))
#'

`sample_names<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.character(value))
  
  if (is_null(names(value))) {
    new_ids <- value
    
  } else {
    stopifnot(!any(duplicated(names(value))))
    stopifnot(all(names(value) %in% sample_names(x)))
    
    new_ids <- sample_names(x)
    new_ids[names(value)] <- unname(value)
  }
  
  stopifnot(!any(duplicated(new_ids)))
  
  rownames(x[['metadata']]) <- new_ids
  colnames(x[['counts']])   <- new_ids
  
  return (x)
}


#' Set the taxa names.
#' 
#' @param x   A \code{BIOM} object, as returned from \link{read_biom}.
#' @param value   A character vector. A named character vector can be used to
#'        only change some of the taxa names.
#' @family setters
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     
#'     head(taxa_names(biom))
#'     
#'     taxa_names(biom) <- sub('Unc', 'Uncultured_', taxa_names(biom))
#'     head(taxa_names(biom))
#'     
#'     taxa_names(biom) <- c('AnmMass2' = 'One', 'PreBivi6' = 'Two')
#'     head(taxa_names(biom))
#'

`taxa_names<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.character(value))
  
  if (is_null(names(value))) {
    new_ids <- value
    
  } else {
    stopifnot(!any(duplicated(names(value))))
    stopifnot(all(names(value) %in% taxa_names(x)))
    
    new_ids <- taxa_names(x)
    new_ids[names(value)] <- unname(value)
  }
  
  stopifnot(!any(duplicated(new_ids)))
  
  if (has_phylogeny(x)) { # Likely ordered differently
    idx <- match(taxa_names(x), x[['phylogeny']][['tip.label']])
    x[['phylogeny']][['tip.label']][idx] <- new_ids
  }
  
  if (has_sequences(x))
    names(x[['sequences']]) <- new_ids
  
  rownames(x[['counts']])   <- new_ids
  rownames(x[['taxonomy']]) <- new_ids
  
  return (x)
}


#' Set the taxa rank names.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read_biom}.
#' @param value  A character vector. A named character vector can be used to
#'        only change some of the taxa rank names.
#' @family setters
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     
#'     taxa_ranks(biom)
#'     
#'     taxa_ranks(biom) <- c("OTU" = "ASV")
#'     taxa_ranks(biom)
#'     
#'     taxa_ranks(biom) <- paste0("Level", 1:7)
#'     taxa_ranks(biom)
#'

`taxa_ranks<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.character(value))
  
  if (is_null(names(value))) {
    new_ids <- value
    
  } else {
    stopifnot(!any(duplicated(names(value))))
    stopifnot(all(names(value) %in% taxa_ranks(x)))
    
    new_ids <- taxa_ranks(x)
    new_ids[names(value)] <- unname(value)
  }
  
  stopifnot(!any(duplicated(new_ids)))
  
  colnames(x[['taxonomy']]) <- new_ids
  
  return (x)
}


#' Set the abundance counts.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read_biom}.
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
#' @param x  A \code{BIOM} object, as returned from \link{read_biom}.
#' @param value A \code{phylo} class object with tip.labels matching 
#'        \code{taxa_names(x)}. If there are more tip.labels than
#'        taxa names, then the tree will be subset.
#' @family setters
#' @export
#' @examples
#' \dontrun{
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     phylogeny(biom) <- read_tree('path/to/newick.tre')
#' }
#' 

`phylogeny<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  
  if (is_scalar_character(value))
    value <- read_tree(value)
  
  stopifnot(is(value, 'phylo'))
  stopifnot(all(taxa_names(x) %in% value$tip.label))
  
  if (length(value$tip.labels) > length(taxa_names(x)))
    value <- subtree(value, taxa_names(x))
  
  x[['phylogeny']] <- value
  
  return (x)
}


#' Set nucleotide sequences associated with each taxonomic identifier.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read_biom}.
#' @param value A named character vector. Names must match \code{taxa_names(x)}.
#' @family setters
#' @export
#' @examples
#' \dontrun{
#'     library(rbiom)
#'     
#'     sequences(biom) <- read_fasta('path/to/sequences.fa')
#'     sequences(biom) <- c(OTU1 = 'ATCGGGTA', OTU2 = 'GGCATTAGC')
#' }
#' 

`sequences<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.character(value))
  
  if (length(value) == 1 && is.null(names(value)))
    value <- read_fasta(value, ids = taxa_names(x))
  
  stopifnot(all(taxa_names(x) %in% names(value)))
  
  x[['sequences']] <- value[taxa_names(x)]
  
  return (x)
}


#' Set a \code{BIOM} object's id.
#' 
#' The BIOM specification includes \code{id} and \code{comment} fields
#' for free-form text.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read_biom}.
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

`id<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is_scalar_character(value))
  stopifnot(!is.na(value))
  
  value <- trimws(value, whitespace = "[\\h\\v]")
  
  x[['info']][['id']] <- value
  
  return (x)
}


#' Set a \code{BIOM} object's comment.
#' 
#' The BIOM specification includes \code{id} and \code{comment} fields
#' for free-form text.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read_biom}.
#' @param value The identifier to add (character vector of length 1).
#' @family setters
#' @export
#' @examples
#'     library(rbiom)
#'     biom <- hmp50
#'     
#'     comments(biom) <- "A description of this study"
#'     comments(biom)
#'

`comments<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is_scalar_character(value))
  stopifnot(!is.na(value))
  
  value <- trimws(value, whitespace = "[\\h\\v]")
  
  x[['info']][['comment']] <- value
  
  return (x)
}
