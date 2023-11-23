
#' Get or set the sample names.
#' 
#' @inherit documentation_default
#' 
#' @family setters
#' 
#' @param value   A character vector of the new sample names.
#' 
#' @return A character vector of the sample names in \code{biom}.
#' 
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

sample_names <- function (biom) {
  validate_biom(clone = FALSE)
  return (colnames(biom[['counts']]))
}


#' @rdname sample_names
#' @export

`sample_names<-` <- function (biom, value) {
  
  validate_biom(clone = TRUE)
  
  lvls <- levels(value)
  if (is.factor(value)) {
    value %<>% as.character()
    lvls  %<>% intersect(value)
  }
  
  stopifnot(is_character(value))
  stopifnot(all(!is.na(value)))
  stopifnot(all(nchar(value) > 0))
  stopifnot(all(!duplicated(value)))
  stopifnot(length(value) == n_samples(biom))
  
  value <- unname(value)
  
  colnames(biom[['counts']]) <- value
  
  if (!is.null(lvls)) { value %<>% factor(levels = lvls)
  } else              { value %<>% as.factor() }
  
  biom[['metadata']][['.sample']] <- value
  
  invalidate_biom()
  return (biom)
}



#' Get or set the OTU names.
#' 
#' @inherit documentation_default
#' 
#' @family setters
#' 
#' @param value   A character vector with the new OTU names.
#' 
#' @return A character vector of the taxa names in \code{biom}.
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     biom <- hmp50
#'     
#'     head(otu_names(biom))
#'     
#'     otu_names(biom) <- sub('Unc', 'Uncultured_', otu_names(biom))
#'     head(otu_names(biom))

otu_names <- function (biom) {
  validate_biom(clone = FALSE)
  return (rownames(biom[['counts']]))
}


#' @rdname otu_names
#' @export

`otu_names<-` <- function (biom, value) {
  
  validate_biom(clone = TRUE)
  
  lvls <- levels(value)
  if (is.factor(value)) {
    value %<>% as.character()
    lvls  %<>% intersect(value)
  }
  
  stopifnot(is_character(value))
  stopifnot(!is.na(value))
  stopifnot(nchar(value) > 0)
  stopifnot(!duplicated(value))
  stopifnot(length(value) == n_otus(biom))
  
  value <- unname(value)
  
  if (has_tree(biom)) { # Likely ordered differently
    idx <- match(otu_names(biom), biom[['phylogeny']][['tip.label']])
    biom[['phylogeny']][['tip.label']][idx] <- value
  }
  
  if (has_sequences(biom))
    names(biom[['sequences']]) <- value
  
  rownames(biom[['counts']]) <- value
  
  
  if (!is.null(lvls)) { value %<>% factor(levels = lvls)
  } else              { value %<>% as.factor() }
  
  biom[['taxonomy']][['.otu']] <- value
  
  invalidate_biom()
  return (biom)
}



#' Get or set the names of the taxonomic ranks.
#' 
#' @inherit documentation_default
#' 
#' @family setters
#' 
#' @param value  A character vector the new taxa rank names.
#' 
#' @return A character vector of the taxa rank names in \code{biom}.
#'        
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     biom <- hmp50
#'     
#'     taxa_ranks(biom)
#'     
#'     taxa_ranks(biom) <- c('.otu', paste0("Level", seq_len(n_ranks(biom) - 1)))
#'     taxa_ranks(biom)

taxa_ranks <- function (biom) {
  validate_biom(clone = FALSE)
  res <- colnames(biom[['taxonomy']])
  return (res)
}


#' @rdname taxa_ranks
#' @export

`taxa_ranks<-` <- function (biom, value) {
  
  validate_biom(clone = TRUE)
  
  stopifnot(is_character(value))
  stopifnot(!is.na(value))
  stopifnot(nchar(value) > 0)
  stopifnot(!duplicated(value))
  stopifnot(length(value) == n_ranks(biom))
  stopifnot(eq(value[[1]], '.otu'))
  stopifnot(!startsWith(value, '.') | value == '.otu')
  
  colnames(biom[['taxonomy']]) <- unname(value)
  
  invalidate_biom()
  return (biom)
}



#' Get or set the metadata column names.
#' 
#' @inherit documentation_default
#' 
#' @family setters
#' 
#' @param value  A character vector of the new metadata column names. Names
#'        cannot start with a \code{.} character. Leading and trailing
#'        whitespace is automatically removed.
#' 
#' @return A character vector of the metadata column names in \code{biom}.
#'        
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     biom <- hmp50
#'     
#'     metadata_names(biom)
#'     
#'     metadata_names(biom) <- sub(' ', '_', tolower(metadata_names(biom)))
#'     metadata_names(biom)

metadata_names <- function (biom) {
  validate_biom(clone = FALSE)
  res <- colnames(biom[['metadata']])
  return (res)
}


#' @rdname metadata_names
#' @export

`metadata_names<-` <- function (biom, value) {
  
  validate_biom(clone = TRUE)
  
  stopifnot(is_character(value))
  stopifnot(!anyNA(value))
  
  value <- trimws(as.vector(value))
  stopifnot(nchar(value) > 0)
  stopifnot(!duplicated(value))
  stopifnot(length(value) == n_metadata(biom))
  stopifnot(eq(value[[1]], '.sample'))
  stopifnot(!startsWith(value, '.') | value == '.sample')
  
  colnames(biom[['metadata']]) <- value
  
  return (biom)
}



#' Get or set the abundance counts.
#' 
#' @inherit documentation_default
#'
#' @family setters
#' 
#' @param value A numeric matrix. Rownames and colnames must be eq to
#'        the current [otu_matrix()] value.
#'        
#' @return A numeric matrix with samples as column names, and OTU identifiers 
#'         as row names.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     taxa <- c('Unc53100', 'Unc00c7g', 'Unc25731')
#'     
#'     otu_matrix(biom)[taxa,1:5]
#'     
#'     biom <- sample_rarefy(biom, 100)
#'     otu_matrix(biom)[taxa,1:5]
#'     
#'     otu_matrix(biom) <- otu_matrix(biom) / 100
#'     otu_matrix(biom)[taxa,1:5]

otu_matrix <- function (biom, sparse = FALSE) {
  validate_biom(clone = FALSE)
  if (isTRUE(sparse)) { biom[['counts']]
  } else              { as.matrix(biom[['counts']]) }
}


#' @rdname otu_matrix
#' @export
`otu_matrix<-` <- function (biom, value) {
  
  validate_biom(clone = TRUE)
  
  if (!slam::is.simple_triplet_matrix(value)) {
    stopifnot(is.matrix(value))
    stopifnot(is.numeric(value))
    value <- slam::as.simple_triplet_matrix(value)
  }
  
  biom[['counts']] <- value
  biom %<>% biom_repair()
  
  invalidate_biom()
  return (biom)
}


#' @rdname otu_matrix
#' @keywords internal
#' @export
dimnames.rbiom <- function (x) { dimnames(x$counts) }

#' @rdname otu_matrix
#' @keywords internal
#' @export
dim.rbiom <- function (x) { dim(x$counts) }




#' Get or set the phylogenetic tree.
#' 
#' @inherit documentation_default
#' 
#' @family setters
#' 
#' @param value A \code{phylo} class object with tip.labels matching 
#'        \code{otu_names(x)}. If there are more tip.labels than
#'        taxa names, then the tree will be subset.
#' 
#' @return A \code{phylo} class object of the tree in \code{biom}.
#' 
#' @export
#' @examples
#' \dontrun{
#'     library(rbiom)
#'     
#'     biom <- hmp50
#'     otu_tree(biom) <- read_tree('path/to/newick.tre')
#' }
#' 

otu_tree <- function (biom) {
  validate_biom(clone = FALSE)
  return (biom[['phylogeny']])
}


#' @rdname otu_tree
#' @export

`otu_tree<-` <- function (biom, value) {
  
  validate_biom(clone = TRUE)
  
  if (!is.null(value)) {
    
    if (is_scalar_character(value))
      value <- read_tree(value)
    
    stopifnot(is(value, 'phylo'))
    stopifnot(all(otu_names(biom) %in% value$tip.label))
    
    if (length(value$tip.labels) > length(otu_names(biom)))
      value <- tree_subset(value, otu_names(biom))
  }
  
  biom[['phylogeny']] <- value
  
  invalidate_biom()
  return (biom)
}



#' Get or set the nucleotide sequences associated with each OTU.
#' 
#' @inherit documentation_default
#' 
#' @family setters
#' 
#' @param value A named character vector. Names must match \code{otu_names(x)}.
#' 
#' @return A named character vector of sequences in \code{biom}. If this data
#'           is not present, then returns \code{NULL}.
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     substr(otu_sequences(hmp50)[1:4], 1, 10)
#'     
#' \dontrun{
#'     # Write to a compressed fasta file in the temporary directory:
#'     seqs <- otu_sequences(hmp50)
#'     conn <- bzfile(file.path(tempdir(), "Sequences.fa.bz2"), "w")
#'     cat(sprintf(">%s\n%s", names(seqs), seqs), file=conn, sep="\n")
#'     close(conn)
#'     
#'     # You can also use the write_fasta function for this task:
#'     write_fasta(hmp50, file.path(tempdir(), "Sequences.fa.gz"))
#'     
#'     # Set/replace the sequences
#'     biom <- read_biom('path/to/file.biom')
#'     otu_sequences(biom) <- read_fasta('path/to/sequences.fa')
#'     otu_sequences(biom) <- c(OTU1 = 'ATCGGGTA', OTU2 = 'GGCATTAGC')
#' }
#' 

otu_sequences <- function (biom) {
  validate_biom(clone = FALSE)
  return (biom[['sequences']])
}


#' @rdname otu_sequences
#' @export

`otu_sequences<-` <- function (biom, value) {
  
  validate_biom(clone = TRUE)
  
  stopifnot(is_null(value) || is_character(value))
  
  if (length(value) == 1 && is.null(names(value)))
    value <- read_fasta(value, ids = otu_names(biom))
  
  stopifnot(all(otu_names(biom) %in% names(value)))
  
  biom[['sequences']] <- value[otu_names(biom)]
  
  invalidate_biom()
  return (biom)
}



#' Get or set an \code{rbiom} object's id or comment.
#' 
#' The BIOM specification includes \code{id} and \code{comment} fields
#' for free-form text.
#' 
#' @inherit documentation_default
#' 
#' @family setters
#' 
#' @param value The identifier to add (character vector of length 1).
#' 
#' @return A length 1 character vector.
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     biom <- hmp50
#'     
#'     biom_id(biom)
#'     biom_id(biom) <- "My new title/id"
#'     biom_id(biom)
#'     
#'     biom_comment(biom)
#'     biom_comment(biom) <- "A description of this study"
#'     biom_comment(biom)
#'

biom_id <- function (biom) {
  validate_biom(clone = FALSE)
  return (biom[['info']][['id']])
}


#' @rdname biom_id
#' @export

`biom_id<-` <- function (biom, value) {
  
  validate_biom(clone = TRUE)
  
  stopifnot(is_scalar_character(value) && !is_na(value))
  
  value <- trimws(value, whitespace = "[\\h\\v]")
  
  biom[['info']][['id']] <- value
  
  invalidate_biom()
  return (biom)
}


#' @rdname biom_id
#' @family setters
#' @export
biom_comment <- function (biom) {
  validate_biom(clone = FALSE)
  return (biom[['info']][['comment']])
}


#' @rdname biom_id
#' @export

`biom_comment<-` <- function (biom, value) {
  
  validate_biom(clone = TRUE)
  
  stopifnot(is_scalar_character(value) && !is_na(value))
  
  value <- trimws(value, whitespace = "[\\h\\v]")
  
  biom[['info']][['comment']] <- value
  
  invalidate_biom()
  return (biom)
}


