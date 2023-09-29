#' Get or set the sample names.
#' 
#' @family setters
#' 
#' @param biom,x   A \code{BIOM} object, as returned from [read_biom()].
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
  stopifnot(is(biom, 'BIOM'))
  return (colnames(biom[['counts']]))
}


#' @rdname sample_names
#' @export

`sample_names<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is_character(value))
  stopifnot(all(!is.na(value)))
  stopifnot(all(nchar(value) > 0))
  stopifnot(all(!duplicated(value)))
  stopifnot(length(value) == n_samples(x))
  
  value <- unname(value)
  
  rownames(x[['metadata']]) <- value
  colnames(x[['counts']])   <- value
  
  return (x)
}



#' Get or set the OTU names.
#' 
#' @family setters
#' 
#' @param biom,x   A \code{BIOM} object, as returned from [read_biom()].
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
  stopifnot(is(biom, 'BIOM'))
  return (rownames(biom[['counts']]))
}


#' @rdname otu_names
#' @export

`otu_names<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is_character(value))
  stopifnot(all(!is.na(value)))
  stopifnot(all(nchar(value) > 0))
  stopifnot(all(!duplicated(value)))
  stopifnot(length(value) == n_otus(x))
  
  value <- unname(value)
  
  if (has_tree(x)) { # Likely ordered differently
    idx <- match(otu_names(x), x[['phylogeny']][['tip.label']])
    x[['phylogeny']][['tip.label']][idx] <- value
  }
  
  if (has_sequences(x))
    names(x[['sequences']]) <- value
  
  rownames(x[['counts']])   <- value
  rownames(x[['taxonomy']]) <- value
  
  return (x)
}



#' Get or set the names of the taxonomic ranks.
#' 
#' @family setters
#' 
#' @param biom,x  A \code{BIOM} object, as returned from [read_biom()].
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
#'     taxa_ranks(biom) <- paste0("Level", seq_len(n_ranks(biom)))
#'     taxa_ranks(biom)

taxa_ranks <- function (biom) {
  stopifnot(is(biom, 'BIOM'))
  return (colnames(biom[['taxonomy']]))
}


#' @rdname taxa_ranks
#' @export

`taxa_ranks<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is_character(value))
  stopifnot(all(!is.na(value)))
  stopifnot(all(nchar(value) > 0))
  stopifnot(all(!duplicated(value)))
  stopifnot(length(value) == n_ranks(x))
  
  value <- unname(value)
  
  colnames(x[['taxonomy']]) <- value
  
  return (x)
}



#' Get or set the abundance counts.
#'
#' @family setters
#' 
#' @param biom,x  A \code{BIOM} object, as returned from [read_biom()].
#' 
#' @param value A numeric matrix. Rownames and colnames must be identical to
#'        the current [otu_matrix()] value.
#'        
#' @return A numeric matrix with samples as row names, and OTU identifiers as 
#'         column names.
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

otu_matrix <- function (biom) {
  stopifnot(is(biom, 'BIOM'))
  return (as.matrix(biom[['counts']]))
}


#' @rdname otu_matrix
#' @export
`otu_matrix<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.matrix(value))
  stopifnot(is.numeric(value))
  
  stopifnot(identical(dimnames(otu_matrix(x)), dimnames(value)))
  
  x[['counts']] <- slam::as.simple_triplet_matrix(value)
  
  return (x)
}




#' Get or set the phylogenetic tree.
#' 
#' @family setters
#' 
#' @param biom,x  A \code{BIOM} object, as returned from [read_biom()].
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
  stopifnot(is(biom, 'BIOM'))
  return (biom[['phylogeny']])
}


#' @rdname otu_tree
#' @export

`otu_tree<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  
  if (is_scalar_character(value))
    value <- read_tree(value)
  
  stopifnot(is(value, 'phylo'))
  stopifnot(all(otu_names(x) %in% value$tip.label))
  
  if (length(value$tip.labels) > length(otu_names(x)))
    value <- tree_subset(value, otu_names(x))
  
  x[['phylogeny']] <- value
  
  return (x)
}



#' Get or set the nucleotide sequences associated with each OTU.
#' 
#' @family setters
#' 
#' @param biom,x  A \code{BIOM} object, as returned from [read_biom()].
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
  stopifnot(is(biom, 'BIOM'))
  return (biom[['sequences']])
}


#' @rdname otu_sequences
#' @export

`otu_sequences<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is_null(value) || is_character(value))
  
  if (length(value) == 1 && is.null(names(value)))
    value <- read_fasta(value, ids = otu_names(x))
  
  stopifnot(all(otu_names(x) %in% names(value)))
  
  x[['sequences']] <- value[otu_names(x)]
  
  return (x)
}



#' Get or set a \code{BIOM} object's id or comment.
#' 
#' The BIOM specification includes \code{id} and \code{comment} fields
#' for free-form text.
#' 
#' @family setters
#' 
#' @param biom,x  A \code{BIOM} object, as returned from [read_biom()].
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
  stopifnot(is(biom, 'BIOM'))
  return (biom[['info']][['id']])
}


#' @rdname biom_id
#' @export

`biom_id<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is_scalar_character(value) && !is_na(value))
  
  value <- trimws(value, whitespace = "[\\h\\v]")
  
  x[['info']][['id']] <- value
  
  return (x)
}


#' @rdname biom_id
#' @family setters
#' @export
biom_comment <- function (biom) {
  stopifnot(is(biom, 'BIOM'))
  return (biom[['info']][['comment']])
}


#' @rdname biom_id
#' @export

`biom_comment<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is_scalar_character(value) && !is_na(value))
  
  value <- trimws(value, whitespace = "[\\h\\v]")
  
  x[['info']][['comment']] <- value
  
  return (x)
}


