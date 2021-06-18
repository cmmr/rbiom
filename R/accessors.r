#' Get the sample names.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A character vector of the sample IDs / names in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     sample.names(biom)
#'

sample.names <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In sample.names(), biom must be a BIOM-class object.'))
  return (colnames(biom[['counts']]))
}


#' Sum the observations in each sample.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A named numeric vector of the number of observations in each 
#'         sample. The names are the sample IDs.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     sample.sums(biom)
#'     
#'     sample.sums(biom) %>% sort() %>% head()
#'     
#'     hist(sample.sums(biom))
#'

sample.sums <- function (biom, long=FALSE, md=FALSE) {
  
  if (!is(biom, 'BIOM'))
    stop (simpleError('In sample.names(), biom must be a BIOM-class object.'))
  
  result <- slam::col_sums(biom[['counts']])
  
  if (isTRUE(long) || !isFALSE(md)) {
    
    #--------------------------------------------------------------
    # Convert to long format
    #--------------------------------------------------------------
    result <- data.frame(
      stringsAsFactors = FALSE,
      'Sample' = names(result),
      'Reads'  = unname(result)
      )
    
    #--------------------------------------------------------------
    # Add Metadata
    #--------------------------------------------------------------
    if (identical(md, TRUE))  md <- colnames(rbiom::metadata(biom))
    if (identical(md, FALSE)) md <- c()
    for (i in unique(md))
      result[[i]] <- metadata(biom, i)[result[['Sample']]]
    
  }
  
  return (result)
  
}


#' Get the taxa names.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A character vector of the taxa IDs / names in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     taxa.names(biom) %>% head()
#'

taxa.names <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In taxa.names(), biom must be a BIOM-class object.'))
  return (rownames(biom[['counts']]))
}


#' Get the taxa ranks.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A character vector of the taxa ranks in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     taxa.ranks(biom)
#'

taxa.ranks <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In taxa.ranks(), biom must be a BIOM-class object.'))
  return (colnames(biom[['taxonomy']]))
}


#' Get the abundance counts.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A numeric matrix of the sample abundance counts in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     counts(biom)[1:4,1:5]
#'

counts <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In counts(), biom must be a BIOM-class object.'))
  return (as.matrix(biom[['counts']]))
}


#' Checks if a BIOM object is rarefied.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return \code{TRUE} if the BIOM object is rarefied, \code{FALSE} otherwise.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     is.rarefied(biom)
#'

is.rarefied <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In is.rarefied(), biom must be a BIOM-class object.'))
  return (length(unique(sample.sums(biom))) == 1)
}


#' The rarefaction depth of a BIOM object.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return The rarefaction depth. If the BIOM object is not rarefied, this will be 
#'         a sorted vector of all the unique depths.
#' @family accessor functions
#' @seealso [sample.sums()]
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     depth(biom) %>% head()
#'     
#'     rarefy(biom) %>% depth()
#'

depth <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In depth(), biom must be a BIOM-class object.'))
  return (sort(unique(sample.sums(biom))))
}


#' Checks if a phylogenetic tree is present.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return \code{TRUE} if a phylogenetic tree is present, \code{FALSE} otherwise.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     has.phylogeny(biom)
#'

has.phylogeny <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In has.phylogeny(), biom must be a BIOM-class object.'))
  return (is(biom[['phylogeny']], 'phylo'))
}


#' Get the phylogenetic tree.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A \code{phylo} class object of the tree in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     summary(phylogeny(biom))
#'

phylogeny <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In phylogeny(), biom must be a BIOM-class object.'))
  return (biom[['phylogeny']])
}


#' Get the sample metadata.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @param field  The name of the metadata field to retrieve (optional).
#' @return A data frame of the metadata in \code{biom}. If \code{field} is
#'    given, will return a named vector of the field's values.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     metadata(biom)[1:4,1:3]
#'

metadata <- function (biom, field=NULL) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In metadata(), biom must be a BIOM-class object.'))
  
  if (is.null(field))
    return (biom[['metadata']])
  
  if (!field %in% names(biom[['metadata']]))
    stop(paste0("Field '", field, "' is not present in the metadata."))
  
  return (setNames(biom[['metadata']][[field]], rownames(biom[['metadata']])))
}


#' Checks if a DNA sequences are present.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return \code{TRUE} if DNA sequences are present, \code{FALSE} otherwise.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     has.sequences(biom)
#'

has.sequences <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In has.sequences(), biom must be a BIOM-class object.'))
  return (!is.null(biom[['sequences']]))
}


#' DNA sequence associated with each taxonomic identifier.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A named character vector of sequences in \code{biom}. If this data
#'           is not present, then returns \code{NULL}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     sequences(biom)[1:4]
#'     
#'     # Write to a compressed fasta file in the temporary directory:
#'     seqs <- sequences(biom)
#'     conn <- bzfile(file.path(tempdir(), "Sequences.fa.bz2"), "w")
#'     cat(sprintf(">%s\n%s", names(seqs), seqs), file=conn, sep="\n")
#'     close(conn)
#'     
#'     # You can also use the write.fasta function for this task:
#'     write.fasta(biom, file.path(tempdir(), "Sequences.fa.gz"))
#'

sequences <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In sequences(), biom must be a BIOM-class object.'))
  return (biom[['sequences']])
}


#' Get biom's misc information.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A data frame of the metadata in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     info(biom)
#'

info <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In info(), biom must be a BIOM-class object.'))
  return (biom[['info']])
}


#' Number of samples in a BIOM.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return The number of samples present.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     nsamples(biom)
#'

nsamples <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In nsamples(), biom must be a BIOM-class object.'))
  return (ncol(biom[['counts']]))
}


#' Number of taxa in a BIOM.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return The number of taxa present.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     ntaxa(biom)
#'

ntaxa <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In ntaxa(), biom must be a BIOM-class object.'))
  return (nrow(biom[['counts']]))
}
