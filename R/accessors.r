#' Get the sample names.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A character vector of the sample IDs / names in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     sample.names(hmp50)
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
#'     sample.sums(hmp50) %>% sort() %>% head()
#'     
#'     hist(sample.sums(hmp50))
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
#'     taxa.names(hmp50) %>% head()
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
#'     taxa.ranks(hmp50)
#'

taxa.ranks <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In taxa.ranks(), biom must be a BIOM-class object.'))
  return (c(colnames(biom[['taxonomy']]), 'OTU'))
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
#'     counts(hmp50)[1:4,1:5]
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
#'     is.rarefied(hmp50)
#'     
#'     rarefy(hmp50, 1000) %>% is.rarefied()
#'

is.rarefied <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In is.rarefied(), biom must be a BIOM-class object.'))
  return (isTRUE(attr(biom, 'rarefaction') > 0))
}


#' The rarefaction depth of a BIOM object.
#' 
#' The returned values reflect changes made after any \code{subset()} commands. 
#' To get the rarefaction level immediately after \code{read.biom()} or
#' \code{rarefy()} use \code{attr(biom, 'rarefaction')} instead.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return The rarefaction depth. If the BIOM object is not rarefied, this will 
#'         be a sorted vector of all the unique depths.
#' @family accessor functions
#' @seealso [sample.sums()]
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     depth(hmp50) %>% head()
#'     rarefy(hmp50) %>% depth()
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
#'     has.phylogeny(hmp50)
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
#'     summary(phylogeny(hmp50))
#'

phylogeny <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In phylogeny(), biom must be a BIOM-class object.'))
  return (biom[['phylogeny']])
}


#' Checks if DNA sequences are present.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return \code{TRUE} if DNA sequences are present, \code{FALSE} otherwise.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     has.sequences(hmp50)
#'

has.sequences <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In has.sequences(), biom must be a BIOM-class object.'))
  return (!is.null(biom[['sequences']]))
}


#' Nucleotide sequences associated with each taxonomic identifier.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A named character vector of sequences in \code{biom}. If this data
#'           is not present, then returns \code{NULL}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     sequences(hmp50)[1:4]
#'     
#'     # Write to a compressed fasta file in the temporary directory:
#'     seqs <- sequences(hmp50)
#'     conn <- bzfile(file.path(tempdir(), "Sequences.fa.bz2"), "w")
#'     cat(sprintf(">%s\n%s", names(seqs), seqs), file=conn, sep="\n")
#'     close(conn)
#'     
#'     # You can also use the write.fasta function for this task:
#'     write.fasta(hmp50, file.path(tempdir(), "Sequences.fa.gz"))
#'

sequences <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In sequences(), biom must be a BIOM-class object.'))
  return (biom[['sequences']])
}


#' Get \code{BIOM} object's miscellaneous information.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A list of the top-level metadata in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     info(hmp50)
#'

info <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In info(), biom must be a BIOM-class object.'))
  return (c(biom[['info']], comment(biom)))
}


#' Get \code{BIOM} object's identifier / title.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A data frame of the metadata in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     id(hmp50)
#'

id <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In id(), biom must be a BIOM-class object.'))
  return (biom[['info']][['id']])
}


#' Number of samples in a BIOM.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return The number of samples present.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     nsamples(hmp50)
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
#'     ntaxa(hmp50)
#'

ntaxa <- function (biom) {
  if (!is(biom, 'BIOM'))
    stop (simpleError('In ntaxa(), biom must be a BIOM-class object.'))
  return (nrow(biom[['counts']]))
}
