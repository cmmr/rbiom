#' Get the sample names.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A character vector of the sample IDs / names in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     sample.names(biom)
#'

sample.names <- function (biom) {
  if (!is(biom, 'BIOM'))
    return (simpleError('In sample.names(), biom must be a BIOM-class object.'))
  return (colnames(biom[['counts']]))
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
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     taxa.names(biom) %>% head()
#'

taxa.names <- function (biom) {
  if (!is(biom, 'BIOM'))
    return (simpleError('In taxa.names(), biom must be a BIOM-class object.'))
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
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     taxa.ranks(biom)
#'

taxa.ranks <- function (biom) {
  if (!is(biom, 'BIOM'))
    return (simpleError('In taxa.ranks(), biom must be a BIOM-class object.'))
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
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     counts(biom)[1:4,1:5]
#'

counts <- function (biom) {
  if (!is(biom, 'BIOM'))
    return (simpleError('In counts(), biom must be a BIOM-class object.'))
  return (as.matrix(biom[['counts']]))
}


#' Get the taxonomy table.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A character matrix of the named taxonomies in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     taxonomy(biom)[1:4,]
#'

taxonomy <- function (biom) {
  if (!is(biom, 'BIOM'))
    return (simpleError('In taxonomy(), biom must be a BIOM-class object.'))
  return (biom[['taxonomy']])
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
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     summary(phylogeny(biom))
#'

phylogeny <- function (biom) {
  if (!is(biom, 'BIOM'))
    return (simpleError('In phylogeny(), biom must be a BIOM-class object.'))
  return (biom[['phylogeny']])
}


#' Get the sample metadata.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A data frame of the metadata in \code{biom}.
#' @family accessor functions
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     metadata(biom)[1:4,1:3]
#'

metadata <- function (biom) {
  if (!is(biom, 'BIOM'))
    return (simpleError('In metadata(), biom must be a BIOM-class object.'))
  return (biom[['metadata']])
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
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     info(biom)
#'

info <- function (biom) {
  if (!is(biom, 'BIOM'))
    return (simpleError('In info(), biom must be a BIOM-class object.'))
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
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     nsamples(biom)
#'

nsamples <- function (biom) {
  if (!is(biom, 'BIOM'))
    return (simpleError('In nsamples(), biom must be a BIOM-class object.'))
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
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     ntaxa(biom)
#'

ntaxa <- function (biom) {
  if (!is(biom, 'BIOM'))
    return (simpleError('In ntaxa(), biom must be a BIOM-class object.'))
  return (nrow(biom[['counts']]))
}
