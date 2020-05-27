#' Generate a matrix of samples by taxa, at the specified taxonomic rank.
#'
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'     object, as returned from \link{read.biom}. For matrices, the rows and 
#'     columns are assumed to be the taxa and samples, respectively.
#' @param rank  The taxonomic rank. E.g. \bold{\dQuote{OTU}}, 
#'     \bold{\dQuote{Phylum}}, etc. May also be given numerically: 0 for
#'     OTU, 1 for the highest level (i.e. Kingdom), and extending to the number
#'     of taxonomic ranks encoded in the original biom file. See example below
#'     to fetch the names of all available ranks.
#' @param map  A character matrix defining the value that each taxa IDs is
#'     assigned for each taxonomic rank. If \code{map=NULL} and \code{biom} is a
#'     \code{BIOM} class object, the map will be automatically loaded from 
#'     \code{biom$taxonomy}. \code{map} must not be \code{null} when \code{biom}
#'     is a \code{matrix} or \code{simple_triplet_matrix}. See the example below
#'     for an example of \code{map}'s structure.
#' @param lineage  Include all ranks in the name of the taxa. For instance,
#'     setting to \code{TRUE} will produce 
#'     \code{Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales}. 
#'     Whereas setting to \code{FALSE} (the default) will return simply
#'     \code{Coriobacteriales}. You want to set this to TRUE if you have
#'     genus names (such as \emph{Incertae_Sedis}) that map to multiple higher
#'     level ranks.
#' @param sparse  If true, returns a sparse matrix as described by 
#'     \code{slam::simple_triplet_matrix}, otherwise returns a normal R
#'     matrix object. Sparse matrices will likely be considerably more
#'     memory efficient in this scenario.
#' @return A numeric matrix with samples as column names, and taxonomic
#'     identifiers as row names.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     colnames(biom$taxonomy)
#'     
#'     phyla <- taxa.rollup(biom, 'Phylum')
#'     phyla[1:4,1:6]
#'     
#'     # Custom matrices should be formatted like so:
#'     counts <- as.matrix(biom$counts)
#'     map    <- biom$taxonomy
#'     
#'     counts[1:3,1:6]
#'     map[1:3,1:4]
#'     
#'     phyla <- taxa.rollup(counts, 'Phylum', map=map)
#'     phyla[1:3,1:6]
#'


taxa.rollup <- function (biom, rank='OTU', map=NULL, lineage=FALSE, sparse=FALSE) {
  
  
  #--------------------------------------------------------------
  # Get the input into a simple_triplet_matrix
  #--------------------------------------------------------------
  
  if (is(biom, "simple_triplet_matrix")) { counts <- biom
  } else if (is(biom, "BIOM"))           { counts <- biom$counts
  } else if (is(biom, "matrix"))         { counts <- slam::as.simple_triplet_matrix(biom)
  } else {
    stop(simpleError("biom must be a matrix, simple_triplet_matrix, or BIOM object."))
  }
  
  if (is.null(map) & is(biom, "BIOM")) map <- biom$taxonomy
  
  
  #--------------------------------------------------------------
  # Sanity Checks
  #--------------------------------------------------------------
  
  if (!is(biom, "BIOM") & is.null(map))
    stop(simpleError("Please provide a taxonomy map."))
  
  if (!is(map, "matrix") | !identical(typeof(map), "character"))
    stop(simpleError("Taxonomy map must be a character matrix."))
  
  if (!all(rownames(counts) %in% rownames(map)))
    stop(simpleError("Taxonomy map does not include all taxa IDs."))
  
  ranks <- colnames(map)
  if (is(rank, "character"))
    rank <- grep(rank, c('OTU', ranks), ignore.case=TRUE) - 1
  
  if (length(rank) != 1 | !is(rank, "numeric"))
      stop(simpleError("Invalid taxonomic rank specified."))
  
  if (rank > ncol(map))
      stop(simpleError("Invalid taxonomic rank specified."))
  
  
  #--------------------------------------------------------------
  # The OTU table is easy
  #--------------------------------------------------------------
  
  if (rank == 0) {
    if (identical(sparse, TRUE)) {
      return (biom[['counts']])
    } else {
      return (as.matrix(biom[['counts']]))
    }
  }
  
  
  #--------------------------------------------------------------
  # compute abundance matrix
  #--------------------------------------------------------------
  
  if (identical(lineage, TRUE)) rank <- 1:rank
  
  res <- try(silent=TRUE, {
    otuNames <- counts$dimnames[[1]]
    otu2taxa <- map[otuNames, rank, drop=FALSE]
    otu2taxa <- apply(otu2taxa, 1L, paste, collapse="; ")
    slam::rollup(counts, 1L, otu2taxa, sum)
  })
  if (!is(res, "simple_triplet_matrix"))
    stop(simpleError("Unable to group by taxonomic level."))
  
  res <- res[order(tolower(rownames(res))), colnames(counts), drop=FALSE]
  
  
  #--------------------------------------------------------------
  # Convert to a regular matrix if requested via sparse=FALSE
  #--------------------------------------------------------------
  
  if (identical(sparse, FALSE))
    res <- as.matrix(res)
  
  
  
  return (res)
}
