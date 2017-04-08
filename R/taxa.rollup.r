#' Generate a matrix of samples by taxa, at the specified taxonomic rank.
#'
#' @param biom  A BIOM object, as returned from \link{read.biom}.
#' @param rank  The taxonomic rank. E.g. \bold{\dQuote{OTU}}, 
#'     \bold{\dQuote{Phylum}}, etc. May also be given numerically: 0 for
#'     OTU, 1 for the highest level (i.e. Kingdom), and extending to the number
#'     of taxonomic ranks encoded in the original biom file. See example below
#'     to fetch the names of all available ranks.
#' @param lineage  Include all ranks in the name of the taxa. For instance,
#'     setting to \code{TRUE} will produce 
#'     \code{Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales}. 
#'     Whereas setting to \code{FALSE} (the default) will return simply
#'     \code{Coriobacteriales}. You want to set this to TRUE if you have
#'     genus names (such as \emph{Incertae_Sedis}) that map to multiple higher
#'     level ranks.
#' @param sparse  If true, returns a sparse matrix as described by 
#'     \link[slam]{simple_triplet_matrix}, otherwise returns a normal R
#'     matrix object. Sparse matrices will likely be considerably more
#'     memory efficient in this scenario.
#' @param progressbar  Whether to display a progress bar and status messages
#'     (logical). Will automatically tie in with \pkg{shiny} if run within a
#'     \pkg{shiny} session.
#' @return A numeric matrix with samples as column names, and taxonomic
#'     identifiers as row names.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     colnames(biom$taxonomy)
#'     
#'     phyla <- taxa.rollup(biom, 'Phylum')
#'     phyla[1:4,1:6]
#'


taxa.rollup <- function (biom, rank='OTU', lineage=FALSE, sparse=FALSE, progressbar=FALSE) {
  
  
  #--------------------------------------------------------------
  # Sanity Checks
  #--------------------------------------------------------------
  
  if (!is(biom, "BIOM"))
      stop(simpleError("Invalid BIOM object."))
  
  ranks <- colnames(biom$taxonomy)
  if (is(rank, "character"))
    rank <- grep(rank, c('OTU', ranks), ignore.case=TRUE) - 1
  
  if (length(rank) != 1 | !is(rank, "numeric"))
      stop(simpleError("Invalid taxonomic level."))
  
  if (rank > ncol(biom$taxonomy))
      stop(simpleError("Invalid taxonomic level."))
  
  
  
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
  
  pb <- progressBar(progressbar)
  pb$set(0, paste('Calculating', ranks[rank], 'abundances'))
  
  if (identical(lineage, TRUE)) rank <- 1:rank
  
  res <- try(silent=TRUE, {
    otuNames <- biom$counts$dimnames[[1]]
    otu2taxa <- biom$taxonomy[otuNames, rank, drop=FALSE]
    otu2taxa <- apply(otu2taxa, 1L, paste, collapse="; ")
    slam::rollup(biom$counts, 1L, otu2taxa, sum)
  })
  if (!is(res, "simple_triplet_matrix"))
    stop(simpleError("Unable to group by taxonomic level."))
  
  res <- res[order(tolower(rownames(res))), rownames(biom$metadata), drop=FALSE]
  
  
  #--------------------------------------------------------------
  # Convert it to a regular matrix before returning
  #--------------------------------------------------------------
  
  if (identical(sparse, FALSE))
    res <- as.matrix(res)
  
  
  pb$close()
  
  return (res)
}
